#' Calculate a single episcore
#'
#' This function calculates one epigenetic score using methylation information.
#'
#' @description
#' This function calculates an epigenetic score based on methylation data. Two inputs are required: a beta matrix with CpG methylation and an external reference with the association and effect values for each CpG
#'
#' @import readxl
#' @import dplyr
#' @import data.table
#' @param thrs.criteria Thresholding criteria for CpG selection.
#' @param beta Beta from sample CpG.
#' @param ewas.file Path and name to the EWAS file.
#' @param missingness Maximum tolerated missing data for CpG sites.
#' @return The calculated episcore.
#' @export
#'
calculate_episcore <- function(thrs.criteria = 0.05,
                               beta.file = "test/betamatrix_test.txt",
                               ewas.file = "test/sumstats_test.txt",
                               missingness = 0.2,
                               path = getwd()) {

  start_time <- Sys.time()

  message("╔══════════════════════════════════════════════════════════════╗")
  message("║                  calculate.episcore START                    ║")
  message("╚══════════════════════════════════════════════════════════════╝")
  message("═= Step 1: beta matrix and load summary statistic file")

  episcore_name <- tools::file_path_sans_ext(basename(ewas.file))
  full_path <- file.path(path, ewas.file)

  file_ext <- tools::file_ext(beta.file)
  beta <- if (file_ext == "xlsx") {
    read_xlsx(beta.file, col_names = TRUE)
  } else {
    read.table(beta.file, header = TRUE)
  }
  message("      Reading beta matrix .", file_ext, " file... Done!")


  file_ext <- tools::file_ext(ewas.file)
  ewas_data <- if (file_ext == "xlsx") {
    read_xlsx(ewas.file, col_names = TRUE)
  } else {
    read.table(ewas.file, header = TRUE)
  }
  message("      Reading summary statistics .", file_ext, " file... Done!")
  message("      p value thresholding is set at ", thrs.criteria)

  message("      Maximum p value provided by EWAS summary statistics is ", max(ewas_data$p),"\n")
  if (max(ewas_data$p)<thrs.criteria){
    message("      **CAUTION** Your thresholding criteria may include CpGs not provided by the EWAS summary statistics")
    message("      Some CpGs may not be included in the episcore\n")
  }

  message("═= Step 2: CpG selection")
  selected_cpgs <- dplyr::filter(ewas_data, p <= thrs.criteria)
  if (nrow(selected_cpgs) == 0) {
    stop("No CpGs meet the threshold criteria.")
  }

  beta_filtered <- beta[beta$cpg %in% selected_cpgs$cpg, ]
  beta_filtered <- na.omit(beta_filtered)

  if (nrow(beta_filtered) < (1 - missingness) * nrow(selected_cpgs)) {
    message("      **CAUTION** Missingness in the beta matrix above ", missingness)
  }

  message("      CpG selection... Ready!\n")


  message("═= Step 3: epigenetic score calculation")
  message("      ", episcore_name, " episcore will be calculated\n")

  valid_indices <- which(selected_cpgs$cpg %in% beta_filtered$cpg)
  scores <- colSums(beta_filtered[-1] * scale(as.numeric(selected_cpgs$beta[valid_indices])))

  message("      Epigenetic score calculation complete\n")

  episcore <- setNames(data.frame(scores), episcore_name)

  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")


  message("═= Step 4: creating log file")

  log_message <- paste0(
    "Episcore Name: ", episcore_name, "\n",
    "\n",
    "Parameters: thrs.criteria = ", thrs.criteria, ", missingness = ", missingness, "\n",
    "Number of CpG in beta file: ", nrow(beta), "\n",
    "Number of CpG in EWAS summary statistics: ", nrow(ewas_data), "\n",
    "Number of CpGs after p-value thresholding: ", nrow(selected_cpgs), "\n",
    "CpG missingness ", round(((nrow(ewas_data))-nrow(beta_filtered)/nrow(ewas_data))/100,2),"%", "\n",
    "\n",
    "Time elapsed ", elapsed, "seconds"
  )
  writeLines(log_message, con = paste0(episcore_name, ".log"))

  message("      log file created\n")

  message("╔══════════════════════════════════════════════════════════════╗")
  message("║                    calculate.episcore END                    ║")
  message("╚══════════════════════════════════════════════════════════════╝")
  return(episcore)
}
