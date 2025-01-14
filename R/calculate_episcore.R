#' Calculate a single episcore
#'
#' This function calculates one epigenetic score using methylation information.
#'
#' @description
#' This function calculates an epigenetic score based on methylation data. Two inputs are required: a beta matrix with CpG methylation and an external reference with the association and effect values for each CpG
#'
#' @import readxl
#' @import dplyr
#' @param thrs.criteria Thresholding criteria for CpG selection.
#' @param beta Beta from sample CpG.
#' @param ewas.path.file Path and name to the EWAS file.
#' @param missingness Maximum tolerated missing data for CpG sites.
#' @return The calculated episcore.
#' @export
#'
calculate_episcore <- function(thrs.criteria = 0.05,
                               beta.file = "",
                               ewas.path.file = "",
                               missingness = 0.2) {  
  
  start_time <- Sys.time()
  episcore_name <- tools::file_path_sans_ext(basename(ewas.path.file))     
  
  message("╔══════════════════════════════════════════════════════════════╗")
  message("║                  calculate.episcore START                    ║")
  message("╚══════════════════════════════════════════════════════════════╝")
  message("═= Step 1: beta matrix and load summary statistic file")                 
  
  if (!exists("beta.file")) stop("Beta file does not exist.")
  if (!dir.exists(ewas.list.path)) stop("EWAS directory does not exist.")
  
  file_ext <- tools::file_ext(ewas.path.file)
  ewas_data <- if (!is.character(ewas.path.file)) {
    stop(paste(
      "Error: You must specify the full path and file name with extension.",
      "For example: your/ewas/sumstats/path/sumstats.extension",
      sep = "\n"
    ))
  } else if (file_ext == "xlsx") {
    read_xlsx(ewas.path.file, col_names = TRUE)
  } else {
    read.table(ewas.path.file, header = TRUE)
  }
  message("      Reading summary statistics .", file_ext, " file... Done!")
  message("      p value thresholding is set at ", thrs.criteria)
  
  message("      Maximum p value provided by EWAS summary statistics is ", max(ewas_data$p),"\n")
  if (max(ewas_data$p)<thrs.criteria){
    message("      **CAUTION** Your thresholding criteria may include CpGs not provided by the EWAS summary statistics")
    message("      Some CpG meeting the p-value thresholding criteria may not be included in the episcore\n")
    log_message_thrs="**CAUTION** Some CpG meeting the p-value thresholding criteria may not be included in the episcore."
  } else{
    log_message_thrs=""
  }
  
  message("═= Step 2: CpG selection")
  selected_cpgs <- dplyr::filter(ewas_data, p <= thrs.criteria)
  if (nrow(selected_cpgs) == 0) {
    stop("No CpGs meet the threshold criteria.")
  }
  
  if ("data.table" %in% class(beta)) {
    setDF(beta)
  }
  
  beta_filtered <- beta.file[beta.file$cpg %in% selected_cpgs$cpg, ]
  beta_filtered <- na.omit(beta_filtered)
  
  if ((nrow(beta_filtered)/nrow(selected_cpgs)) > missingness) {
    message("      **CAUTION** Missingness in the beta matrix above ", missingness)
    log_message_miss="**CAUTION** Missingness exceeds defined parameter."
  } else {
    log_message_miss=""
  }
  message("      CpG selection... Ready!\n")
  
  message("═= Step 3: epigenetic score calculation")
  
  valid_indices <- which(selected_cpgs$cpg %in% beta_filtered$cpg)
  scores <- colSums(beta_filtered[,-1] * scale(as.numeric(selected_cpgs$beta[valid_indices])))
  message("      Epigenetic score calculation complete\n")
  
  episcore <- setNames(data.frame(scores), episcore_name)
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")
  
  message("═= Step 4: creating log file")
  
  logs_dir <- file.path(getwd(), "logs")
  if (!dir.exists(logs_dir)) {
    dir.create(logs_dir)
  }
  log_path <- file.path(logs_dir, paste0(episcore_name, ".log"))

  log_message <- paste0(
    episcore_name, " calculaton Log\n",
    "=========================\n",
    "\n",
    "Episcore Name: ", episcore_name, "\n",
    "Timestamp: ", Sys.time(), "\n",
    "\n",
    "Parameters:\n",
    "  - Thresholding criteria (p-value <=): ", thrs.criteria, "\n",
    "  - CpG missingness: ", round(nrow(beta_filtered) / nrow(selected_cpgs), 2), "\n",
    log_message_thrs, "\n",
    log_message_miss, "\n",
    "\n",
    "CpG Count Summary:", "\n",
    "Number of CpG in beta file: ", nrow(beta.file), "\n",
    "Number of CpG in EWAS summary statistics: ", nrow(ewas_data), "\n",
    "Number of CpG in EWAS summary after p-value thresholding: ", nrow(selected_cpgs), "\n",
    "Number of CpG in episcore: ", length(valid_indices), "\n",
    "\n",
    "Time elapsed ", round(elapsed,2), " seconds\n",
    "\n",
    "Log File Path: ", log_path, "\n"
  )
  writeLines(log_message, con = log_path)
  
  message("      log file created\n")
  
  message("╔══════════════════════════════════════════════════════════════╗")
  message("║                    calculate.episcore END                    ║")
  message("╚══════════════════════════════════════════════════════════════╝")
  return(episcore)
}

  
