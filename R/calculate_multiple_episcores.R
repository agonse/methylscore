#' Calculate multiple episcores
#'
#' This function calculates multiple epigenetic scores using methylation information.
#' @description
#' This function calculates multiple epigenetic score based on methylation data. Two inputs are required: a beta matrix with CpG methylation and a list of external references with the association and effect values for each CpG
#'
#' @import readxl
#' @import dplyr
#' @param thrs.criteria Thresholding criteria for CpG selection.
#' @param beta Beta from sample CpG.
#' @param ewas.file Path and name to the list of EWAS files.
#' @param missingness Maximum tolerated missing data for CpG sites.
#' @return The calculated episcore.
#' @export
#'
calculate_multiple_episcores <- function(thrs.criteria = 0.05,
                                            beta.file = "test/betamatrix_test.txt",
                                            ewas.list.path = "test/sumstats/",
                                            missingness = 0.2,
                                            path = getwd()){

  sumstatslist <- list.files(ewas.list.path, full.names = T, pattern = "\\.txt$|\\.xlsx$")

  file_ext <- tools::file_ext(beta.file)
  beta <- if (file_ext == "xlsx") {
    read_xlsx(beta.file, col_names = TRUE)
  } else {
    read.table(beta.file, header = TRUE)
  }

  epi <- data.frame(id = colnames(beta)[-1])  # Assuming 'beta' is already a data.frame

  for (file in sumstatslist) {
    tryCatch({
      epi0 <- calculate_episcore(thrs.criteria, beta.file = beta.file, ewas.file = file, missingness = missingness, path = path)
      epi <- cbind(epi, setNames(epi0, tools::file_path_sans_ext(basename(file))))
    }, error = function(e) {
      error_log <- paste("Failed to process", basename(file), ":", e$message, "\n")
      all_logs <- c(all_logs, error_log)
    })
  }

  # Save all error logs to a single file
  writeLines(all_logs, file.path(path, "error.log"))
  return(epi)
}
