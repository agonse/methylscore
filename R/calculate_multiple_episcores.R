#' Calculate multiple episcores
#'
#' This function calculates multiple epigenetic scores using methylation information.
#' @description
#' This function calculates multiple epigenetic score based on methylation data. It requires a beta matrix and a directory containing EWAS files.
#'
#' @import readxl
#' @import dplyr
#' @param thrs.criteria Thresholding criteria for CpG selection.
#' @param beta.file A preloaded beta matrix (data.frame) with CpG methylation values.
#' @param ewas.list.path Path to the directory of EWAS files.
#' @param missingness Maximum tolerated missing data for CpG sites.
#' @return A data.frame with calculated episcores for all EWAS files.
#' @exportewas.list.pathewas.list.path
#'
calculate_multiple_episcores <- function(thrs.criteria = 0.05,
                                         beta.file = "",
                                         ewas.list.path = "",
                                         missingness = 0.2){
  
  if (!dir.exists(ewas.list.path)) stop("Error: EWAS directory does not exist.")

  sumstatslist <- list.files(ewas.list.path, full.names = T, pattern = "\\.txt$|\\.xlsx$")
  
  epi <- data.frame(id = colnames(beta.file)[-1])
  all_logs <- character()
  
  for (file in sumstatslist) {
    tryCatch({
      epi0 <- calculate_episcore(thrs.criteria = thrs.criteria,
                                 beta.file = beta.file,
                                 ewas.path = ewas.list.path,
                                 ewas.file = basename(file),
                                 missingness = missingness)
      episcore_name=tools::file_path_sans_ext(basename(file))
      epi <- cbind(epi, epi0[[1]])
      colnames(epi)[ncol(epi)] <- episcore_name
    }, error = function(e) {
      error_log <- paste("Failed to process", basename(file), ":", e$message, "\n")
      all_logs <- c(all_logs, error_log)
    })
  }
  
  if (length(all_logs) > 0) {
    writeLines(all_logs, file.path(path, "error.log"))
  }
  return(epi)
}
