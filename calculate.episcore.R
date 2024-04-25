#' Calculate episcore
#'
#' This function calculates epigenetic scores using methylation information.
#'
#' @import readxl
#' @import dplyr
#' @param thrs.criteria Thresholding criteria for CpG selection.
#' @param beta Beta from sample CpG.
#' @param ewas.file Path and name to the EWAS file.
#' @param missingness Maximum tolerated missing data for CpG sites.
#' @return The calculated episcore.
#' @export
calculate.episcore=function(thrs.criteria=0.05,
                            beta="methylscore/betamatrix_test.txt",
                            ewas.file="methylscore/sumstats_test.txt",
                            missingess=0.20)
{
  message("╔══════════════════════════════════════════════════════════════╗")
  message("║                  calculate.episcore START                    ║")
  message("╚══════════════════════════════════════════════════════════════╝")
  message("═= Step 1: load summary statistic file")

  file_ext <- tools::file_ext(ewas.file)
  if (file_ext == "xlsx") {
    message("      Reading ", ewas.file)
    cpgms <- read_xlsx(ewas.file,col_names = TRUE)
  } else {
    message("      Reading ", ewas.file)
    cpgms <- read.table(ewas.file,header = TRUE)
  }
  message("      Reading .", file_ext, " file... Done!")
  message("      p value thresholding is set at ", thrs.criteria)
  message("      Maximum p value provided by EWAS summary statistics is ", max(cpgms$p),"\n")
  if (max(cpgms$p)<thrs.criteria){
    message("      **CAUTION** Your thresholding criteria can include CpGs not provided by the EWAS summary statistics")
    message("      Some CpGs may not be included in the episcore\n")
  }

  message("═= Step 2: CpG selection")

  cpgms=cpgms%>%subset(p<=thrs.criteria)

  if (!is.data.frame(beta)){
    stop("      Error: The beta matrix must be provided as a data frame.")
  }

  beta0=beta[beta$cpg %in% cpgms$cpg,]
  beta0=na.omit(beta0)

  episcore.name=tools::file_path_sans_ext(basename(ewas.file))

  if (nrow(beta0)<1 | nrow(cpgms)<1){
    scorestring=rep("NA",times=ncol(beta0)-1)
    episcdummy <- data.frame(setNames(list(scorestring), episcore.name))
    epi=cbind.data.frame(epi,episcdummy)
    return(epi)
    stop("      Error: No available CpGs to create episcore")
  } else if ((nrow(beta0)/nrow(cpgms)*100)<(1-missingess)){
    message("      **CAUTION** CpG missigness exceeds ", missingess)
  }

  cpgms0=cpgms[cpgms$cpg %in% beta0$cpg,]
  cpgms0=cpgms0[match(beta0$cpg,cpgms0$cpg),]
  cpgms0$beta0=scale(cpgms0$beta)

  message("      CpG selection... Ready!\n")

  message("═= Step 3: epigenetic score calculation")
  message("      ", episcore.name, " episcore will be calculated\n")
  scorestring = numeric()
  for (j in 2:(ncol(beta0))){
    singleid=0
    for (i in 1:nrow(cpgms0)){
      singlecpg=cpgms0$beta[i]*beta0[[i,j]]
      singleid=singleid+singlecpg
    }
    scorestring=c(scorestring,singleid)
  }

  episc <- data.frame(setNames(list(scorestring), episcore.name))
  episcore=cbind.data.frame(epi,episc)
  colnames(episcore)=c(colnames(episcore)[1],episcore.name)

  message("      Epigenetic score calculation complete\n")


  message("                 *****   Final Report      *****         ")
  message("  episcore ",episcore.name," calculation successful.\n")
  message(paste0("  There are ", nrow(cpgms)," CpGs with p value <= ", thrs.criteria, " in the reference.\n"),
                  "  ",nrow(beta0), "(", round(nrow(beta0)/nrow(cpgms)*100,2), "%) CpGs are present in the beta matrix.\n")

  message("╔══════════════════════════════════════════════════════════════╗")
  message("║                    calculate.episcore END                    ║")
  message("╚══════════════════════════════════════════════════════════════╝")

  return(episcore)
}
