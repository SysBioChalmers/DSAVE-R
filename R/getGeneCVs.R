#' getGeneCVs
#'
#' Help function for DSAVECalcBTMScore. Calculates coefficient of variation (CV).
#'
#' @param ds numeric matrix
#' @param toLog default FALSE, if TRUE calculation will be done on log transformed data
#' @param logTPMAddon default 1, when using log transformed data, this number is added to the data before transforming to avoid log(0)
#' @author Juan Inda, <inda@@chalmers.se>
#' @return a vector
getGeneCVs <- function(ds, logTPMAddon, toLog=FALSE){
  ds <- tpmDSAVE(ds)
  if(toLog){
    totset <- log(ds + logTPMAddon)
    avgRefExpr <- rowMeans(totset)
    sd <- apply(totset, 1, sd)
    logcv <- sd / (avgRefExpr + 0.05)
  } else{
    totset <- ds
    avgRefExpr <- rowMeans(totset)
    sd <- apply(totset, 1, sd)
    logcv <- log( sd / (avgRefExpr + 0.05) +1)
  }

  return(logcv)
}



