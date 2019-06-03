#' GetGeneCVs
#'
#' GetGeneCVs
#'
#' Calculate coefficient of variation
#'
#' @param data numeric matrix
#' @param toLog default FALSE, if TRUE calculation will be done on log transformed data
#' @param logTPMAddon default 1, when using log transformed data, this number is added to the data before transforming to avoid log(0)
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return a vector
#' @examples
#' \dontrun{
#' }
GetGeneCVs <- function(ds, logTPMAddon, toLog=FALSE){
  ds_red <- tpmDSAVE(ds)
  if(toLog){
    totset <- log(ds + logTPMAddon)
  } else{
    totset <- ds
  }
  avgRefExpr <- rowMeans(totset)
  sd <- apply(totset, 1, sd)
  logcv <- sd / (avgRefExpr + 0.01)
  return(logcv)
}



