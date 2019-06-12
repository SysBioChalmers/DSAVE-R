#' averageIntoBins
#'
#' Help function for the DSAVE Score calculation. Places genes into
#' gene expression range bins and then calculates the mean for each
#' bin.
#'
#' averageIntoBins
#'
#' @param TPMData numeric matrix, tpmdata
#' @param logCV numeric vector, Gene Coefficient of Variation
#' @param templInfo list with 5 elements
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return a list with CV and meanGeneExpr

averageIntoBins <- function(TPMData, logcv, templInfo){
  avgRefExpr <- rowMeans(TPMData)
  numBins <- length(templInfo$binningInfo$binningInfo.means)
  cv <- rep(0,numBins)
  meanGeneExpr <- templInfo$binningInfo$binningInfo.means

  for(i in 1:numBins){
    #select the genes within the expression range
    sel <- avgRefExpr >= templInfo$binningInfo$binningInfo.lbs[i] & avgRefExpr <= templInfo$binningInfo$binningInfo.ubs[i]
    #y value in the graph
    if(sum(sel) == 0){
      cv[i] <- 0
      meanGeneExpr[i] <- 0
    }else{
      cv[i] <- mean(logcv[sel])
      meanGeneExpr[i] <- 2^mean(log2(avgRefExpr[sel]+0.05)) - 0.05
    }
  }
  return(list(cv = cv, meanGeneExpr = meanGeneExpr))
}





