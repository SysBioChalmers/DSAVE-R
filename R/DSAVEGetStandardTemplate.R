#' DSAVEGetStandardTemplate
#'
#' Gets the standard template that is distributed with the package, requiring > 2000 cells.
#' There are variants for 2000 (standard), 1000 and 500 cells.
#'
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return the template

DSAVEGetStandardTemplate <- function() {
  result = list()
  result$binningInfo <- DSAVE::binInf2000
  result$UMIDistr <- DSAVE::UMIDistr2000
  result$geneSet <- DSAVE::genesForTemplate
  result$fractionUpperOutliers <- 0.025
  result$fractionLowerOutliers <- 0.025

  return(result)
}

#' DSAVEGetStandardTemplate1000
#'
#' Gets the standard template that is distributed with the package, requiring > 1000 cells.
#'
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return the template

DSAVEGetStandardTemplate1000 <- function() {
  result = list()
  result$binningInfo <- DSAVE::binInf1000
  result$UMIDistr <- DSAVE::UMIDistr1000
  result$geneSet <- DSAVE::genesForTemplate
  result$fractionUpperOutliers <- 0.025
  result$fractionLowerOutliers <- 0.025

  return(result)
}

#' DSAVEGetStandardTemplate500
#'
#' Gets the standard template that is distributed with the package, requiring > 500 cells.
#'
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return the template
DSAVEGetStandardTemplate500 <- function() {
  result = list()
  result$binningInfo <- DSAVE::binInf500
  result$UMIDistr <- DSAVE::UMIDistr500
  result$geneSet <- DSAVE::genesForTemplate
  result$fractionUpperOutliers <- 0.025
  result$fractionLowerOutliers <- 0.025

  return(result)
}
