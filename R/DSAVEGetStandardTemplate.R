#' DSAVEGetStandardTemplate
#'
#' Gets the standard template that is redistributed with the package
#'
#' @importFrom graphics hist
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return the template
#' @examples
#' \dontrun{ DSAVEGetStandardTemplate()
#' }

DSAVEGetStandardTemplate <- function() {
  result = list()
  result$binningInfo <- DSAVE::bc2t_binningInfo
  result$UMIDistr <- as.numeric(DSAVE::bc2t_UMIdistribution)
  result$geneSet <- DSAVE::genesForTemplate
  result$fractionUpperOutliers <- 0.025
  result$fractionLowerOutliers <- 0.025

  return(result)
}
