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

  load("data/bc2t_binningInfo.rda")
  load("data/bc2t_UMIdistribution.rda")
  load("data/genesForTemplate.rda")

  result = list()
  result$binningInfo = bc2t_binningInfo
  result$UMIDistr = as.numeric(bc2t_UMIdistribution)
  result$geneSet = genesForTemplate
  result$fractionUpperOutliers = 0.025
  result$fractionLowerOutliers = 0.025

  return(result)
}
