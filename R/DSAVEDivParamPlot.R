#' DSAVEDivParamPlot
#'
#' Shows number of UMIs as function of divergence. Useful for determining thresholds for discarding bad cells.
#'
#' @param divData The output from the function DSAVEGetSingleCellDivergence
#' @param paramData Vector containing the param value for each cell.
#' @importFrom methods is as
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return -
#' @examples
#' \dontrun{DSAVEDivParamPlot(data, divData)}
#'

DSAVEDivParamPlot <- function(divData, paramData, paramName = "Param"){
  df = data.frame(div=divData$divs, param = paramData)
  ind = sort(paramData, index.return=T)$ix
  df = df[ind,]

  lf = loess(div~param, df)
  pred = predict(lf)

  ggplot(df, aes(x=div, y=param, colour = "Individual cell")) +
    geom_point() +
    ggtitle(paste0(paramName, " vs Cell Divergence")) +
    xlab("Divergence") + ylab(paramName) +
    theme(legend.justification = c(1, 1),
          legend.title=element_blank(),
          legend.position = c(1, 1),
          legend.text = element_text(size = 8),
          axis.text.x = element_text(size = 10, angle = 90),
          axis.text.y = element_text(size = 10)) +
    geom_path(aes(pred, param, colour="Mean divergence"))

}
