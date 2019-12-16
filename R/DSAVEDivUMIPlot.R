#' DSAVEDivUMIPlot
#'
#' Shows number of UMIs as function of divergence. Useful for determining thresholds for discarding bad cells.
#'
#' @param data The UMI count data
#' @param divData The output from the function DSAVEGetSingleCellDivergence
#' @importFrom methods is as
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return nothing
#' @examples
#' \dontrun{DSAVEDivUMIPlot(data, divData)}
#'

DSAVEDivUMIPlot <- function(data, divData){
  numUMIs <- colSums(as.matrix(data))
  df = data.frame(div=divData$lls, numUMIs = numUMIs)
  ind = sort(numUMIs, index.return=T)$ix
  df = df[ind,]

  lf = loess(div~numUMIs, df)
  pred = predict(lf)

  ggplot(df, aes(x=div, y=numUMIs, colour = "Individual cell")) +
    geom_point() +
    ggtitle("UMI Counts vs Cell Divergence") +
    xlab("Log-likelihood") + ylab("UMI counts") +
    theme(legend.justification = c(1, 1),
          legend.title=element_blank(),
          legend.position = c(1, 1),
          legend.text = element_text(size = 8),
          axis.text.x = element_text(size = 10, angle = 90),
          axis.text.y = element_text(size = 10)) +
    geom_path(aes(pred, df$numUMIs, colour="Mean log-likelihood"))

}
