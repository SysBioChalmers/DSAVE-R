#' DSAVEPlotTotalVariation
#'
#' Plots the total variation as a function of pool size, with bulk variation as reference.
#'
#'
#' @param dsData Output from DSAVEGetTotalVariation or a list of such outputs
#' @param dsnames The names of the dataset(s), either a string or a list of strings
#' @importFrom methods is as
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return nothing
#' @examples
#' \dontrun{DSAVEPlotTotalVariation(dsData, dsNames)}
#'

DSAVEPlotTotalVariation <- function(dsData, dsNames){
  #detect if in param is a list - since the dsData is a list in itself, check for the Rs column...
  if (!(is.list(dsData) & is.null(dsData$Rs))) {
    dsData = list(dsData)
  }
  if (!is.list(dsNames)) {
    dsNames = list(dsNames)
  }
  if (length(dsData) != length(dsNames)) {
    stop("Mismatch in length between in params");
  }

  bulkMean1Vs1 <- rep(mean(bulkTotalVar1vs1[[1]][[2]]), 2)
  bulkMean4Vs4 <- rep(mean(bulkTotalVar4vs4[[1]][[2]]), 2)

  minVals = rep(0,length(dsData));
  maxVals = rep(0,length(dsData));

  sc = list();
  for(i in 1:length(dsData)) {
    sc = rbind(sc, cbind(dsData[[i]]$poolSizes, dsData[[i]]$Rs, dsNames[[i]], "sc"))
    minVals[i] = min(dsData[[i]]$poolSizes)
    maxVals[i] = max(dsData[[i]]$poolSizes)
  }
  minps = min(minVals)
  maxps = max(maxVals)

  #add the two bulk lines:
  sc = rbind(sc, cbind(c(minps,maxps), bulkMean1Vs1, "Single bulk sample", "bulk"))
  sc = rbind(sc, cbind(c(minps,maxps), bulkMean4Vs4, "Mean of 4 bulk samples", "bulk"))

  sc <- as.data.frame(sc)
  colnames(sc) <- c("PoolSize", "TotalVariation", "Dataset", "DataType")
  sc$PoolSize <- as.numeric(as.character(sc$PoolSize))
  sc$TotalVariation <- as.numeric(as.character(sc$TotalVariation))
  sc$Dataset = unlist(sc$Dataset)
  sc$DataType = unlist(sc$DataType)

  ggplot(sc, aes(x=PoolSize, y=TotalVariation, group = Dataset)) +
    ggtitle("Variation per Cell Pool Size, CPM > 0.5") +
    geom_line(aes(color = Dataset, linetype = DataType)) +
    ylim(0, max(sc$TotalVariation)) +
    theme(legend.title = element_text(size = 10),
          legend.justification = c(1, 1),
          legend.position = c(1, 1),
          legend.text = element_text(size = 8),
          axis.text.x = element_text(size = 10, angle = 90),
          axis.text.y = element_text(size = 10)) +
    #ylab(bquote('Variation ( ~ R[mean])')) +
    labs(y=expression(paste(Variation (R[mean] ))),
         x="Pool size (number of cells)") +
    guides(linetype = FALSE)

}
