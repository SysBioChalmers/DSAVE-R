#' DSAVEPlotDivergence
#'
#' Creates an interactive plotly plot showing cell divergence and gene divergence per cell
#'
#' Calculates the average pairwise total variation between two pools of cells with specific size,
#' with a TPM filtration on the genes.
#'
#' @param data The UMI count data
#' @param divData The output from the function DSAVEGetSingleCellDivergence
#' @importFrom plotly plot_ly layout
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return nothing
#' @examples
#' \dontrun{DSAVEPlotDivergence(data, divData)}
#'

DSAVEPlotDivergence <- function(data, divData){

  if (dynutils::is_sparse(data)) {
    y = sparse_Sums(data, rowSums=F)
  } else {
    y = colSums(as.matrix(data))
  }

  x=divData$divs

  #Generate suitable text for hover
  numCells = length(divData$divs)
  ids = 1:numCells
  templGene <- rep("", numCells)
  templVal <- rep(0,numCells)

  t5N = data.frame(gene1=templGene,gene2=templGene,gene3=templGene, gene4=templGene, gene5=templGene, stringsAsFactors = F)
  t5V = data.frame(val1=templVal, val2=templVal, val3=templVal, val4=templVal, val5=templVal)

  #row.names(divData$geneDivs) = as.character(row.names(divData$geneDivs))
  #fill using a loop
  #for (i in ids) {
  #  a = sort(divData$geneDivs[,i], decreasing = T)
  #  t5N[i,] = names(a[1:5])
  #  t5V[i,] = a[1:5]
  #}
  t5N = t(divData$geneDivGenes)
  t5V = t(divData$geneDivVals)

  p = plot_ly(x=x, y=y, type="scatter", mode="markers",
              text = ~paste(ids, '<br>',
                            t5N[,1], " ", t5V[,1], "<br>",
                            t5N[,2], " ", t5V[,2], "<br>",
                            t5N[,3], " ", t5V[,3], "<br>",
                            t5N[,4], " ", t5V[,4], "<br>",
                            t5N[,5], " ", t5V[,5], "<br>")) %>%
    layout(xaxis = list(title="Cell divergence"), yaxis = list(title="UMI counts"), title="Interactive Divergence Plot")
  p

  # p = plot_ly(data = df, x=~V1, y=~V2, type="scatter", mode="markers",
  #             text = ~paste(ids, '<br>',
  #                           t5N[,1], " ", t5V[,1], "<br>",
  #                           t5N[,2], " ", t5V[,2], "<br>",
  #                           t5N[,3], " ", t5V[,3], "<br>",
  #                           t5N[,4], " ", t5V[,4], "<br>",
  #                           t5N[,5], " ", t5V[,5], "<br>")) %>%
  #   layout(xaxis = list(title="Cell divergence"), yaxis = list(title="UMI counts"), title="Interactive Divergence Plot")
  # p

}
