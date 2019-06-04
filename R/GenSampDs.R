#' GenSampDs
#'
#' GenSampDs
#'
#' GenSampDs
#'
#' @param data numeric matrix, the input dataset (cell population)
#' @param countsPerGene counts per gene
#' @param edges edges
#' @importFrom graphics hist
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return list(logCV, variances)
#' @examples
#' \dontrun{
#' }
#'
GenSampDs <- function(data, countsPerGene, edges){
  new_data <- data
  new_data[,] <- 0
  numGenes <- dim(data)[1]
  for(i in 1:numGenes){
    if(countsPerGene[i] > 0){
      r <- runif(countsPerGene[i])
      a <- hist(r,edges, plot = F)$counts
      new_data[i,] <- a
    }
  }
  return(new_data)
}

