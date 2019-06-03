#' SET UMI DISTRIBUTION
#'
#' Set UMI distribution
#'
#' Creates a UMI distribution
#'
#' @param data Numeric matrix
#' @param numUMI average number of molecules per cell
#' @importFrom graphics hist
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return a template object
#' @examples
#' \dontrun{
#' }
setUMIDistribution = function(data = NULL, numUMI = 750){
  stopifnot(is.numeric(data), is.matrix(data))
  stopifnot(is.integer(numUMI), length(numUMI) == 1)
  totalUMI <- sum(data)
  sumUMI <- colSums(data)
  totalTargetUMI <- dim(data)[2] * numUMI
  toRemove <- totalUMI - totalTargetUMI
  UMIDistr <- sumUMI
  if(toRemove < 0){
    stop("Could not reach target UMIs per cell due to that the UMIs in the template dataset were to few")
  } else {
      #cat(paste("removing", toRemove, "\n"))
      #cat(paste("sampling \n"))
      indToRem <- sample(totalUMI, toRemove)
      #cat(paste("forming edges \n"))
      edges <- c(0,cumsum(UMIDistr) + 0.01)
      #cat(paste("substracting \n"))
      subtr <- hist(indToRem, breaks = edges, plot = F)$counts
      UMIDistr <- UMIDistr - subtr
      UMIDistr <- sort(UMIDistr)
  }
  return(UMIDistr)
}
