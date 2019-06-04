#' BINNING INFORMATION
#'
#' Samples object elements and methods
#'
#' Sets the binning information
#'
#' @param data Numeric matrix
#' @param lbnonlog lower bound in original scale
#' @param ubnonlog upper bound in original scale
#' @param numPoints number of points for the distribution
#' @param poolSize number of genes to pool at each point in the distribution
#' @param step increas in lower and upper bound to get enough genes at each point in the distribution
#' @param toplim maximum number of iterations to obtain the frequency at each point in the distribution
#' @param reScale maximum number of iterations to obtain the frequency at each point in the distribution
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return a template object
#' @examples
#' \dontrun{
#' }


setBinningInfo = function(data = NULL, lbnonlog = NULL, ubnonlog = NULL,
                          reScale = FALSE, numPoints = 1000L, poolSize = 500L,
                          step = 5e-4, toplim = 1000L){
  stopifnot(is.numeric(data), is.matrix(data))
  stopifnot(is.integer(lbnonlog), length(lbnonlog) == 1)
  stopifnot(is.integer(ubnonlog), length(ubnonlog) == 1)
  stopifnot(lbnonlog < ubnonlog)
  stopifnot(is.integer(numPoints), length(numPoints) == 1)
  stopifnot(is.integer(poolSize), length(poolSize) == 1)
  stopifnot(is.integer(toplim), length(toplim) == 1)
  stopifnot(is.numeric(step), length(step) == 1)

  meansLog <- seq(log10(lbnonlog), log10(ubnonlog), length.out =numPoints)
  templInfo.binningInfo.means <- 10^meansLog

  data0 <- data
  if(reScale){
    cat("Re-scaling gene means \n")
    data <- tpmDSAVE(data)
    }

  gm <- as.matrix(rowMeans(data), ncol = 1)
  #gm <- tpmDSAVE(as.matrix(rowMeans(data0), ncol = 1))
  gmLog <- log10(gm)

  lbs <- rep(0, numPoints)
  ubs <- rep(0, numPoints)

  cat("Forming the distributions \n")

  for(i in 1:numPoints){
    if(i %in% seq(floor(numPoints / 10), numPoints, by = floor(numPoints / 10))) {
      cat(paste0("progress: ", round(i / numPoints, 2) * 100, "% \n"))
    }
    lb <- meansLog[i] - step
    ub <- meansLog[i] + step

    for(s in 1:toplim){
      sel <- gmLog >= lb & gmLog <= ub
      numGenes <- sum(sel)
      if(numGenes >= poolSize) {break}

      # increase bounds

      if(numGenes == 0) { meanOfCaughtGenes <- 0} else { meanOfCaughtGenes <- mean(gmLog[sel])}
      if(meanOfCaughtGenes > meansLog[i]){
        lb <- lb - step
      } else {
        ub <- ub + step
      }
    }
    lbs[i] <- 10^lb
    ubs[i] <- 10^ub
  }
  output <- list()
  output[["binningInfo.means"]] <- templInfo.binningInfo.means
  output[["binningInfo.lbs"]] <- lbs
  output[["binningInfo.ubs"]] <- ubs
  return(output)
}
