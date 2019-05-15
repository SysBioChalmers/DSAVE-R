#' DSAVEGetTotalVariationPoolSize
#'
#' Calculate total variation with a specific pool size
#'
#' Calculates the average pairwise total variation between two pools of cells with specific size,
#' with a TPM filtration on the genes.
#'
#' @param sample object sample to get the data from
#' @param poolSize numeric value, it should be maximum half of the total size of the data
#' in the sample object
#' @param upperBoundTPM filters out the genes with mean expression higher than this value
#' @param lowerBoundTPM filters out the genes with mean expression lower than this value
#' @param na.rm logical, sends na.rm value to the rowMeans function
#' @param rescale logical, determine if the data should be rescaled
#' @param seed positive integer for to set the seed of the samples
#' @param repetitionPerSize positive integer how many times the variation is calculated at each
#' pool size
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return numeric vector
#' @examples
#' \dontrun{
#' DSAVEGetTotalVariationPoolSize(samples, poolSize = 50, upperBoundTPM = 1000,
#' lowerBoundTPM = 0.01, na.rm = TRUE, seed = 1, repetitionPerSize = 30, 0.5)
#' }
#'

DSAVEGetTotalVariationPoolSize <- function(sample, poolSize = 4, upperBoundTPM = 1e5,
                                           lowerBoundTPM = 5e-1, na.rm = TRUE, seed = NULL,
                                           repetitionPerSize = 30L, rescale = TRUE){
  #print("Control parameters")
  stopifnot(class(sample)[1] == "Samples",
            is.numeric(poolSize), length(poolSize) == 1,
            is.numeric(upperBoundTPM), length(upperBoundTPM) == 1,
            is.numeric(lowerBoundTPM), length(lowerBoundTPM) == 1,
            is.logical(na.rm), length(na.rm) == 1,
            is.logical(rescale), length(rescale) == 1,
            (is.null(seed) | is.integer(seed)),
            is.integer(repetitionPerSize), repetitionPerSize > 0)

  if(poolSize > floor(sample$numberSamples/2)){
    stop("Cannot make non-overlapping pool of cells. Choose a smaller poolSize.")
  } else {
    #print("Re-scaling")
    if(rescale) sample$reScale()
    #print("Filtering genes")
    row_means <- rowMeans(sample$data, na.rm = na.rm)
    genesToKeep <- names(which(row_means >= lowerBoundTPM & row_means <= upperBoundTPM))
    sample <- sample$geneSubset(genesToKeep = genesToKeep);
    numSamp <- sample$numberSamples
    numGenes <- sample$numberGenes

    #print("Creating combinations")
    ind <- as.integer(1:numSamp)
    ix <- 1
    combs.list <- list()
    ind.tmp <- ind
    if(!is.null(seed)) set.seed(seed)
    while(ix <= repetitionPerSize ){#  & ix <= maxComb){
      if(length(ind.tmp) >= 2*poolSize){
        a <- sample(ind.tmp, poolSize)
        ind.tmp <- ind.tmp[-a]
        b <- sample(ind.tmp, poolSize)
        ind.tmp <- ind.tmp[-b]
        combs.list[[ix]] <- list(a = a, b = b)
        ix <- ix + 1
      } else {
        ind.tmp <- ind
      }
    }
    #print("Calculating variation")
    mean_vector <- sapply(combs.list, function(v){
      a <- rowMeans(sample$data[,v$a], na.rm = na.rm)
      b <- rowMeans(sample$data[,v$b], na.rm = na.rm)
      mean(abs(log((a + 0.05) / (b + 0.05))))
    })

  }
  return(mean_vector)
}
