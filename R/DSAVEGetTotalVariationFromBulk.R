#' DSAVEGetTotalVariationFromBulk
#'
#' Calculate total variation from bulk sample
#'
#' Calculates the average pairwise total variation between a list of bulk samples, with a TPM
#' filtration on the genes. If pool4samples is true, it will compare the mean
#' of 4 randomly selected samples with the mean of 4 others. This expects a
#' list of 8 samples.
#'
#' @param data a numericl matrix
#' @param pool4samples logic value, if TRUE, the function will compare the mean
#' of 4 randomly selected samples with the mean of 4 others
#' @param upperBoundTPM filters out the genes with mean expression higher than this value
#' @param lowerBoundTPM filters out the genes with mean expression lower than this value
#' @param na.rm logical, sends na.rm value to the rowMeans function
#' @param rescale logical, determine if the data should be rescaled
#' @param seed positive integer for to set the seed of the samples
#' @param nComb maximum number of combinations of sub samples of 4
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return numeric vector
#' @examples
#' \dontrun{
#' DSAVEGetTotalVariationFromBulk(sample, pool4samples = FALSE, upperBoundTPM = 100000,
#' lowerBoundTPM = 0.5)
#' }
#'

DSAVEGetTotalVariationFromBulk <- function(data, pool4samples, upperBoundTPM = 1e5,
                                           lowerBoundTPM = 5e-1, na.rm = TRUE,
                                           nComb = 1000L, seed = NULL, rescale = TRUE){
  stopifnot(is.matrix(data),
            is.logical(pool4samples), length(pool4samples) == 1,
            is.numeric(upperBoundTPM), length(upperBoundTPM) == 1,
            is.numeric(lowerBoundTPM), length(lowerBoundTPM) == 1,
            is.logical(na.rm), length(na.rm) == 1,
            is.logical(rescale), length(rescale) == 1,
            (is.null(seed) | is.integer(seed)), is.integer(nComb))

  if(rescale) {data <- tpmDSAVE(data)}
  row_means <- rowMeans(data, na.rm = na.rm)
  genesToKeep <- which(row_means >= lowerBoundTPM & row_means <= upperBoundTPM)
  sample <- data[genesToKeep,]
  numSamp <- dim(data)[2]
  numGenes <- dim(data)[1]
  diffs <- matrix(0.0, nrow = numGenes, ncol = numSamp*(numSamp-1)/2)


  if(!pool4samples){
    combs.list <- list()
    ix <- 1
    for(i in 1:(numSamp-1)){
      for(j in (i+1):numSamp){
        combs.list[[ix]] <- list(a = i, b = j)
        ix = ix + 1
      }
    }

    mean_vector <- sapply(combs.list, function(v){
      a <- data[,v$a]
      b <- data[,v$b]
      mean(abs(log((a + 0.05) / (b + 0.05))))
    })


  } else{
    if(numSamp < 8){
      stop("Variation 4 samples on 4 more only works with min samples in object = 8")
    } else {
      #here the average of 4 samples is compared to the average of another 4
      #get all permutations of four

      ind <- as.integer(1:numSamp)
      k <- 4
      maxComb <- choose(numSamp, k)
      ix <- 1
      combs.list <- list()
      ind.tmp <- ind
      if(!is.null(seed)) set.seed(seed)
      while(ix <= nComb & ix <= maxComb){
        if(length(ind.tmp) >= 2*k){
          a <- sample(ind.tmp, k)
          ind.tmp <- ind.tmp[-a]
          b <- sample(ind.tmp, k)
          ind.tmp <- ind.tmp[-b]
          combs.list[[ix]] <- list(a = a, b = b)
          ix <- ix + 1
        } else {
          ind.tmp <- ind
        }
      }

      mean_vector <- sapply(combs.list, function(v){
        a <- rowMeans(data[,v$a], na.rm = na.rm)
        b <- rowMeans(data[,v$b], na.rm = na.rm)
        mean(abs(log((a + 0.05) / (b + 0.05))))
      })
    }
  }
  return(mean_vector)
}
