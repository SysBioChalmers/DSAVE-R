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
#' @param upperBound filters out the genes with mean expression higher than this value
#' @param lowerBound filters out the genes with mean expression lower than this value
#' @param na.rm logical, sends na.rm value to the rowMeans function
#' @param nComb ?
#' @param rescale logical, determine if the data should be rescaled to TPM
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return numeric vector
#' @examples
#' DSAVEGetTotalVariationFromBulk(sample, pool4samples = FALSE, upperBound = 100000,
#' lowerBound = 0.5)
#'

DSAVEGetTotalVariationFromBulk <- function(data, pool4samples, upperBound = 1e5,
                                           lowerBound = 5e-1, na.rm = TRUE,
                                           nComb = 1000L, rescale = TRUE){
  stopifnot(is.matrix(data),
            is.logical(pool4samples), length(pool4samples) == 1,
            is.numeric(upperBound), length(upperBound) == 1,
            is.numeric(lowerBound), length(lowerBound) == 1,
            is.logical(na.rm), length(na.rm) == 1,
            is.logical(rescale), length(rescale) == 1)

  if(rescale) {data <- tpmDSAVE(data)}
  row_means <- rowMeans(data, na.rm = na.rm)
  genesToKeep <- which(row_means >= lowerBound & row_means <= upperBound)
  data <- data[genesToKeep, , drop=FALSE]
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
    if(numSamp != 8){
      stop("Variation 4 samples on 4 more is only implemented for 8 samples in total")
    } else {
      #here the average of 4 samples is compared to the average of another 4
      #get all permutations of four

      comb = combn(8,4);
      numComb = dim(comb)[2];
      combs.list <- list();
      indices = 1:8;
      for (i in 1:numComb) {
        a = comb[,i];
        b = indices[-a];
        combs.list[[i]] <- list(a = a, b = b)
      }

      mean_vector <- sapply(combs.list, function(v){
        a <- rowMeans(data[,v$a, drop=FALSE], na.rm = na.rm)
        b <- rowMeans(data[,v$b, drop=FALSE], na.rm = na.rm)
        mean(abs(log((a + 0.05) / (b + 0.05))))
      })
    }
  }
  return(mean(mean_vector))
}
