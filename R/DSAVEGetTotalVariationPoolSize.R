#' DSAVEGetTotalVariationPoolSize
#'
#' Calculates the total variation with a specific pool size.
#'
#' Calculates the average pairwise total variation between two pools of cells with specific size,
#' with a TPM filtration on the genes.
#'
#' @param data a numeric matrix
#' @param upperBound filters out the genes with mean expression higher than this value.
#' Defaults to 100,000, which in practice includes almost all genes.
#' @param lowerBound filters out the genes with mean expression lower than this value. Defaults to 0.5.
#' @param poolSizes list of numeric pool size values, all should be maximum half of the total
#' size of the data in the data object
#' @param repetitionPerSize Positive integer describing how many iterations that are used.
#' Defaults to 30.
#' @param rescale logical, determine if the data should be rescaled. Should normally be TRUE.
#' pool size
#' @param allowSubSampling [optional] If true, the function will reduce the dataset used to
#' 2 times the max pool size used. This is to save memory usage for large datasets. Defaults to FALSE.
#' @param silent (optional) If true, no progress bar is shown. Defaults to FALSE
#' @importFrom textTinyR sparse_Means
#' @export
#' @author Juan Inda, <inda@@chalmers.se>, Johan Gustafsson, <gustajo@@chalmers.se>
#' @return a list(poolSizes = numeric vector, Rs = numeric vector (R metric for each pool size))
#' @examples
#' \dontrun{DSAVEGetTotalVariationPoolSize(data, poolSize = 50, upperBound = 1000,
#' lowerBound = 0.01, na.rm = TRUE, seed = 1, repetitionPerSize = 30, 0.5)}
#'

DSAVEGetTotalVariationPoolSize <- function(data, upperBound = 1e5, lowerBound = 5e-1,
                                           poolSizes = NA,
                                           repetitionPerSize = 30L, rescale = TRUE, allowSubSampling = FALSE,
                                           silent=FALSE){
  #print("Control parameters")
  stopifnot((is.matrix(data) |  is(data, 'sparseMatrix')),
            is.numeric(upperBound), length(upperBound) == 1,
            is.numeric(lowerBound), length(lowerBound) == 1,
            is.logical(rescale), length(rescale) == 1,
            is.integer(repetitionPerSize), repetitionPerSize > 0)

  if (is.na(poolSizes)) {
    poolSizes = c(100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 10000)
    poolSizes = poolSizes[poolSizes <= ((dim(data)[2])/2)] #prefilter here to hide warning
  } else {
    stopifnot(is.numeric(poolSizes))
  }
  #remove too large pool sizes
  sel = poolSizes <= ((dim(data)[2])/2)
  poolSizes = poolSizes[sel]
  if ((!all(sel)) & !silent) {
    print("Warning: Removed pool sizes - they were larger than half the dataset size")
  }
  numPS = length(poolSizes)
  if (length(poolSizes) == 0) {
    stop("Poolsizes has length zero.")
  }


  if (!silent) {
    pb <- progress_bar$new(format = "Calculating total variation [:bar] :percent eta: :eta",
                           total = numPS + 1, clear = FALSE)
    pb$tick();
  }


  if (!is.matrix(data)) {
    #Sparse matrix, special code for this

    #If the dataset has more cells than 2 times the largest pool size, reduce the dataset to save memory when converting to a full matrix
    if (allowSubSampling & (ncol(data) > 2*max(poolSizes))) {
      data = data[, sample(ncol(data), 2*max(poolSizes)), drop=FALSE]
    }
    rm = sparse_Means(data, rowMeans = T)
    sel = rm != 0
    data = data[sel, , drop=FALSE]
    #Normalize to TPM/CPM
    if(rescale) {
      colsms = sparse_Sums(data, rowSums = FALSE)
      data = Matrix::t(Matrix::t(data)*10^6/colsms)
    }
    #Test:
    #sparse_Sums(data, rowSums = FALSE) - ok

    results = rep(0, length(poolSizes))

    #Filter genes
    #row_means <- rowMeans(data)
    row_means <- sparse_Means(data, rowMeans = TRUE)
    sel = (row_means >= lowerBound) & (row_means <= upperBound)
    data <- data[sel, , drop=FALSE]
    numCells = ncol(data)
    ind = 1:ncol(data)


    #loop through all sizes
    for (i in 1:length(poolSizes)) {
      meanSum = 0
      for (j in 1:repetitionPerSize) {
        #Don't use sample with a vector, it behaves differently when the vector is of size 1!
        a <- sample(numCells, poolSizes[i])
        indNoA <- ind[-a, drop=FALSE]
        b = indNoA[sample(numCells-poolSizes[i], poolSizes[i]),drop=FALSE]
        a <- sparse_Means(data[,a,drop=FALSE], rowMeans = TRUE)
        b <- sparse_Means(data[,b,drop=FALSE], rowMeans = TRUE)
        meanSum = meanSum + mean(abs(log((a + 0.05) / (b + 0.05))))
      }
      results[i] = meanSum/repetitionPerSize;

      if (!silent) {
        pb$tick()
      }
    }
  } else {
    #Normalize to TPM/CPM
    if(rescale) {
      data <- tpmDSAVE(data)
    }

    #Test:
    #colSums(data) - ok

    results = rep(0, length(poolSizes))

    #Filter genes
    row_means <- rowMeans(data)
    sel = (row_means >= lowerBound) & (row_means <= upperBound)
    data <- data[sel, , drop=FALSE]
    numCells = ncol(data)
    ind = 1:ncol(data)


    #loop through all sizes
    for (i in 1:length(poolSizes)) {
      meanSum = 0
      for (j in 1:repetitionPerSize) {
        #Don't use sample with a vector, it behaves differently when the vector is of size 1!
        a <- sample(numCells, poolSizes[i])
        indNoA <- ind[-a, drop=FALSE]
        b = indNoA[sample(numCells-poolSizes[i], poolSizes[i]),drop=FALSE]
        a <- rowMeans(data[,a,drop=FALSE])
        b <- rowMeans(data[,b,drop=FALSE])
        meanSum = meanSum + mean(abs(log((a + 0.05) / (b + 0.05))))
      }
      results[i] = meanSum/repetitionPerSize;

      if (!silent) {
        pb$tick()
      }
    }
  }


  if (!silent) {
    pb$terminate()
  }

  return(list(poolSizes = poolSizes, Rs = results))
}
