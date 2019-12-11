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
#' @param na.rm logical, sends na.rm value to the rowMeans function, should normally be true.
#' @param seed positive integer for to set the seed of the datas
#' @param repetitionPerSize Positive integer describing how many iterations that are used.
#' Defaults to 30.
#' @param rescale logical, determine if the data should be rescaled. Should normally be TRUE.
#' pool size
#' @param silent (optional) If true, no progress bar is shown. Defaults to FALSE
#' @importFrom methods is as
#' @export
#' @author Juan Inda, <inda@@chalmers.se>, Johan Gustafsson, <gustajo@@chalmers.se>
#' @return numeric vector
#' @examples
#' \dontrun{DSAVEGetTotalVariationPoolSize(data, poolSize = 50, upperBound = 1000,
#' lowerBound = 0.01, na.rm = TRUE, seed = 1, repetitionPerSize = 30, 0.5)}
#'

DSAVEGetTotalVariationPoolSize <- function(data, upperBound = 1e5, lowerBound = 5e-1,
                                           poolSizes = NA, na.rm = TRUE, seed = NULL,
                                           repetitionPerSize = 30L, rescale = TRUE,
                                           silent=FALSE){
  #print("Control parameters")
  stopifnot((is.matrix(data) |  is(data, 'sparseMatrix')),
            is.numeric(upperBound), length(upperBound) == 1,
            is.numeric(lowerBound), length(lowerBound) == 1,
            is.logical(na.rm), length(na.rm) == 1,
            is.logical(rescale), length(rescale) == 1,
            (is.null(seed) | is.integer(seed)),
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


  #if sparse, prefilter genes to save some memory and then convert to non-sparse
  if (!is.matrix(data)) {
    rm = sparse_Means(data, rowMeans = T)
    sel = rm != 0
    data = data[sel, , drop=FALSE]
    data = as.matrix(data)
  }

  results = rep(0, length(poolSizes))

  #Normalize to TPM/CPM
  if(rescale) data <- tpmDSAVE(data)
  #Filter genes
  row_means <- rowMeans(data, na.rm = na.rm)
  sel = (row_means >= lowerBound) & (row_means <= upperBound)
  data <- data[sel, , drop=FALSE]

  #loop through all sizes
  for (i in 1:length(poolSizes)) {

    #Create combinations
    ind <- as.integer(1:dim(data)[2])
    ix <- 1
    numCells = dim(data)[2]
    combs.list <- list()
    if(!is.null(seed)) set.seed(seed)

    while(ix <= repetitionPerSize ) {
      #Don't use sample with a vector, it behaves differently when the vector is of size 1!
      ind.tmp <- ind
      a <- sample(numCells, poolSizes[i])
      ind.tmp <- ind[-a, drop=FALSE]
      b <- ind.tmp[sample(numCells-poolSizes[i], poolSizes[i]),drop=FALSE]
      combs.list[[ix]] <- list(a = a, b = b)
      ix <- ix + 1

    }
    #Calculate variation
    mean_vector <- sapply(combs.list, function(v){
      a <- rowMeans(data[,v$a,drop=FALSE], na.rm = na.rm)
      b <- rowMeans(data[,v$b,drop=FALSE], na.rm = na.rm)
      mean(abs(log((a + 0.05) / (b + 0.05))))
    })

    results[i] = mean(mean_vector);
    if (!silent) {
      pb$tick()
    }
  }

  if (!silent) {
    pb$terminate()
  }

  return(list(poolSizes = poolSizes, Rs = results))
}
