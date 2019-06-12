#' DSAVEGetTotalVariationPoolSize
#'
#' Calculates the total variation with a specific pool size.
#'
#' Calculates the average pairwise total variation between two pools of cells with specific size,
#' with a TPM filtration on the genes.
#'
#' @param data a numeric matrix
#' @param poolSize numeric value, it should be maximum half of the total size of the data
#' in the data object
#' @param upperBound filters out the genes with mean expression higher than this value
#' @param lowerBound filters out the genes with mean expression lower than this value
#' @param na.rm logical, sends na.rm value to the rowMeans function
#' @param rescale logical, determine if the data should be rescaled
#' @param seed positive integer for to set the seed of the datas
#' @param repetitionPerSize positive integer how many times the variation is calculated at each
#' pool size
#' @importFrom methods is as
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return numeric vector
#' @examples
#' DSAVEGetTotalVariationPoolSize(data, poolSize = 50, upperBound = 1000,
#' lowerBound = 0.01, na.rm = TRUE, seed = 1, repetitionPerSize = 30, 0.5)
#'

DSAVEGetTotalVariationPoolSize <- function(data, poolSize = 4, upperBound = 1e5,
                                           lowerBound = 5e-1, na.rm = TRUE, seed = NULL,
                                           repetitionPerSize = 30L, rescale = TRUE){
  #print("Control parameters")
  stopifnot((is.matrix(data) |  is(data, 'sparseMatrix')),
            is.numeric(poolSize), length(poolSize) == 1,
            is.numeric(upperBound), length(upperBound) == 1,
            is.numeric(lowerBound), length(lowerBound) == 1,
            is.logical(na.rm), length(na.rm) == 1,
            is.logical(rescale), length(rescale) == 1,
            (is.null(seed) | is.integer(seed)),
            is.integer(repetitionPerSize), repetitionPerSize > 0)
  data <- as.matrix(data)

  if(poolSize > floor(dim(data)[2]/2)){
    stop("Cannot make non-overlapping pool of cells. Choose a smaller poolSize.")
  } else {
    #print("Re-scaling")
    if(rescale) data <- tpmDSAVE(data)
    #print("Filtering genes")
    row_means <- rowMeans(data, na.rm = na.rm)
    genesToKeep <- names(which(row_means >= lowerBound & row_means <= upperBound))
    data <- data[genesToKeep, , drop=FALSE]

    #print("Creating combinations")
    ind <- as.integer(1:dim(data)[2])
    ix <- 1
    numCells = dim(data)[2]
    combs.list <- list()
    if(!is.null(seed)) set.seed(seed)

    while(ix <= repetitionPerSize ) {
      #Don't use sample with a vector, it behaves differently when the vector is of size 1!
      ind.tmp <- ind
      a <- sample(numCells, poolSize)
      ind.tmp <- ind[-a, drop=FALSE]
      b <- ind.tmp[sample(numCells-poolSize, poolSize),drop=FALSE]
      combs.list[[ix]] <- list(a = a, b = b)
      ix <- ix + 1

    }
    #print("Calculating variation")
    mean_vector <- sapply(combs.list, function(v){
      a <- rowMeans(data[,v$a,drop=FALSE], na.rm = na.rm)
      b <- rowMeans(data[,v$b,drop=FALSE], na.rm = na.rm)
      mean(abs(log((a + 0.05) / (b + 0.05))))
    })

  }
  return(mean(mean_vector))
}
