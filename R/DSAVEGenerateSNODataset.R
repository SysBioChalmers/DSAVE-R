#' DSAVEGenerateSNODataset
#'
#' Generates a Sampling Noise Only (SNO) dataset
#'
#' This function uses a template dataset to compute the probability p that
#' a molecule from a certain gene should be picked each time a new UMI is
#' found. This is calculated from the counts value, i.e. counts/sum of all
#' counts. This is based on a multinomial distribution, and it will select
#' exactly the same number of UMIs that are in the original set using the
#' probabilities for the genes from the mean of the dataset sent in.
#' Important that this is really UMI counts in this case, not TPM!
#' Random multiplicative noise will be added (0 == no noise)
#' noiseLevel should be 0 or greater. A standard random normal distributed
#' noise multiplied by noiseLevel will be multiplied to the probabilities.
#' templDSForProfile -
#'
#' @param templDs The input dataset (cell population), numeric matrix
#' @param numCells (optional) Can be used to specify the number of cells. Defaults to the number of cells in the input dataset.
#' @param noiseLevel (optional) The noise level to add; defaults to 0 (no noise)
#' @param templDSForProfile (optional) can be used if you want to generate data from a different cell type -
#' defaults to templDs - the genes in the datasets need to be synchronized
#' @importFrom graphics hist
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return a matrix
#' @examples
#' \dontrun{ DSAVEGenerateSNODataset(ds)
#' }

DSAVEGenerateSNODataset <- function(templDs, numCells=NULL, noiseLevel=0, templDSForProfile=NULL){
  stopifnot(is.numeric(templDs), is.matrix(templDs))
  stopifnot((is.integer(numCells) & length(numCells) == 1) | is.null(numCells))
  stopifnot(is.numeric(noiseLevel), length(noiseLevel) == 1, noiseLevel >= 0)
  if(is.null(templDSForProfile)){ templDSForProfile <- templDs}
  if(is.null(numCells)) {numCells <- dim(templDs)[2]}

  #create empty dataset
  ds <- templDs
  ds[,] <- 0

  #generate probabilities
  tmpTempl <- tpmDSAVE(templDSForProfile)
  meanRefExpr <- rowMeans(tmpTempl)
  #meanRefExpr(isnan(meanRefExpr)) = 0;
  if(sum(meanRefExpr) == 0){prob <- 0} else{
    prob = meanRefExpr / sum(meanRefExpr)
  }
  #generate a vector with the sum of probabilities up to this gene (including this gene)

  ###   probSum <- rep(0, length(prob))  #just create a vector of the right size
  ###   numGenes <- dim(tmpTempl)[1]
  ###    for(i in 2:numGenes){
  ###     probSum[i] <- probSum[i-1] + prob[i-1]
  ###   }
  probSum <- cumsum(prob)
  #add right edge; make sure it is not smaller than the previous due to roundoff problems.

  edges <- c(0, probSum)


  #create vector of number of UMIs
  if(dim(templDs)[2] == numCells){
    UMIs <- colSums(templDs)
  } else {
    UMIs <- rep(0, numCells)
    sumUMI <- colSums(templDs)
    ind <- 1
    cellsLeft <- numCells
    cellsInTemplate <- dim(templDs)[2]
    while(cellsLeft > 0){
      cellsThisRound <- min(cellsLeft, cellsInTemplate)
      sel <- sample(cellsInTemplate, cellsThisRound)
      UMIs[, ind:(ind+cellsThisRound-1)] <- sumUMI[sel]
      ind <- ind + cellsThisRound
      cellsLeft <- cellsLeft - cellsThisRound
    }
    UMIs <- round(UMIs) # in case the dataset has been scaled or something
  }

  #generate cells
  if(noiseLevel == 0){# no noise, will go faster
    for(i in 1:numCells){
      if(UMIs[i] > 0){#handle the fact that cells sometimes can contain 0 UMIs,
        #which doesn't work with the code below. The data is automatically 0 in that case.
        r <- runif(UMIs[i])
        ds[,i] <- hist(r, breaks = edges, plot = F)$counts
      }
    }
  } else{ # with noise
    numGenes <- dim(prob)[1]
    for( i in 1:numCells){
      X <- rnorm(numGenes) * noiseLevel
      noise <- 2^X
      probTmp <- prob * noise
      #scale probabilities back to the sum of 1
      probTmp <- probTmp / sum(probTmp)
      #have to recalculate edges
      probSum <- cumsum(prob)
      edges <- c(0, probSum)
      r <- runif(UMIs[i])
      ds[,i] <- hist(r, breaks = edges, plot = F)$counts
    }
  }
  return(ds)
}


