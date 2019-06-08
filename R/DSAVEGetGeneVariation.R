#' DSAVEGetGeneVariation
#'
#' Calculate Gene Variation
#'
#' Calculates the DSAVE gene-wise BTM variation metric.
#'
#' @param data numeric matrix, the input dataset (cell population)
#' @param lb (optional) TPM lower bound, genes below this will not be investigated. Defaults to 10 TPM/CPM.
#' @param iterations (optional) The number of SNO datasets to generate.
#' Recommended value is 100 if no p values are needed,
#' 10,000 - 100,000 if p-values are of interest. Defaults to 100,000.
#' @param maxNumCells (optional) ds is reduced to this number of cells if it
#' contains more, to save computation time. Defaults to 2,000.
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return list(genes, logCVDifference, pVals, SNOVariances, SNOCountsPerGene)
#' @examples
#' \dontrun{templ <- DSAVEGetGeneVariation(ds, 1, 10000, 1000)
#' }
DSAVEGetGeneVariation <- function(data, lb=NULL, iterations = NULL, maxNumCells=NULL){
  stopifnot(is.numeric(data), is.matrix(data))
  stopifnot(iterations == round(iterations), length(iterations) == 1)
  stopifnot(is.numeric(lb),  length(lb)==1)
  stopifnot(maxNumCells == round(maxNumCells), length(maxNumCells) == 1)

  #cut down to 2000 cells to get somewhat reasonable computation times
  if(is.null(maxNumCells)){maxNumCells <- 2000 }
  if(is.null(iterations)){iterations <- 100}
  if(is.null(lb)){iterations <- 10}

  numCells <- dim(data)[2]
  if(numCells > maxNumCells){
    id <- sample(numCells, maxNumCells)
    data <- data[,id]
    #keep this line if we start using numCells below later
    numCells <- maxNumCells
  }


  data_tpm <- tpmDSAVE(data)
  dstpm <- rowMeans(data_tpm)
  #### dstpm <- tpmDSAVE(rowMeans(data))

  #skip all below threshold, and all with 0 or 1 counts, doesn't make sense to look at those
  sel <- (dstpm >= lb) & (dstpm != 0) & (rowSums(data) != 1)

  #filter genes to avoid problems with division by zero, etc, and to save compilation time.
  #No point in looking at too lowly expressed ones anyway, they will not become significant.
  stopifnot(sum(sel) > 0)

  data <- data[sel,,drop=FALSE]
  numGenes <- dim(data)[1]

  #generate SNO TPMs
  #convert the TPMs into counts
  SNOUMIsPerCell <- colSums(data)
  SNOCountsPerGene <- unique(rowSums(data))
  #skip 0 and 1
  id <- SNOCountsPerGene>=2
  stopifnot(sum(id)>0)
  SNOCountsPerGene <- SNOCountsPerGene[id]

  numCountVals <- length(SNOCountsPerGene)

  SNOLogCVS <- matrix(0, nrow=numCountVals, ncol=iterations)
  SNOVariances <- SNOLogCVS
  VarAndLog <- GetVarAndLogCV(data, colSums(data))
  logCVDS <- VarAndLog[["logCV"]]
  varianceDS <- VarAndLog[["variances"]]

  #precalc things for generating SNO datasets
  #generate probabilities
  prob <- SNOUMIsPerCell/sum(SNOUMIsPerCell)

  #generate a vector with the sum of probabilities up to this cell (including this cell)
  probSum <-cumsum(prob) #just create a vector of the right size


  #add right edge; make sure it is not smaller than the previous due to roundoff problems.
  SNOEdges <- c(0, probSum)

  #preallocate data matrix that will be reused and overwritten in each iteration,
  #to avoid copying of data
  SNOdata <- matrix(0, nrow = numCountVals, ncol = numCells)

  for( it in 1:iterations){
    SNOdata <- GenSampDs(SNOdata, SNOCountsPerGene, SNOEdges)
    #so, use the UMIs per cell from the original dataset when TPM:ing
    varAndLogCV <- GetVarAndLogCV(SNOdata, SNOUMIsPerCell)
    SNOLogCVS[,it] <- varAndLogCV[["logCV"]]
    SNOVariances[,it] <- varAndLogCV[["variances"]]
  }



  genes <-  rownames(data)

  #Calculate p value with a non-parametric method:

  countsPerGene <- rowSums(data)
  SNOCountsPerGene <- SNOCountsPerGene

  #map counts to SNO counts
  #should only find one index per row

  indices <- match(countsPerGene, SNOCountsPerGene)
  #indices <- arrayfun( @(x)( find(SNOCountsPerGene==x) ), countsPerGene);

  sz <- dim(SNOVariances)[2]

  pVals <- rowSums(SNOVariances[indices,] >= varianceDS) / sz


  #so subtraction of CVs
  logCVSNOm <- rowMeans(SNOLogCVS)

  logCVDifference <- logCVDS - logCVSNOm[indices]

  return(list(genes = genes, logCVDifference = logCVDifference,
              pVals = pVals, SNOVariances = SNOVariances,
              SNOCountsPerGene = SNOCountsPerGene))
}
