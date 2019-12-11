#' DSAVEGetSingleCellDivergence
#'
#' Calculates the DSAVE cell-wise variation metric.
#'
#' The divergence is the
#' log-likelihood for getting the observed counts' distribution for each
#' cell when sampling counts from the mean dataset gene expression. The
#' values are negative, and the lower (i.e. more negative) the value, the
#' more divergent the cell is.
#'
#' @param data numeric matrix (can be sparse), the input dataset (cell population)
#' @param minUMIsPerCell All cells are downsampled to this number of counts, or.
#' higher if possible without losing cells. The result is NA for cells with lower counts.
#' @param iterations (optional) The number of times to iterate. Defaults to 15.
#' @param tpmLowerBound (optional) TPM lower bound, genes below this will not
#' be included. Defaults to 0 (meaning all genes are included).
#' @param silent (optional) If true, no progress bar is shown. Defaults to FALSE
#' @importFrom stats median
#' @importFrom combinat dmnom
#' @importFrom stats dmultinom
#' @importFrom progress progress_bar
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return a list (lls = vector of divergence, one value per cell, geneLls = )
#' @examples
#' \dontrun{a = DSAVEGetSingleCellDivergence(data)}
DSAVEGetSingleCellDivergence <- function(data, minUMIsPerCell = 200, tpmLowerBound = 0, iterations = 15, silent=FALSE) {
  stopifnot(is.matrix(data) | is(data, 'sparseMatrix'))
  stopifnot(minUMIsPerCell == round(minUMIsPerCell), length(minUMIsPerCell) == 1)
  stopifnot(is.numeric(tpmLowerBound),  length(tpmLowerBound)==1)
  stopifnot(iterations == round(iterations), length(iterations) == 1)

  if (!silent) {
    pb <- progress_bar$new(format = "Calculating cell divergence [:bar] :percent eta: :eta",
                           total = iterations + 1, clear = FALSE)
    pb$tick();
  }


  data = as.matrix(data);

  #remove all lowly expressed genes
  me = tpmDSAVE(matrix(rowMeans(data), ncol = 1));
  data = data[(me >= tpmLowerBound) & (me != 0), , drop=FALSE];

  numCells = dim(data)[2];
  numGenes = dim(data)[1];

  #Figure out what to downsample to. If we have cells with fewer UMIs than
  #minUMIsPerCell, those can not be evaluated
  origUMIs = as.numeric(colSums(data));
  UMIsPerCell = origUMIs;
  targetUMIsPerCell = max(c(min(origUMIs),minUMIsPerCell));

  #Get mean expression and create probabilities
  expr = tpmDSAVE(data);
  meanRefExpr = apply(expr,1,mean);
  prob = as.numeric(meanRefExpr / sum(meanRefExpr));

  allLls = matrix(0,iterations, numCells);

  allGeneLls = array(0, dim = c(iterations, numGenes, numCells));


  logFacs = GetLogFactorials(targetUMIsPerCell);

  #run several times since there is randomness in downsampling
  for (it in 1:iterations) {
    #downsample to the right number of UMIs per cell
    dsd = data;
    toRemUMIs = UMIsPerCell - targetUMIsPerCell;
    toRemUMIs[toRemUMIs < 0] = 0;
    for (i in 1:numCells) {
      if (toRemUMIs[i] > 0) {
        #first, randomly select a number of indexes from the total UMIs to
        #remove
        indexesToRem = sample(origUMIs[i], toRemUMIs[i], replace = FALSE);
        #Then create index edges for each gene, so if a gene has 5 UMIs the
        #edges to the left and right will differ 5.
        cellData = dsd[,i];
        edges = c(0,cumsum(cellData)+0.1);#we need to add 0.1 to the edges since histcounts checks if edge(k) <= index < edge(k+1). Otherwise, index 1 would not end up between 0 and 1
        #Now get the number of index hits you get within the edge range for
        #each gene
        subtr <- hist(indexesToRem, breaks = edges, plot = F)$counts
        dsd[,i] <- dsd[,i] - subtr
      }
    }


    #loop through all cells and calculate log likelihood
    for (i in 1:numCells) {
      if (UMIsPerCell[i] < targetUMIsPerCell) {
        allLls[it,i] = NaN;
        allGeneLls[it,,i] = NaN;
      } else {
        #calculate log likelihood from the
        #multinomial distribution
        allLls[it,i] = dmultinom(dsd[,i], sum(dsd[,i]), prob = prob, log = TRUE)
        #Also fill in gene binomial pdf
        allGeneLls[it,,i] = LogBinomialPDF(dsd[,i], prob = prob, logFacs);
      }
    }
    if (!silent) {
      pb$tick()
    }
  }


  #take the median of all runs
  lls = apply(allLls,2,median);

  #replace all nan in the gene-wise with 0. So, nan means that the gene is not
  #expressed, and then we are certain of the outcome -> PDF = 1
  # ->log(PDF) = 0
  allGeneLls[,is.nan(prob),] = 0;
  geneLls = t(apply(allGeneLls, c(3,2), median))

  if (!silent) {
    pb$terminate()
  }

  return(list(lls = lls, geneLls = geneLls))
}
