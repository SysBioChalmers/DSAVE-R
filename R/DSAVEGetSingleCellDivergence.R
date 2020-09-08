#' DSAVEGetSingleCellDivergence
#'
#' Calculates the DSAVE cell-wise variation metric.
#'
#' The divergence is the negative
#' log-likelihood for getting the observed counts' distribution for each
#' cell when sampling counts from the mean dataset gene expression. The
#' values are positive, and the higher the value, the
#' more divergent the cell is.
#'
#' @param data numeric matrix (can be sparse), the input dataset (cell population)
#' @param minUMIsPerCell All cells are downsampled to this number of counts, or.
#' higher if possible without losing cells. The result is NA for cells with lower counts.
#' @param numGenesGW The number of gene specific divergence values to store per cell. Will store the most divergent ones
#' @param iterations (optional) The number of times to iterate. Defaults to 15.
#' @param tpmLowerBound (optional) TPM lower bound, genes below this will not
#' be included. Defaults to 0 (meaning all genes are included).
#' @param silent (optional) If true, no progress bar is shown. Defaults to FALSE
#' @importFrom stats median
#' @importFrom robustbase colMedians
#' @importFrom combinat dmnom
#' @importFrom stats dmultinom
#' @importFrom progress progress_bar
#' @importFrom textTinyR sparse_Sums
#' @importFrom dynutils is_sparse
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return a list (divs = vector of divergence, one value per cell, geneDivGenes = The top divergent genes per cell, a matrix, geneDivVals = The corresponding divergence values for the divergent genes per cell)
#' @examples
#' \dontrun{a = DSAVEGetSingleCellDivergence(data)}
DSAVEGetSingleCellDivergence <- function(data, minUMIsPerCell = 200, numGenesGW=5, tpmLowerBound = 0, iterations = 15, silent=FALSE) {
  stopifnot(is.matrix(data) | is(data, 'sparseMatrix'))
  stopifnot(minUMIsPerCell == round(minUMIsPerCell), length(minUMIsPerCell) == 1)
  stopifnot(is.numeric(tpmLowerBound),  length(tpmLowerBound)==1)
  stopifnot(iterations == round(iterations), length(iterations) == 1)

  #data = as.matrix(data);
  numCells = dim(data)[2];

  if(dynutils::is_sparse(data)) {
    origUMIs = as.numeric(textTinyR::sparse_Sums(data, rowSums=F))
    #half of TPM summation calculation
    #totExpr = as.numeric(textTinyR::sparse_Sums(data, rowSums=T))
  } else {
    origUMIs = as.numeric(colSums(data))
    #totExpr = as.numeric(rowSums(data))
  }

  #calc CPM
  totExpr = rep(0,dim(data)[1])
  for(i in 1:numCells) {
    totExpr = totExpr + data[,i]/origUMIs[i]
  }
  expr = totExpr*10^6/sum(totExpr);
  #remove all lowly expressed genes (if that is specified)
  sel = (expr >= tpmLowerBound) & (expr != 0)
  data = data[sel, , drop=FALSE];
  expr = expr[sel]

  numGenes = dim(data)[1];

  if (!silent) {
    pb <- progress_bar$new(format = "Calculating cell divergence [:bar] :percent eta: :eta",
                           total = numCells + 1, clear = FALSE)
    pb$tick();
  }

  #Figure out what to downsample to. If we have cells with fewer UMIs than
  #minUMIsPerCell, those can not be evaluated
  UMIsPerCell = origUMIs;
  targetUMIsPerCell = max(c(min(origUMIs),minUMIsPerCell));

  #Get mean expression and create probabilities

  #meanRefExpr = apply(expr,1,mean);
  prob = as.numeric(expr / sum(expr));#don't divide by 10^6, since the expression can be filtered

  lls = rep(NA, numCells)#matrix(0,iterations, numCells);
  iterLls = rep(NA, iterations)

  iterGeneLls = matrix(NA, nrow = numGenes, ncol = iterations)


  logFacs = GetLogFactorials(targetUMIsPerCell);

  toRemUMIs = UMIsPerCell - targetUMIsPerCell;
  toRemUMIs[toRemUMIs < 0] = 0;

  scDivGenes = matrix("", nrow=numGenesGW, ncol=numCells)
  scDivVals = matrix(NA, nrow=numGenesGW, ncol=numCells)

  for (i in 1:numCells) {
    #run several times since there is randomness in downsampling
    for (it in 1:iterations) {
      dsCellData = data[,i]
      if (toRemUMIs[i] > 0) {
        #first, randomly select a number of indexes from the total UMIs to
        #remove
        indexesToRem = sample(origUMIs[i], toRemUMIs[i], replace = FALSE);
        #Then create index edges for each gene, so if a gene has 5 UMIs the
        #edges to the left and right will differ 5.
        edges = c(0,cumsum(dsCellData)+0.1);#we need to add 0.1 to the edges since histcounts checks if edge(k) <= index < edge(k+1). Otherwise, index 1 would not end up between 0 and 1
        #Now get the number of index hits you get within the edge range for
        #each gene
        subtr <- hist(indexesToRem, breaks = edges, plot = F)$counts
        dsCellData <- dsCellData - subtr
      }
      if (UMIsPerCell[i] < targetUMIsPerCell) {
        iterLls[it] = NA;
      } else {
        #calculate log likelihood from the
        #multinomial distribution
        iterLls[it] = dmultinom(dsCellData, sum(dsCellData), prob = prob, log = TRUE)
        #Also fill in gene binomial pdf
        iterGeneLls[,it] = LogBinomialPDF(dsCellData, prob = prob, logFacs);
      }
    }
    lls[i] = median(iterLls)
    iterGeneLls[is.nan(prob),] = 0;
    tmp = matrixStats::rowMedians(iterGeneLls)
    sel = tmp != 0
    tmp = tmp[sel]
    genes = row.names(data)[sel]
    srt = sort(tmp, index.return=T)
    scDivGenes[,i] = genes[srt$ix[1:numGenesGW]]
    scDivVals[,i] = srt$x[1:numGenesGW]
    if (!silent) {
         pb$tick()
    }
  }

  if (!silent) {
    pb$terminate()
  }

  return(list(divs = -lls, geneDivGenes = scDivGenes, geneDivVals = -scDivVals))
}
