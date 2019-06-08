#' DSAVEAlignDataset
#'
#' Align a data set to a template
#'
#' Aligns the dataset to the template, by removing genes, removing cells and down-sampling the data to match the template as good as possible.
#' All datasets successfully aligned to the same template will have almost
#' identical sampling noise.
#'
#' @param data matrix or sparse matrix, the type is preserved in the output
#' @param templInfo list with one element called UMIDistr containing template information
#' @importFrom graphics hist
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return a matrix
#' @examples
#' \dontrun{ DSAVEAlignDataset(data, templInfo)
#' }

DSAVEAlignDataset <- function(data, templInfo){
  stopifnot(is.matrix(data) | is(data, 'sparseMatrix'))
  stopifnot(length(templInfo) == 5, sum(!names(templInfo) %in%
                                          c("UMIDistr", "geneSet",
                                            "fractionUpperOutliers",
                                            "fractionLowerOutliers",
                                            "binningInfo")) == 0)
  # create the adapted dataset of the right size and with the right genes
  nCells <- dim(data)[2]
  stopifnot(nCells >=  length(templInfo$UMIDistr), nCells > 0)
  #stopifnot(!is.null(names(templInfo$UMIDistr)))
  stopifnot(!is.null(rownames(data)))

  if(is.null(colnames(data))){ colnames(data) <- 1:nCells}
  gd <- rownames(data) %in% templInfo$geneSet
  if(length(gd) == 0){stop("there are not ovelapping genes between template and dataset")}

  id <- sample(1:nCells, length(templInfo$UMIDistr))
  data_sub <- data[gd,id]
  ds <- data_sub

  #ds.name = ['Aligned dataset from ds ' inDs.name];
  numCells <- dim(data_sub)[2]
  ds[,] <- 0

  #unsparsify if needed
  wasSparse = is(data_sub, 'sparseMatrix')
  if (wasSparse) {
    data_sub = as.matrix(data_sub);
  }

  #then do downsampling
  origUMIs <- colSums(data_sub) #sum of UMIs for each cell
  matchUMIs <- templInfo$UMIDistr #sum of UMIs for each cell
  #we need create matching UMIs from the template. The template does not
  #necessarily have the same number of cells, so we need to handle that.
  #First duplicate the template until there are fewer cells left than in the
  #template; then randomly select from the template to fill out the rest.

  #sort both UMI vectors and try to downsample to match the match set
  #as good as possible
  iOrig <- order(origUMIs, decreasing = F)
  iMatchUMIs <- order(matchUMIs, decreasing = F)
  idealUMIs <- rep(0, numCells)
  idealUMIs[iOrig] <- matchUMIs[iMatchUMIs]

  tooSmall <- origUMIs < idealUMIs
  newUMIs <- pmin(origUMIs, idealUMIs)
  UMIsLost <- sum(idealUMIs[tooSmall] - newUMIs[tooSmall])
  UMIsToSpend <- UMIsLost
  toRemUMIs <- origUMIs - newUMIs

  #for some cells, the UMIs
  #spread the lost UMIs over the other cells
  i <- 1

  #take care of the fact that the match dataset may
  #have more UMIs than the ds to downsample
  if(UMIsToSpend > sum(toRemUMIs)){
    cat("Warning: Failed to downsample due to that the template had more UMIs than the target dataset! \n")
    toRemUMIs = rep(0, numCells)
  } else{
    while(UMIsToSpend > 0){
      if(toRemUMIs[i] > 0){
        toRemUMIs[i] <- toRemUMIs[i] - 1;
        UMIsToSpend <- UMIsToSpend - 1;
      }
      i <- i+1
      if(i > numCells){i <- 1}
    }
  }

  for(i in 1:numCells){
    #first, randomly select a number of indexes from the total UMIs to
    #remove
    indexesToRem <- sample(origUMIs[i], toRemUMIs[i], replace = FALSE)
    #Then create index edges for each gene, so if a gene has 5 UMIs the
    #edges to the left and right will differ 5.
    cellData <- data_sub[,i]
    edges <- c(0, cumsum(cellData) + 0.1)
    #we need to add 0.1 to the edges since histcounts checks if edge(k) <= index < edge(k+1). Otherwise, index 1 would not end up between 0 and 1
    #Now get the number of index hits you get within the edge range for
    #each gene
    subtr <- hist(indexesToRem, breaks = edges, plot = F)$counts
    ds[,i] <- data_sub[,i] - subtr
  }

  if (wasSparse) {
    ds = as(ds, "sparseMatrix");
  }

  return(ds)
}
