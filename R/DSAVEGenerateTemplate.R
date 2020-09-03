#' DSAVEGenerateTemplate
#'
#' Creates a template from an existing dataset and a gene list
#'
#' @param ds a dataset, numerical matrix with genes as row names
#' @param genesToUse (optional) List of genes to filter on in the template; defaults to the genes in ds
#' @param numCells (optional) Number of cells to include in the template. Defaults to 2000.
#' @param numUMIs (optional)  Average number of UMIs per cell to aim for in the template
#' Defaults to 750.
#' @param fractionUpperOutliers the fraction of upper outliers that should be discarded.
#' Defaults to 0.025.
#' @param fractionLowerOutliers the fraction of lower outliers that should be discarded.
#' Defaults to 0.025.
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return list<...>


DSAVEGenerateTemplate <- function(ds, genesToUse = NULL, numCells = 2000, numUMIs = 750,
                                  fractionUpperOutliers = 0.025, fractionLowerOutliers = 0.025){
  stopifnot((is.matrix(ds) | is(ds, 'sparseMatrix')),
            numCells == round(numCells), length(numCells) == 1,
            numUMIs == round(numUMIs), length(numUMIs) == 1,
            dim(ds)[2] >= numCells,
            is.numeric(fractionUpperOutliers), length(fractionUpperOutliers) == 1,
            is.numeric(fractionLowerOutliers), length(fractionLowerOutliers) == 1)

  if (is.null(genesToUse)) {
    genesToUse = row.names(ds);
  } else {
    #Remove the genes that are filtered out from the dataset
    genes <- intersect(row.names(ds),genesToUse)
    indices = match(genes, row.names(ds))
    ds = ds[indices, , drop=FALSE]

  }

  #select numCells cells:
  ds = ds[, sample(dim(ds)[2], numCells)]

  m = as.matrix(ds)

  result = list()
  result$binningInfo <- generateBinningInfo(m);
  result$UMIDistr <- generateUMIDistribution(m, numUMIs);
  result$geneSet <- genesToUse;
  result$fractionUpperOutliers <- fractionUpperOutliers
  result$fractionLowerOutliers <- fractionLowerOutliers


  return (result)
}
