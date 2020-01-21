#' extractSeuratData
#'
#' Extracts the data from a Seurat object. Optionally, a cluster id can be sent in, to extract
#' a subset of the data
#'
#' @param seuratObj The Seurat object
#' @param clusterId [optional] The id of the cluster, for example 6
#' @importFrom graphics hist
#' @importFrom Seurat Read10X
#' @importFrom utils untar download.file
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return the template

extractSeuratData <- function(seuratObj, clusterId = NA) {

  if (is.na(clusterId)) {
    data = seuratObj@assays$RNA@counts
  } else {
    data = seuratObj@assays$RNA@counts[,seuratObj@active.ident == clusterId]
  }

  return (data)
}

