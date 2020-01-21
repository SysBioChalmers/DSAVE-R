#' extractSeuratData
#'
#' Extracts the data from a Seurat object. Optionally, a cluster id can be sent in, to extract
#' a subset of the data
#'
#' @importFrom graphics hist
#' @importFrom Seurat Read10X
#' @importFrom utils untar download.file
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return the template

extractSeuratData <- function(seuratObj, clusterId = NA) {

  if (is.na(clusterId)) {
    data = so@assays$RNA@counts
  } else {
    data = so@assays$RNA@counts[,so@active.ident == clusterId]
  }

  return (data)
}

