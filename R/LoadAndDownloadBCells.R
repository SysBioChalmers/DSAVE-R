#' LoadAndDownloadBCells
#'
#' Gets the standard template that is redistributed with the package. This function uses Seurat.
#' You need to run the following two lines to use it:
#' install.packages('Seurat')
#' library(Seurat)
#'
#' @importFrom graphics hist
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return the template
#' @examples
#' \dontrun{ LoadAndDownloadBCells()
#' }

LoadAndDownloadBCells <- function() {

  #download, untar and load the bcell data and then store it in a .rda file
  if (!file.exists("data/BCells.rda")) {
    download.file("http://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz",
                  "data/BCells.tar.gz", "internal")
    untar("data/BCells.tar.gz", exdir = "data/BCells")
    BCellsData <- Read10X(data.dir = "data/BCells/filtered_matrices_mex/hg19")
    save(BCellsData, file="data/BCells.rda")
    file.remove("data/BCells.tar.gz")
    unlink("data/BCells", recursive = TRUE) #deletes the directory
  }
  else {
    load("data/BCells.rda")
  }


  return(BCellsData)
}


