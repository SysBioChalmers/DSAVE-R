#' LoadAndDownloadBCells
#'
#' Gets the b10k dataset from 10x Genomics. The function tries to save the dataset as an
#' .rds file in the package data folder to avoid having to download it every time it is
#' requested. If that folder is not writable for some reason, it will use the temp folder,
#' which is cleared when you exit R, meaning it will be downloaded again. In such cases, we
#' recommend that you save the object separately.
#' This function uses Seurat.
#'
#' @importFrom graphics hist
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return the template
#' @examples
#' \dontrun{ LoadAndDownloadBCells()
#' }

LoadAndDownloadBCells <- function() {
  #solve the problem that the current directory is sometimes the package
  #root and sometimes in the testthat folder
  #root = dirname(is_r_package$find_file("DSAVE_ROOT_IDENTIFIER.txt"))
  packageRoot = paste0(path.package("DSAVE"),"/data");
  #packageRoot = "H:/"; #add this line to test that save to temp folder instead works
  a = TryLoadBCells(packageRoot)
  if (is.null(a)) {
    tmpdir = tempdir(check = TRUE)
    tmpdir = gsub('\\','/',tmpdir, fixed=TRUE)
    a = TryLoadBCells(tmpdir)
  }

  return(a)
}

TryLoadBCells <- function(root) {

  RDAfilename = paste0(root,"/BCells.rds")
  gzfilename = paste0(root,"/BCells.tar.gz")
  extrDir = paste0(root,"/BCells")
  dataDir = paste0(root,"/BCells/filtered_matrices_mex/hg19")
  BCellsData = NULL

  out <- tryCatch(
    {
      if (!file.exists(RDAfilename)) {
        #download, untar and load the bcell data and then store it in a .rda file
        download.file("http://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz",
                      gzfilename, "internal")
        untar(gzfilename, exdir = extrDir)
        BCellsData <- Read10X(data.dir = dataDir)
        saveRDS(BCellsData, file=RDAfilename)
        file.remove(gzfilename)
        unlink(extrDir, recursive = TRUE) #deletes the directory
      }
      else {
        BCellsData = readRDS(RDAfilename)
      }
    },
    error=function(cond) {
    },
    warning=function(cond) {
    },
    finally={
    }
  )

  return (BCellsData)

}

