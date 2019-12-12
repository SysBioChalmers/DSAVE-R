#' downloadData
#'
#' Gets the b10k dataset from 10x Genomics. The function tries to save the dataset as an
#' .rds file in the package data folder to avoid having to download it every time it is
#' requested. If that folder is not writable for some reason, it will use the temp folder,
#' which is cleared when you exit R, meaning it will have to be be downloaded again. In such cases, we
#' recommend that you save the object separately.
#' This function uses Seurat.
#'
#' @importFrom graphics hist
#' @importFrom Seurat Read10X
#' @importFrom utils untar download.file
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return the template

downloadData <- function(url, filenameBase) {
  tmpdir = tempdir(check = TRUE)
  tmpdir = gsub('\\','/',tmpdir, fixed=TRUE)
  base = paste0(tmpdir, "/", filenameBase)

  RDAfilename = paste0(base,".rds")
  gzfilename = paste0(base,".tar.gz")
  extrDir = base
  BCellsData = NULL

  out <- tryCatch(
    {
      #download, untar and load the bcell data and then store it in a .rda file
      download.file(url, gzfilename, "internal")
      untar(gzfilename, exdir = extrDir)
      file.remove(gzfilename)
    },
    error=function(cond) {
    },
    warning=function(cond) {
    },
    finally={
    }
  )

  return (extrDir)
}

