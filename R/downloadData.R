#' downloadData
#'
#' Downloads a dataset. The function tries to save the dataset as an
#' .rds file in the package data folder to avoid having to download it every time it is
#' requested. If that folder is not writable for some reason, it will use the temp folder,
#' which is cleared when you exit R, meaning it will have to be be downloaded again. In such cases, we
#' recommend that you save the object separately.
#' This function uses Seurat.
#'
#' @param url The URL to download
#' @param filenameBase The dataset name to use as filename
#' @param locDir The directory in which the data can be stored on your disk
#' @importFrom graphics hist
#' @importFrom Seurat Read10X
#' @importFrom downloader download
#' @importFrom utils untar download.file
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return the template

downloadData <- function(url, filenameBase, locDir = tempdir(check = TRUE)) {
  locDir = gsub('\\','/',locDir, fixed=TRUE)
  base = paste0(locDir, "/", filenameBase)

  RDAfilename = paste0(base,".rds")
  gzfilename = paste0(base,".tar.gz")
  extrDir = base
  BCellsData = NULL

  out <- tryCatch(
    {
      #download, untar and load the bcell data and then store it in a .rda file
      downloader::download(url, destfile = gzfilename)
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

