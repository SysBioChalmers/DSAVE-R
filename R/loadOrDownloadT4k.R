#' loadOrDownloadT4k
#'
#' Gets the t4k dataset from 10x Genomics. The function tries to save the dataset as an
#' .rds file in the package data folder to avoid having to download it every time it is
#' requested. If that folder is not writable for some reason, it will use the temp folder,
#' which is cleared when you exit R, meaning it will have to be downloaded again. In such cases, we
#' recommend that you save the object separately.
#' This function uses Seurat.
#'
#' @importFrom graphics hist
#' @importFrom utils untar download.file
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#' @return the template

loadOrDownloadT4k <- function() {
  #solve the problem that the current directory is sometimes the package
  #root and sometimes in the testthat folder
  #root = dirname(is_r_package$find_file("DSAVE_ROOT_IDENTIFIER.txt"))
  #packageRoot = paste0(path.package("DSAVE"),"/tempData");
  #packageRoot = "H:/"; #add this line to test that save to temp folder instead works
  #a = TryLoadT4k(packageRoot)
  #if (is.null(a)) {
    tmpdir = tempdir(check = TRUE)
    tmpdir = gsub('\\','/',tmpdir, fixed=TRUE)
    a = TryLoadT4k(tmpdir)
  #}

  return(a)
}

TryLoadT4k <- function(root) {

  RDAfilename = paste0(root,"/T4k.rds")
  gzfilename = paste0(root,"/T4k.tar.gz")
  extrDir = paste0(root,"/T4k")
  dataDir = paste0(root,"/T4k/filtered_gene_bc_matrices/GRCh38")
  data = NULL

  out <- tryCatch(
    {
      if (!file.exists(RDAfilename)) {
        #download, untar and load the t4k data and then store it in a .rda file
        download.file("http://cf.10xgenomics.com/samples/cell-exp/2.1.0/t_4k/t_4k_filtered_gene_bc_matrices.tar.gz",
                      gzfilename, "internal")
        untar(gzfilename, exdir = extrDir)
        data <- Read10X(data.dir = dataDir)
        saveRDS(data, file=RDAfilename)
        file.remove(gzfilename)
        unlink(extrDir, recursive = TRUE) #deletes the directory
      }
      else {
        data = readRDS(RDAfilename)
      }
    },
    error=function(cond) {
    },
    warning=function(cond) {
    },
    finally={
    }
  )

  return (data)

}

