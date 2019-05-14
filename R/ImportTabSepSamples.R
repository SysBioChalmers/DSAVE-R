#' Text file to Sample object
#'
#' Imports a txt file and creates a Sample object
#'
#' Imports a txt tab separated, takes first column as rows' names and first row as columns' name
#'
#' @import utils
#' @export
#' @param filename path to the file
#' @param sampleName name of the Sample object
#' @author Juan Inda, <inda@@chalmers.se>
#' @examples
#' \dontrun{
#' samp <- ImportTabSepSamples("~/Downloads/dsaver/tcellCD4ProfilesTMMNormalized.txt",
#' "tcellCD4Profiles")
#' }
ImportTabSepSamples <- function(filename, sampleName){
  samp  <- as.matrix(read.table(filename, header = T, sep = "\t", row.names = 1))
  tcellCD4Profiles <- Samples$new(sampleName, samp)
}


