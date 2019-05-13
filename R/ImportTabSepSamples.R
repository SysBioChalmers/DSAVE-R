#samp <- ImportTabSepSamples("~/Downloads/dsaver/tcellCD4ProfilesTMMNormalized.txt", "tcellCD4Profiles")

ImportTabSepSamples <- function(filename, sampleName){
  samp  <- as.matrix(read.table(filename, header = T, sep = "\t", row.names = 1))
  tcellCD4Profiles <- Samples$new(sampleName, samp)
}


