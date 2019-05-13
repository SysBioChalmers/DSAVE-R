#' Normalizes all the columns to a 1e6 count.
#' 
#' @param count A matrix.
tpmDSAVE <- function(counts){
  tmp.tpm <- t(t(counts)/colSums(counts)) * 1e6
  tmp.tpm[is.na(tmp.tpm)] <- 0
  return(tmp.tpm)
}
