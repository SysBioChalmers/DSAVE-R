#' tpmDSAVE 
#' 
#' Standarizes by column making it sum to 1e6
#' 
#' Standarizes by column making it sum to 1e6
#' 
#' @param count A matrix.
#' @return A matrix
#' @author Juan Inda, <inda@@chalmers.se>
tpmDSAVE <- function(counts){
  tmp.tpm <- t(t(counts)/colSums(counts)) * 1e6
  tmp.tpm[is.na(tmp.tpm)] <- 0
  return(tmp.tpm)
}
