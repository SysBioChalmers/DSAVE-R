#' tpmDSAVE
#'
#' Standarizes by column making it sum to 1e6
#'
#' Standarizes by column making it sum to 1e6
#'
#' @param counts A matrix.
#' @return A matrix
#' @author Juan Inda, <inda@@chalmers.se>, Johan Gustafsson, <gustajo@@chalmers.se>
tpmDSAVE <- function(counts){
  for(i in 1:dim(counts)[2]){
    counts[,i] <- counts[,i] * 1e6 / sum(counts[,i])
  }
  counts[is.na(counts)] <- 0
  return(counts)
}
