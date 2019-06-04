#' tpmDSAVE
#'
#' Standarizes by column making it sum to 1e6
#'
#' Standarizes by column making it sum to 1e6
#'
#' @param counts A matrix.
#' @return A matrix
#' @author Juan Inda, <inda@@chalmers.se>
tpmDSAVE <- function(counts){
  counts2 <- counts
  for(i in 1:dim(counts)[2]){
    #s <- seq(from=1, to=dim(counts)[2], by = round(dim(counts)[2]/10))
    #if(i %in% s) print(paste("Doing: ", round(i/dim(counts)[2] * 100), "%"))
    nc <- counts[,i] / sum(counts[,i])
    nc[is.na(nc)] <- 0
    counts2[,i] <- nc * 1e6
  }
  return(counts2)
  #if(dim(counts)[2] > 1000)
  #tmp.tpm <- t(t(counts)/colSums(counts)) * 1e6
  #tmp.tpm[is.na(tmp.tpm)] <- 0
  #return(tmp.tpm)
}
