#' GetVarAndLogCV
#'
#' GetVarAndLogCV
#'
#' GetVarAndLogCV
#'
#' @param data numeric matrix, the input dataset (cell population)
#' @param UMIsPerCell number of UMIs to normalize by.
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return list(logCV, variances)
#' @examples
#' \dontrun{
#' }

GetVarAndLogCV <- function(dsdata, UMIsPerCell){
  #tpm with specified UMIs per cell
  #There's no good way to do this without a loop in R
  for(i in 1:dim(dsdata)[2]){
    dsdata[,i] <- dsdata[,i] * 1e6 / UMIsPerCell[i]
  }
  avgRefExpr <- rowMeans(dsdata)
  variances <- apply(dsdata, 1, var)
  sds <- sqrt(variances)
  #Coefficient of Variation = std/mean. Adding 0.05, a neglectably small number,
  #to handle too lowly expressed genes
  cv <- sds / (avgRefExpr + 0.05)
  #the + 1 says that no variance -> 0 value
  logCV = log(cv + 1)
  return(list(logCV=logCV, variances = variances))
}
