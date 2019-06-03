#' Template object
#'
#' Samples object elements and methods
#'
#' Creates Samples object, updates its values and contains the functions of
#' this object
#'
#' @param name Name of the object
#' @param data Numeric matrix
#' @param genesToKeep genes to include in the template
#' @param numCells number of cells in the template
#' @param numUMI average number of molecules per cell
#' @param fractionUpperOutliers fraction of upper outliers
#' @param fractionLowerOutliers fraction of lower outliers
#' @param binningInfo.lbs distribution
#' @param binningInfo.ubs distribution
#' @param UMIDistr distribution of UMIs
#' @param template.built logical
#' @param lbnonlog lower bound in original scale
#' @param ubnonlog upper bound in original scale
#' @param numPoints number of points for the distribution
#' @param poolSize number of genes to pool at each point in the distribution
#' @param step increas in lower and upper bound to get enough genes at each point in the distribution
#' @param toplim maximum number of iterations to obtain the frequency at each point in the distribution
#' @import R6
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return a template object
#' @examples
#' \dontrun{
#' }


# templInfo <- R6Class("templInfo", inherit = Samples, list(
#   fractionUpperOutliers = NULL,
#   fractionLowerOutliers = NULL,
#   binningInfo.lbs = NULL,
#   binningInfo.ubs = NULL,
#   template.built = FALSE,
#   UMIDistr = NULL,
#
#   initialize = function(name, data = NULL, genes = NULL, sampleIds = NULL, ...) {
#     super$initialize(name, data = data, genes = genes, sampleIds = sampleIds, ...)
#   },
#
#   print = function(...){
#     cat("Template: \n")
#     cat("  name: ", self$name, "\n", sep = "")
#     cat("  data:  ", self$data[1,1:5], "... \n", sep = "")
#     if(!is.null(self$genesToKeep)) cat("  genes to keep:  ", self$genesToKeep[1:5], "...\n", sep = "")
#     cat("  number of cells in data:  ", dim(self$data)[2], "\n", sep = "")
#     #if(!is.null(self$numCells)) cat("  number of cells to form template:  ", self$numCells, "\n", sep = "")
#     #if(!is.null(self$numUMI))    cat("  number of UMIs to form template:  ", self$numUMI, "\n", sep = "")
#     if(!is.null(self$fractionLowerOutliers)) cat("  fractionLowerOutliers to form template:  ", self$fractionLowerOutliers, "\n", sep = "")
#     if(!is.null(self$fractionUpperOutliers)) cat("  fractionUpperOutliers to form template:  ", self$fractionUpperOutliers, "\n", sep = "")
#     if(!is.null(self$UMIDistr))    cat("  UMI distribution:  ", self$UMIDistr[1:5], "\n", sep = "")
#     if(!is.null(self$binningInfo.lbs))    cat("  UMI distribution:  ", self$binningInfo.lbs[1:5], "\n", sep = "")
#     if(!is.null(self$binningInfo.ubs))    cat("  UMI distribution:  ", self$binningInfo.ubs[1:5], "\n", sep = "")
#     invisible(self)
#   },
#
#   setBinningInfo = function(lbnonlog = NULL, ubnonlog = NULL,
#                             numPoints = 1000L, poolSize = 500L, step = 5e-4, toplim = 1000L, ...){
#     stopifnot(is.integer(lbnonlog), length(lbnonlog) == 1)
#     stopifnot(is.integer(ubnonlog), length(ubnonlog) == 1)
#     stopifnot(lbnonlog < ubnonlog)
#     stopifnot(is.integer(numPoints), length(numPoints) == 1)
#     stopifnot(is.integer(poolSize), length(poolSize) == 1)
#     stopifnot(is.integer(toplim), length(toplim) == 1)
#     stopifnot(is.numeric(step), length(step) == 1)
#
#     meansLog <- seq(log10(lbnonlog), log10(ubnonlog), length.out =numPoints)
#     templInfo.binningInfo.means <- 10^meansLog
#     cat("Re-scaling gene means \n")
#     tpmData <- tpmDSAVE(self$data)
#     #self$data <- tpmDSAVE(self$data)
#
#     gm <- as.matrix(rowMeans(tpmData), ncol = 1)
#     gmLog <- log10(gm)
#
#     lbs <- rep(0, numPoints)
#     ubs <- rep(0, numPoints)
#
#     cat("Forming the distributions \n")
#
#     for(i in 1:numPoints){
#       if(i %in% seq(floor(numPoints / 10), numPoints, by = floor(numPoints / 10))) {
#         cat(paste0("progress: ", round(i / numPoints, 2) * 100, "% \n"))
#       }
#       lb <- meansLog[i] - step
#       ub <- meansLog[i] + step
#
#       for(s in 1:toplim){
#         sel <- gmLog >= lb & gmLog <= ub
#         numGenes <- sum(sel)
#         if(numGenes >= poolSize) break
#
#         # increase bounds
#         if(numGenes == 0) { meanOfCaughtGenes <- 0} else { meanOfCaughtGenes <- mean(gmLog[sel])}
#         if(meanOfCaughtGenes > meansLog[i]){
#           lb <- lb - step
#         } else {
#           ub <- ub + step
#         }
#       }
#       lbs[i] <- 10^lb
#       ubs[i] <- 10^ub
#     }
#     self$binningInfo.lbs <- lbs
#     self$binningInfo.ubs <- ubs
#     self$template.built = TRUE
#     invisible(self)
#   },
#
#   subGenes = function(genesToKeep = NULL,...){
#     s2 <- self$geneSubset(genesToKeep)
#     s2 <- templInfo$new(name = paste(self$name,"gene subset"), data = s2$data)
#     return(s2)
#   },
#
#   subCells = function(numCells = NULL,...){
#     stopifnot(is.integer(numCells), length(numCells) == 1)
#     if(numCells > self$numberSamples){
#       stop("The number of desired cells is larger than the number of cells in the dataset")
#     } else {
#       samplesToKeep <- sort(sample(self$numberSamples, numCells, replace = F))
#       s2 <- self$sampleSubset(self$sampleIds[samplesToKeep])
#       s2 <- templInfo$new(name = paste(self$name,"sample subset"), data = s2$data)
#       return(s2)
#     }
#   },
#
#   setUMIDistr = function(numUMI = NULL, ...){
#     stopifnot(is.integer(numUMI), length(numUMI) == 1)
#     totalUMI <- sum(self$data)
#     sumUMI <- colSums(self$data)
#     totalTargetUMI <- self$numberSamples * numUMI
#     toRemove <- totalUMI - totalTargetUMI
#     UMIDistr <- sumUMI
#     if(toRemove < 0){
#       stop("Could not reach target UMIs per cell due to that the UMIs in the template dataset were to few")
#     } else {
#       if(toRemove > 0){
#         indToRem <- sample(totalUMI, toRemove)
#         edges <- c(0,cumsum(UMIDistr)+0.01)
#         subtr <- hist(indToRem, breaks = edges, plot = F)$counts
#         UMIDistr <- UMIDistr - subtr
#       }
#       self$UMIDistr <- sort(UMIDistr)
#     }
#   },
#
#   setFractionOutliers = function(fractionLowerOutliers = NULL, fractionUpperOutliers= NULL, ...){
#     stopifnot(is.numeric(fractionLowerOutliers), length(fractionLowerOutliers) == 1)
#     stopifnot(is.numeric(fractionUpperOutliers), length(fractionUpperOutliers) == 1)
#     self$fractionLowerOutliers <- fractionLowerOutliers
#     self$fractionUpperOutliers <- fractionUpperOutliers
#   }
#
# )
# )
