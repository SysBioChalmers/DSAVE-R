#' DSAVECalcBTMScore
#'
#' Calculate BTM Score
#'
#' Calculate BTM Score
#'
#' @param data numeric matrix (can be sparse), the input dataset (cell population)
#' @param templInfo template information
#' @param skipAlignment logical, default FALSE. If FALSE, alignment will be done against the template data
#' @param iterations number of iterations; defaults to 15
#' @param useLogTransform default FALSE, if TRUE calculation will be done on log transformed data
#' @param logTPMAddon default 1, when using log transformed data, this number is added to the data before transforming to avoid log(0)
#' @param silent (optional) If true, no progress bar is shown. Defaults to FALSE
#' @importFrom stats approxfun
#' @export
#' @author Juan Inda, <inda@@chalmers.se>
#' @return list
#' @examples
#' \dontrun{
#' }

DSAVECalcBTMScore <- function(data, templInfo, skipAlignment=FALSE, iterations = 15,
                              useLogTransform=FALSE, logTPMAddon=1, silent=FALSE){

  stopifnot(is.matrix(data) | is(data, 'sparseMatrix'))
  stopifnot(iterations == round(iterations), length(iterations) == 1)
  stopifnot(is.numeric(logTPMAddon) & length(logTPMAddon) == 1)
  stopifnot(is.logical(skipAlignment), length(skipAlignment) == 1)
  stopifnot(is.logical(useLogTransform), length(useLogTransform) == 1)
  stopifnot(length(templInfo) == 5, sum(!names(templInfo) %in%
                                          c("UMIDistr", "geneSet",
                                            "fractionUpperOutliers",
                                            "fractionLowerOutliers",
                                            "binningInfo")) == 0)


  lbnonlog <- 10
  ubnonlog <- 1000
  numPoints <- 1000 # use the same number of points to simplify for data allocation
  meansLog <- seq(log10(lbnonlog), log10(ubnonlog), length.out = numPoints)
  xes <-  10^meansLog

  numXVals <- length(xes)
  alignedCVs <- matrix(0, iterations, numXVals) #zeros(iterations,numXVals)
  samplingCVs <- alignedCVs #zeros(iterations,numXVals);
  differenceCVs <- alignedCVs #zeros(iterations,numXVals);

  if (!silent) {
    pb <- progress_bar$new(format = "Calculating DSAVE score [:bar] :percent eta: :eta",
                           total = iterations + 1, clear = FALSE)
    pb$tick();
  }

  for(it in 1:iterations){
    if(skipAlignment){
      aligned <- data
    } else {
      aligned <- DSAVEAlignDataset(data, templInfo) #%use a different subset of cells in each loop iteration
    }

    aligned = as.matrix(aligned); # unsparse if needed

    SNO <- DSAVEGenerateSNODataset(aligned)

    alignedGeneCVs <- GetGeneCVs(aligned, logTPMAddon, toLog=useLogTransform);
    SNOGeneCVs <- GetGeneCVs(SNO, logTPMAddon, toLog=useLogTransform)


    #throw away the most and least variable genes
    numGenes <- dim(aligned)[1]
    numToDiscardUpper <- round(templInfo$fractionUpperOutliers * numGenes)
    numToDiscardLower <- round(templInfo$fractionLowerOutliers * numGenes)
    difference <- alignedGeneCVs - SNOGeneCVs
    gi <- order(difference, decreasing = FALSE)
    discard <- unique(c(head(gi,numToDiscardLower), tail(gi,numToDiscardUpper)))
    anyToDiscard <- sum(discard) > 0

    alignedGeneCVsRem <- alignedGeneCVs
    SNOGeneCVsRem <- SNOGeneCVs
    if(anyToDiscard){
      alignedGeneCVsRem <- alignedGeneCVsRem[-discard]
      SNOGeneCVsRem <- SNOGeneCVsRem[-discard]
    }

    alData <- tpmDSAVE(aligned)
    SNOData <- tpmDSAVE(SNO)
    if(length(discard)>0){
      alData <- alData[-discard,]
      SNOData <- SNOData[-discard,]
    }

    alCVs <- AverageIntoBins(alData, alignedGeneCVsRem, templInfo)
    alXes <- alCVs[[2]]
    alCVs <- alCVs[[1]]

    saCVs <- AverageIntoBins(SNOData, SNOGeneCVsRem, templInfo)
    saXes <- saCVs[[2]]
    saCVs <- saCVs[[1]]

    #remove any NaNs (let linear interpolation fill in)
    #now, recalculate the x:es and cvs to the Xes in the template using linear interpolation;
    #otherwise it will be difficult to take the mean later

    #first remove any points with duplicate x; these will
    #otherwise mess up the linear interpolation.
    ia <- !duplicated(alXes)
    alXes <- alXes[ia]
    alCVs = alCVs[ia]
    alXes <- alXes[!is.na(alXes)]
    alCVs <- alCVs[!is.na(alCVs)]

    ia <- !duplicated(saXes)
    saXes <- saXes[ia]
    saCVs <- saCVs[ia]
    saXes <- saXes[!is.na(saXes)]
    saCVs <- saCVs[!is.na(saCVs)]

    #then use linear interpolation
    #if the code fails on either of these two lines, that is because there are no genes that fits some bins
    #in that case, try using a template that discards fewer outliers
    alignedCVs[it,] <- approx(alXes,alCVs,xes, rule = 2)$y
    samplingCVs[it,] <- approx(saXes,saCVs,xes, rule = 2)$y
    differenceCVs[it,] <- alignedCVs[it,] - samplingCVs[it,]

    if (!silent) {
      pb$tick()
    }

  }
  if (!silent) {
    pb$terminate()
  }

  results = list();
  results[["alignedCVs"]] <- colMeans(alignedCVs)
  results[["samplingCVs"]] <- colMeans(samplingCVs)
  results[["differenceCVs"]] <- colMeans(differenceCVs)
  #mean over all points ranging over different TPM
  results[["DSAVEScore"]] <- mean(results[["differenceCVs"]])
  results[["tpms"]] <- xes

  return(results)
}






