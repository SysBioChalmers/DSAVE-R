#' createPackageVariables
#'
#' Sets up all package data except the standard template, which was done manually.
#' Make sure the current directory is the package root before calling this function.
#'
#' @importFrom graphics hist
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#'

createPackageVariables <- function() {


#public dataset scores
#The numbers come from the publication, figure 3A.
datasetScoresHuman <- list(list("BC", 0.1378),
                           list("OC", 0.0389),
                           list("LC", 0.0283),
                           list("LIVC", 0.1102),
                           list("PBMC68k", 0.0352),
                           list("B10k", 0.0194),
                           list("CD4TMEM", 0.0139),
                           list("HCA CB", 0.0163),
                           list("CD8T", 0.0321))


usethis::use_data(datasetScoresHuman, overwrite = T)


#bulk variation
#This requires access to the data, which is not provided in the package
samp  <- as.matrix(read.table("data/tcellCD4ProfilesTMMNormalized.txt", header = T, sep = "\t", row.names = 1))
samp <- samp[,32:39]

#scale so the average count per sample is 10^6
normFact = 10^6 / mean(colSums(samp))
samp = samp * normFact



bulkMean1Vs1_05_100k <- DSAVEGetTotalVariationFromBulk(samp,
                                                       pool4samples = FALSE, upperBound = 100000,
                                                       lowerBound = 0.5, rescale = FALSE)
bulkMean4Vs4_05_100k <- DSAVEGetTotalVariationFromBulk(samp,
                                                      pool4samples = TRUE, upperBound = 100000, na.rm = TRUE,
                                                      lowerBound = 0.5, nComb = 10000L, rescale =FALSE)

bulkMean1Vs1_05_2 <- DSAVEGetTotalVariationFromBulk(samp,
                                                       pool4samples = FALSE, upperBound = 2,
                                                       lowerBound = 0.5, rescale = FALSE)
bulkMean4Vs4_05_2 <- DSAVEGetTotalVariationFromBulk(samp,
                                                      pool4samples = TRUE, upperBound = 2, na.rm = TRUE,
                                                      lowerBound = 0.5, nComb = 10000L, rescale =FALSE)

bulkMean1Vs1_100_100k <- DSAVEGetTotalVariationFromBulk(samp,
                                                       pool4samples = FALSE, upperBound = 100000,
                                                       lowerBound = 100, rescale = FALSE)
bulkMean4Vs4_100_100k <- DSAVEGetTotalVariationFromBulk(samp,
                                                      pool4samples = TRUE, upperBound = 100000, na.rm = TRUE,
                                                      lowerBound = 100, nComb = 10000L, rescale =FALSE)

bulkMean1Vs1_05_50 <- DSAVEGetTotalVariationFromBulk(samp,
                                                    pool4samples = FALSE, upperBound = 50,
                                                    lowerBound = 0.5, rescale = FALSE)
bulkMean4Vs4_05_50 <- DSAVEGetTotalVariationFromBulk(samp,
                                                   pool4samples = TRUE, upperBound = 50, na.rm = TRUE,
                                                   lowerBound = 0.5, nComb = 10000L, rescale =FALSE)

bulkMean1Vs1_2_100 <- DSAVEGetTotalVariationFromBulk(samp,
                                                    pool4samples = FALSE, upperBound = 100,
                                                    lowerBound = 2, rescale = FALSE)
bulkMean4Vs4_2_100 <- DSAVEGetTotalVariationFromBulk(samp,
                                                   pool4samples = TRUE, upperBound = 100, na.rm = TRUE,
                                                   lowerBound = 2, nComb = 10000L, rescale =FALSE)

bulkTotalVar1vs1 <- list(list("PseudoTPM: 0.5-100k", bulkMean1Vs1_05_100k),
                          list("PseudoTPM: 0.5-2", bulkMean1Vs1_05_2),
                          list("PseudoTPM: 100-100k", bulkMean1Vs1_100_100k),
                          list("PseudoTPM: 0.5-50", bulkMean1Vs1_05_50),
                          list("PseudoTPM: 2-100", bulkMean1Vs1_2_100))

bulkTotalVar4vs4 <- list(list("PseudoTPM: 0.5-100k", bulkMean4Vs4_05_100k),
                          list("PseudoTPM: 0.5-2", bulkMean4Vs4_05_2),
                          list("PseudoTPM: 100-100k", bulkMean4Vs4_100_100k),
                          list("PseudoTPM: 0.5-50", bulkMean4Vs4_05_50),
                          list("PseudoTPM: 2-100", bulkMean4Vs4_2_100))

usethis::use_data(bulkTotalVar1vs1, overwrite = T)
usethis::use_data(bulkTotalVar4vs4, overwrite = T)









#devtools::use_data(OBJECT_NAME, overwrite = T)
}
