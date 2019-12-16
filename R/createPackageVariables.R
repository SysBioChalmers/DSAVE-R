#' createPackageVariables
#'
#' Sets up all package data except the standard template, which was done manually.
#' Make sure the current directory is the package root before calling this function.
#'
#' @importFrom graphics hist
#' @importFrom utils read.table
#' @importFrom R.matlab readMat
#' @export
#' @author Johan Gustafsson, <gustajo@@chalmers.se>
#'

createPackageVariables <- function() {


#public dataset scores
#The numbers come from the publication, figure 3A.
datasetScoresHuman <- list(list("BC", 0.1371),
                           list("OC", 0.0385),
                           list("LC", 0.0270),
                           list("LIVC", 0.1096),
                           list("PBMC68k", 0.0351),
                           list("B10k", 0.0190),
                           list("CD4TMEM", 0.0137),
                           list("HCA CB", 0.0160),
                           list("CD8T", 0.0326))


usethis::use_data(datasetScoresHuman, overwrite = T)

datasetScoresHuman1000 <- list(list("BC", 0.1361),
                           list("OC", 0.0381),
                           list("LC", 0.0279),
                           list("LIVC", 0.1083),
                           list("PBMC68k", 0.0352),
                           list("B10k", 0.0194),
                           list("CD4TMEM", 0.0142),
                           list("HCA CB", 0.0174),
                           list("CD8T", 0.0327))

usethis::use_data(datasetScoresHuman1000, overwrite = T)

datasetScoresHuman500 <- list(list("BC", 0.1311),
                           list("OC", 0.0371),
                           list("LC", 0.0283),
                           list("LIVC", 0.1061),
                           list("PBMC68k", 0.0345),
                           list("B10k", 0.0193),
                           list("CD4TMEM", 0.0161),
                           list("HCA CB", 0.0180),
                           list("CD8T", 0.0319))

usethis::use_data(datasetScoresHuman500, overwrite = T)

#bulk variation
#This requires access to the data, which is not provided in the package
samp  <- as.matrix(read.table("data/scaledTMMMatrix.txt", header = T, sep = "\t", row.names = 1))


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

#bulkTotalVar1vs1[[1]][[2]] = bulkTotalVar1vs1[[1]][[2]]*5
usethis::use_data(bulkTotalVar1vs1, overwrite = T)
usethis::use_data(bulkTotalVar4vs4, overwrite = T)

#now save the pbmc68k cell types
#first, read the file and the data
pbmc68kmetadata  <- as.matrix(read.table("data/PBMC68000PatAFresh/filtered_matrices_mex/hg19/68k_pbmc_barcodes_annotation.tsv", header = T, sep = "\t"))
#test that the order in both files is the same. It is...
#library(Seurat)
#pbmc68kdata <- Read10X(data.dir = "data/PBMC68000PatAFresh/filtered_matrices_mex/hg19/")
#cn1 = colnames(pbmc68kdata)
#cn2 = pbmc68kmetadata[,3]
#sum(cn1 == cn2)
#sum(cn1 != cn2)

ctPbmc68k = as.factor(pbmc68kmetadata[,4])
usethis::use_data(ctPbmc68k, overwrite = T)


#import template stuff from matlab files
a = readMat("templInfo2000BinningInfo.mat")
tmp = a$c[,,1]
binInf2000 = list(binningInfo.means = tmp$means[1,], binningInfo.lbs = tmp$lbs[1,], binningInfo.ubs = tmp$ubs[1,])
usethis::use_data(binInf2000, overwrite = T)

a = readMat("templInfo1000BinningInfo.mat")
tmp = a$a[,,1]
binInf1000 = list(binningInfo.means = tmp$means[1,], binningInfo.lbs = tmp$lbs[1,], binningInfo.ubs = tmp$ubs[1,])
usethis::use_data(binInf1000, overwrite = T)

a = readMat("templInfo500BinningInfo.mat")
tmp = a$b[,,1]
binInf500 = list(binningInfo.means = tmp$means[1,], binningInfo.lbs = tmp$lbs[1,], binningInfo.ubs = tmp$ubs[1,])
usethis::use_data(binInf500, overwrite = T)

a = readMat("templInfo2000UMIDistr.mat")
UMIDistr2000 = a$c[1,]
usethis::use_data(UMIDistr2000, overwrite = T)

a = readMat("templInfo1000UMIDistr.mat")
UMIDistr1000 = a$a[1,]
usethis::use_data(UMIDistr1000, overwrite = T)

a = readMat("templInfo500UMIDistr.mat")
UMIDistr500 = a$b[1,]
usethis::use_data(UMIDistr500, overwrite = T)





#devtools::use_data(OBJECT_NAME, overwrite = T)
}
