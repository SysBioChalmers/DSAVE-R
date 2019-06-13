## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, fig.width = 7,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(ggplot2)
library(DSAVE)
library(progress)

## ------------------------------------------------------------------------
bcells <- loadOrDownloadB10k()
tcells <- loadOrDownloadT4k()

## ------------------------------------------------------------------------
scTotalVariation_bcells <- list()
scTotalVariation_tcells <- list()
pools <- c(100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)
pools2 <- c(100, 200, 500, 750, 1000, 1500, 2250)
pb <- progress_bar$new(format = "Calculating [:bar] :percent eta: :eta",
                       total = length(pools) + length(pools2), clear = FALSE)

for(poolSize in pools){
  scTotalVariation_bcells[[as.character(poolSize)]] <- DSAVEGetTotalVariationPoolSize(bcells,
        poolSize = poolSize, upperBound = 1e5, lowerBound = 5e-1)
  pb$tick();
}

for(poolSize in pools2){
  scTotalVariation_tcells[[as.character(poolSize)]] <- DSAVEGetTotalVariationPoolSize(tcells,
          poolSize = poolSize, upperBound = 1e5, lowerBound = 5e-1)
  pb$tick();
}

pb$terminate()


## ---- fig.show = 'asis'--------------------------------------------------
bulkMean1Vs1 <- rep(mean(bulkTotalVar1vs1[[1]][[2]]), 2)
bulkMean4Vs4 <- rep(mean(bulkTotalVar4vs4[[1]][[2]]), 2)
sc <- rbind(
  cbind(pools, sapply(scTotalVariation_bcells, 
                      function(v) mean(v)), "B10k", "sc"), 
  cbind(pools2, sapply(scTotalVariation_tcells, 
                       function(v) mean(v)), "T4k", "sc"), 
  cbind(c(pools[1],pools[length(pools)]), 
        bulkMean1Vs1, "single bulk sample", "bulk"),
  cbind(c(pools[1],pools[length(pools)]), 
        bulkMean4Vs4, "mean of 4 bulk samples", "bulk"))

sc <- as.data.frame(sc)
colnames(sc) <- c("PoolSize", "TotalVariation", "Dataset", "DataType")
sc$PoolSize <- as.numeric(as.character(sc$PoolSize))
sc$TotalVariation <- as.numeric(as.character(sc$TotalVariation))

ggplot(sc, aes(x=PoolSize, y=TotalVariation, group = Dataset)) +
  ggtitle("Variation per Cell Pool Size, CPM > 0.5") +
  geom_line(aes(color = Dataset, linetype = DataType)) +
  ylim(0, max(sc$TotalVariation)) +
  theme(legend.title = element_text(size = 10),
        legend.justification = c(1, 1),
        legend.position = c(1, 1),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10)) + 
  #ylab(bquote('Variation ( ~ R[mean])')) +
  labs(y=expression(paste(Variation (R[mean] ))),
       x="Pool size (number of cells)") +
  guides(linetype = FALSE)


## ------------------------------------------------------------------------
templInfo <- DSAVEGetStandardTemplate()

## ------------------------------------------------------------------------
tcells_btm <- DSAVECalcBTMScore(tcells, templInfo, skipAlignment=FALSE,
                                iterations = 15, useLogTransform=FALSE, 
                                logTPMAddon=1, silent=FALSE)

## ------------------------------------------------------------------------
data(datasetScoresHuman)

## ---- fig.show = 'asis'--------------------------------------------------
scores <- data.frame(t(sapply(datasetScoresHuman, 
                              function(x) c(x[[1]], x[[2]]))),
                     stringsAsFactors = F)

scores <- rbind.data.frame(scores,c("T4k", tcells_btm$DSAVEScore))
scores$X2 <- as.numeric(scores$X2)

ggplot(scores, aes( X1)) +
  geom_bar(aes(weight = X2)) +
  coord_flip() +
  ggtitle("Cell-to-cell Variation for Different Datasets") +
  ylab("DSAVE BTM score") + xlab("") +
  theme(axis.text.y = element_text(size = 10))


## ------------------------------------------------------------------------
bcells_gene_variation <- DSAVEGetGeneVariation(as.matrix(bcells), lb=10, iterations = 100, maxNumCells=2000, silent=FALSE)

## ---- fig.show = 'asis'--------------------------------------------------
id <- order(bcells_gene_variation$logCVDifference, decreasing = T)
df <- as.data.frame(cbind(1:length(id), 
                          bcells_gene_variation$logCVDifference[id]))

ggplot(df, aes(x=V1, y=V2, colour = "B10k")) + 
  geom_line() +
  ggtitle("Variation per Gene") +
  xlab("Gene index") + ylab("BTM variation") +
  theme(legend.justification = c(1, 1),
        legend.title=element_blank(),
        legend.position = c(1, 1),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10))

## ---- fig.show = 'asis'--------------------------------------------------
pvals <- bcells_gene_variation$pVals
id <- order(pvals, decreasing = T)

df <- data.frame(gene = bcells_gene_variation$genes[id], 
                 index = 1:length(pvals[id]), pvals = pvals[id])
df$gene <- as.character(df$gene)

ggplot(df, aes(x=index, y=pvals, colour = "B10k")) + 
  geom_line() +
  ggtitle("p values per gene") +
  xlab("Gene index") + ylab("p value") +
  theme(legend.justification = c(1, 1),
        legend.title=element_blank(),
        legend.position = c(1, 1),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10))


## ------------------------------------------------------------------------
bcells_id <- sample(1:dim(bcells)[2], 1000)
bcells_cell_divergence <- DSAVEGetSingleCellDivergence(bcells[,bcells_id],
                              minUMIsPerCell = 200, tpmLowerBound = 0,
                              iterations = 15,silent=FALSE)

## ------------------------------------------------------------------------

numUMIshcatSub2 <- colSums(as.matrix(bcells[,bcells_id]))
linevalYs <- seq(from =  500, to = 10500, by = 250)
linevalXes <- rep(0, length(linevalYs))
for(i in 1:length(linevalYs)){
  ii <- linevalYs[i]
  lb <- ii - 250
  ub <- ii + 250
  sel <- numUMIshcatSub2 >= lb & numUMIshcatSub2 <= ub
  if(sum(sel)!=0){
    linevalXes[i] <- mean(bcells_cell_divergence[sel])
  }else{
    linevalXes[i] <- NA
  }
}
linevalYs <- c(linevalYs, 
               rep(NA, length(numUMIshcatSub2) - length(linevalYs)))
linevalXes <- c(linevalXes, 
                rep(NA, length(numUMIshcatSub2) - length(linevalXes)))
linevalYs[is.na(linevalXes)] <- NA
id <- order(linevalYs, decreasing = F)
linevalXes <- linevalXes[id]
linevalYs <- linevalYs[id]

## ---- fig.show = 'asis'--------------------------------------------------
df <- as.data.frame(cbind(bcells_cell_divergence, colSums(as.matrix(bcells[,bcells_id]))))

ggplot(df, aes(x=bcells_cell_divergence, y=V2, colour = "Individual cell")) + 
  geom_point() +
  ggtitle("UMI Counts vs Cell Divergence") +
  xlab("Log-likelihood") + ylab("UMI counts") +
  theme(legend.justification = c(1, 1),
        legend.title=element_blank(),
        legend.position = c(1, 1),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10)) +
  geom_path(aes(linevalXes, linevalYs, colour="Mean log-likelihood"))

## ---- fig.show = 'asis'--------------------------------------------------

id <- order(bcells_cell_divergence, decreasing = F)
df <- as.data.frame(cbind(1:length(bcells_cell_divergence), bcells_cell_divergence[id]))
ggplot(df, aes(x=V1, y=V2, colour = "B10k")) + 
  geom_line() +
  ggtitle("Cell Divergence") +
  xlab("Cell index") + ylab("Log-likelihood") +
  theme(legend.justification = c(1, 0),
        legend.title=element_blank(),
        legend.position = c(1, 0),
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10))

