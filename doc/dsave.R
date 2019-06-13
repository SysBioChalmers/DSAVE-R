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


