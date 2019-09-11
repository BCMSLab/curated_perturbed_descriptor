# load library
library(tidyverse)
library(Biobase)
library(sva)
library(plot3D)

# load data
## select none treated, MDI induced
genetic <- read_rds('data/genetic_perturbation.rds')
gse <- genetic[, genetic$treatment == 'none' & genetic$media == 'MDI']
gmat <- exprs(gse)

## remove na, and low intensities
ind <- apply(gmat, 1, function(x) sum(is.na(x)))
gmat <- gmat[ind == 0,]

gmat[gmat < 0] <- 0
gmat <- gmat[rowMeans(gmat) > 10,]

## normalize values
gmat1 <- limma::normalizeBetweenArrays(log2(gmat+1))

## remove gpl, series batch
gmat2 <- ComBat(gmat1,
               batch = paste0(gse$gpl, gse$series_id))

## visulaize

gpcs <- map(list(gmat1, gmat2), function(x) prcomp(dist(t(x))))


## select none treated, MDI induced
pharma <- read_rds('data/pharmacological_perturbation.rds')
pse <- pharma[, pharma$treatment == 'none' & pharma$media == 'MDI']
pmat <- exprs(pse)

## remove na, and low intensities
ind <- apply(pmat, 1, function(x) sum(is.na(x)))
pmat <- pmat[ind == 0,]

pmat[pmat < 0] <- 0
pmat <- pmat[rowMeans(pmat) > 10,]

## normalize values
pmat1 <- limma::normalizeBetweenArrays(log2(pmat+1))

## remove gpl, series batch
pmat2 <- ComBat(pmat1,
                batch = paste0(pse$gpl, pse$series_id))

## visulaize

ppcs <- map(list(pmat1, pmat2), function(x) prcomp(dist(t(x))))

png(filename = 'manuscript/figures/none_pca.png',
    height = 16, width = 21, units = 'cm', res = 300)
par(mfrow = c(2,2),
    mar = c(0, 1, 0, 5))

gpl <- as.factor(gse$gpl)

scatter3D(
    x = gpcs[[1]]$x[, 1],
    y = gpcs[[1]]$x[, 2],
    z = gpcs[[1]]$x[, 3],
    pch = 19,
    xlab = 'PC1', ylab = 'PC2', zlab = 'PC3',
    colvar = as.numeric(gpl),
    colkey = list(at = unique(as.numeric(gpl)),
                  labels = levels(gpl),
                  length = .5,
                  cex.axis = .8))
mtext('A', 3, 0, padj = 2, adj = 0)

time <- as.numeric(gse$time)
scatter3D(
    x = gpcs[[2]]$x[, 1],
    y = gpcs[[2]]$x[, 2],
    z = gpcs[[2]]$x[, 3],
    pch = 19,
    xlab = 'PC1', ylab = 'PC2', zlab = 'PC3',
    colvar = as.numeric(gse$time),
    colkey = list(at = c(-20, 40, 140, 500),
                  labels = c('None', 'Early', 'Intermediate', 'Late'),
                  length = .5,
                  cex.axis = .8))
mtext('B', 3, 0, padj = 2, adj = 0)

gpl <- as.factor(pse$gpl)

scatter3D(
    x = ppcs[[1]]$x[, 1],
    y = ppcs[[1]]$x[, 2],
    z = ppcs[[1]]$x[, 3],
    pch = 19,
    xlab = 'PC1', ylab = 'PC2', zlab = 'PC3',
    colvar = as.numeric(gpl),
    colkey = list(at = unique(as.numeric(gpl)),
                  labels = levels(gpl),
                  length = .5,
                  cex.axis = .8))
mtext('C', 3, 0, padj = 2, adj = 0)

time <- as.numeric(pse$time)
scatter3D(
    x = ppcs[[2]]$x[, 1],
    y = ppcs[[2]]$x[, 2],
    z = ppcs[[2]]$x[, 3],
    pch = 19,
    xlab = 'PC1', ylab = 'PC2', zlab = 'PC3',
    colvar = as.numeric(pse$time),
    colkey = list(at = c(-20, 40, 140, 500),
                  labels = c('None', 'Early', 'Intermediate', 'Late'),
                  length = .5,
                  cex.axis = .8))
mtext('D', 3, 0, padj = 2, adj = 0)
dev.off()
