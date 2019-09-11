library(tidyverse)
library(Biobase)
library(sva)

pharma <- read_rds('data/pharmacological_perturbation.rds')

ind <- unique(pharma$series_id[pharma$treatment_target == 'insulin response'])

se <- pharma[, pharma$series_id %in% ind]
se <- se[, order(se$treatment_type)]
pd <- as_tibble(pData(se))
mat <- exprs(se)

## remove na, and low intensities
ind <- apply(mat, 1, function(x) sum(is.na(x)))

mat <- mat[ind == 0,]

mat[mat < 0] <- 0
dim(mat)

ind <- colMeans(mat, na.rm = TRUE) > 10

mat[, ind] <- log2(mat[, ind] + 1)
dim(mat)
## normalize values
boxplot(mat)

mat1 <- limma::normalizeBetweenArrays(mat)
boxplot(mat1)

table(se$treatment_type, se$treatment_duration)

mod <- model.matrix(~time+treatment_duration+treatment_type, data = pd)
dim(mod)
n <- num.sv(mat1, mod = mod)

View(md)
gse <- c('GSE12929', 'GSE14004', 'GSE64075', 'GSE26207')
length(unique(md$series_id))
