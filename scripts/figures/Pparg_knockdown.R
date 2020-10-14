# load library
library(tidyverse)
library(reshape2)
library(cowplot)
library(SummarizedExperiment)
library(sva)
library(limma)
library(org.Mm.eg.db)
library(clusterProfiler)

# load data
genetic <- read_rds('data/genetic_perturbation.rds')

## select pparg knockdown
genetic <- genetic[, genetic$series_id %in% c('GSE14004', 'GSE12929')]
genetic$treatment <- factor(genetic$treatment, levels = c('none', 'knockdown'))

## remove na
genetic <- genetic[apply(assay(genetic), 1, function(x) sum(is.na(x))) == 0,]

## add stage variable
genetic$stage <- cut(genetic$time,
                     breaks = c(-100, 0, 48, 140, 1000),
                     labels = c('None', 'Early', 'Intermediate', 'Late'))

## normalize and remove batch effects
assay(genetic, 'normalized') <- normalizeBetweenArrays((assay(genetic)))
assay(genetic, 'nobatch') <- ComBat(assay(genetic, 'normalized'), genetic$series_id)

## make a model matrix
mod <- model.matrix(~treatment*time, data = colData(genetic))

## apply regression models
fit <- lmFit(log2(assay(genetic, 'normalized') + 1), mod)
fit <- eBayes(fit)

res <- map_df(c(knockdown = 2,
                time = 3,
                knockdown_time = 4), 
              function(x) {
                  topTable(fit, 
                           coef = x,
                           number = Inf,
                           genelist = rownames(genetic))
              },
              .id = 'contrast') %>%
    as_tibble()

targets <- readxl::read_excel('data/PPAR gamma target genes.xls') %>%
    filter(`Tissue / Cell` == 'adipose tissue') %>%
    pull(`Gene Symbol`)

length(intersect(rownames(genetic), targets))

p1 <- res %>%
    filter(contrast == 'knockdown', ID %in% c('Pparg', targets)) %>%
    mutate(lab = ifelse(abs(logFC) > 1 & P.Value < .05, ID, '')) %>%
    ggplot(aes(x = logFC, y = -log10(P.Value))) +
    geom_hline(yintercept = -log10(.05), lty = 2, color = 'red') +
    geom_vline(xintercept = c(-1, 1), lty = 2, color = 'red') +
    geom_point(size = 3) +
    geom_text(aes(y = -log10(P.Value) + .7, label = lab)) +
    #lims(x = c(-4.8, 4.9)) +
    labs(x = 'Fold-change (log_2)',
         y = 'P-value (-log_10)')

term2gene <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = keys(org.Mm.eg.db, 'SYMBOL'),
                                   columns = 'GO',
                                   'SYMBOL') %>%
    dplyr::select(term = GO, gene = SYMBOL)
term2name <- Term(GO.db::GOTERM) %>%
    melt(value.name = 'name') %>%
    rownames_to_column('term') %>%
    filter(term %in% unique(term2gene$term))

go_res <- map_df(unique(res$contrast), 
                 function(x) {
                     deg <- filter(res, contrast == x) %>% pull(logFC)
                     names(deg) <- filter(res, contrast == x) %>% pull(ID)
                     deg <- sort(deg, decreasing = TRUE)
                     
                     set.seed(123)
                     
                     GSEA(deg,
                          TERM2GENE = term2gene,
                          TERM2NAME = term2name,
                          seed = TRUE,
                          pAdjustMethod = 'fdr',
                          pvalueCutoff = 1)@result %>%
                         mutate(contrast = x)
                 })

p2 <- go_res %>%
    filter(ID %in% c('GO:0050873', 'GO:0045600', 'GO:0033993', 'GO:0030374', 'GO:0019395')) %>%
    filter(contrast == 'knockdown') %>%
    ggplot(aes(x = ID, y = enrichmentScore, fill = -log10(pvalue))) +
    geom_col() +
    coord_flip() +
    theme(legend.position = 'top') +
    labs(x = '', y = 'Enrichment Score', fill = 'P-value (-log_10)')

plot_grid(p1,p2,
          ncol = 2,
          labels = 'AUTO',
          scale = .9,
          label_fontface = 'plain') %>%
ggsave(plot = .,
       filename = 'manuscript/figures/Pparg_knockdown.png',
       height = 12, width = 24, units = 'cm')
