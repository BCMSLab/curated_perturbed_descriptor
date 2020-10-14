# load library
library(tidyverse)
library(reshape2)
library(cowplot)
library(SummarizedExperiment)
library(sva)
library(limma)

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

## make a model matrix
mod <- model.matrix(~treatment, data = colData(genetic))

## normalize and remove batch effects
assay(genetic, 'normalized') <- normalizeBetweenArrays((assay(genetic)))
assay(genetic, 'nobatch') <- ComBat(assay(genetic, 'normalized'),
                                    genetic$series_id,
                                    mod)

gpcs <- list(before = 'exprs',
             after = 'nobatch') %>%
    map(function(x) prcomp(dist(t(assay(genetic, x)))))
    
# variance explained

p1 <- map_df(gpcs, function(x) {
    tibble(
        comp = factor(colnames(x$x), colnames(x$x)),
        variance = (x$sdev ^ 2) / sum(x$sdev ^ 2)
    )
}, .id = 'when') %>%
    filter(comp %in% paste0('PC', 1:3)) %>%
    mutate(comp = as.factor(comp)) %>%
    ggplot(aes(x = comp, y = variance, color = when, group = when)) +
    geom_line(size = 1.2) +
    theme(legend.position = 'none') +
    scale_y_continuous(
        limits = c(0, 1),
        labels = scales::number_format(accuracy = 0.1,
                                       decimal.mark = '.')) +
    annotate('text', x = 1.5, y = .8, label = 'before') +
    annotate('text', x = 1, y = .4, label = 'after') +
    labs(x = 'Principal Components',
         y = 'Variance')

p2 <- list(stage = genetic$stage,
     study = genetic$series_id) %>%
    map_df(function(v) {
        gpcs %>%
            map_df(function(g) {
                as_tibble(g$x) %>%
                    map_df(function(x) {
                        if (!is.numeric(v)) {
                            kruskal.test(x, v)$p.value
                        } else {
                            cor.test(x, v)$p.value
                        }
                    })
            }, .id = 'when')
    }, .id = 'variable') %>%
    gather(comp, pval, starts_with('PC')) %>%
    filter(comp %in% paste0('PC', 1:3)) %>%
    ggplot(aes(x = comp, y = -log10(pval))) +
    geom_col() +
    facet_grid(variable ~ when, scales = 'free') +
    labs(x = 'Principal Components',
         y = 'Correlation -log_10(p-value)')

pca_df <- map_df(gpcs, function(x) {
    as_tibble(melt(x$x, varnames = c('sample_id', 'comp')))
}, .id = 'when') %>%
    spread(comp, value) %>%
    left_join(as_tibble(colData(genetic)))

p3 <- pca_df %>%
    filter(when == 'before') %>%
    ggplot(aes(x = scale(PC1), y = scale(PC2), color = series_id)) +
    geom_point(size = 3) +
    theme(legend.position = 'top') +
    theme(legend.position = 'none') +
    labs(x = 'Principal Component 1',
         y = 'Principal Component 2')

p4 <- pca_df %>%
    filter(when == 'after') %>%
    ggplot(aes(x = scale(PC3), y = scale(PC2), color = stage)) +
    geom_point(size = 3) +
    theme(legend.position = 'top') +
    theme(legend.position = 'none') +
    labs(x = 'Principal Component 3',
         y = 'Principal Component 2')

plot_grid(p1,p2,p3,p4,
          ncol = 2,
          labels = 'AUTO',
          scale = .9,
          label_fontface = 'plain') %>%
    ggsave(plot = .,
           filename = 'manuscript/figures/remove_batch.png',
           height = 16, width = 16, units = 'cm')

