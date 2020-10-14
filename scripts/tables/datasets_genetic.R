# load library
library(tidyverse)
library(SummarizedExperiment)
library(bib2df)
library(xtable)

# load data
genetic <- read_rds('data/genetic_perturbation.rds')

pd <- as_tibble(colData(genetic)) 

bib <- bib2df('data/studies.bib') %>%
    dplyr::select(pmid = PMID,
                  texkey = BIBTEXKEY)
pd %>%
    filter(perturbation_type == 'genetic') %>%
    filter_at(vars(time, treatment, treatment_type, treatment_target),
              function(x) !grepl('none', x)) %>%
    group_by(series_id) %>%
    summarise_at(vars(time, treatment, treatment_target, pmid),
                 ~paste(unique(.x), collapse = '/ ')) %>%
    left_join(bib) %>%
    dplyr::select(-pmid) %>%
    mutate(texkey = paste0('\\cite{', texkey, '}')) %>%
    arrange(treatment_target) %>%
    setNames(c('GEO ID', 'Time (hr)', 'Treatment', 'Target', 'Ref.')) %>%
    xtable(align = 'clp{.15\\textwidth}lp{.2\\textwidth}c') %>%
    print(floating = FALSE,
          include.rownames = FALSE,
          booktabs = TRUE,
          sanitize.text.function = identity,
          comment = FALSE,
          file = 'manuscript/tables/datasets_genetic.tex')
