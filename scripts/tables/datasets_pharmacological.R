# load library
library(tidyverse)
library(SummarizedExperiment)
library(bib2df)
library(xtable)

# load data
pharma <- read_rds('data/pharmacological_perturbation.rds')

pd <- as_tibble(colData(pharma)) 

bib <- bib2df('data/studies.bib') %>%
    select(pmid = PMID,
           texkey = BIBTEXKEY)

pd %>%
    filter_at(vars(time, treatment, treatment_type, treatment_target, pmid),
              function(x) !grepl('none', x)) %>%
    group_by(series_id) %>%
    summarise_at(vars(time, treatment_target, treatment_type, pmid),
                 ~paste(unique(.x), collapse = '/ ')) %>%
    left_join(bib) %>%
    dplyr::select(-pmid) %>%
    mutate(texkey = paste0('\\cite{', texkey, '}'),
           treatment_target = ifelse(series_id == 'GSE56688', 'insulin response', treatment_target)) %>%
    arrange(treatment_target) %>%
    setNames(c('GEO ID', 'Time (hr)', 'Target', 'Treatment', 'Ref.')) %>%
    xtable(align = 'clcp{.2\\textwidth}p{.35\\textwidth}c') %>%
    print(floating = FALSE,
          include.rownames = FALSE,
          booktabs = TRUE,
          sanitize.text.function = identity,
          comment = FALSE,
          file = 'manuscript/tables/datasets_pharmacological.tex')
