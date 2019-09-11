# load library
library(tidyverse)
library(Biobase)
library(bib2df)
library(xtable)

# load data
pharma <- read_rds('data/pharmacological_perturbation.rds')
pd <- as_tibble(pData(pharma)) 

bib <- bib2df('data/studies.bib') %>%
    select(pmid = PMID,
           texkey = BIBTEXKEY)

pd %>%
    filter_at(vars(time, treatment, treatment_type, treatment_target, pmid),
              function(x) !grepl('none', x)) %>%
    group_by(series_id) %>%
    summarise_at(vars(time, treatment, treatment_target, pmid),
                 ~paste(unique(.x), collapse = '/ ')) %>%
    left_join(bib) %>%
    mutate(texkey = paste0('\\cite{', texkey, '}')) %>%
    setNames(c('GEO ID', 'Time (hr)', 'Treatment', 'Target', 'PMID', 'Ref.')) %>%

    xtable(align = 'clcp{.15\\textwidth}p{.2\\textwidth}lc') %>%
    print(floating = FALSE,
          include.rownames = FALSE,
          booktabs = TRUE,
          sanitize.text.function = identity,
          comment = FALSE,
          file = 'manuscript/tables/datasets_pharmacological.tex')
