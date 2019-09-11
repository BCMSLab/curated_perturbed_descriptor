# load library
library(tidyverse)
library(Biobase)
library(bib2df)
library(xtable)

# load data
genetic <- read_rds('data/genetic_perturbation.rds')

pd <- as_tibble(pData(genetic)) 

bib <- bib2df('data/studies.bib') %>%
    select(pmid = PMID,
           texkey = BIBTEXKEY)
pd %>%
    filter(perturbation_type == 'genetic') %>%
    filter_at(vars(time, treatment, treatment_type, treatment_target, pmid),
              function(x) !grepl('none', x)) %>%
    group_by(series_id) %>%
    summarise_at(vars(time, treatment, treatment_target, pmid),
                 ~paste(unique(.x), collapse = '/ ')) %>%
    left_join(bib) %>%
    mutate(texkey = paste0('\\cite{', texkey, '}')) %>%
    setNames(c('GEO ID', 'Time (hr)', 'Treatment', 'Target', 'PMID', 'Ref.')) %>%
    xtable(align = 'clp{.15\\textwidth}lp{.2\\textwidth}lc') %>%
    print(floating = FALSE,
          include.rownames = FALSE,
          booktabs = TRUE,
          sanitize.text.function = identity,
          comment = FALSE,
          file = 'manuscript/tables/datasets_genetic.tex')

# ids <- read_tsv('../curatedAdipoArray/inst/scripts/curated_metadata.tsv') %>%
#     pull(pmid) %>%
#     unique()
# 
# pm <- entrez_summary(db = 'pubmed', ids)
# 
# map_df(pm, function(x) {
#     tibble(
#         category = 'ARTICLE',
#         title = x$title,
#         year = str_sub(x$pubdate, 0, 4),
#         journal = x$fulljournalname,
#         author = paste(x$authors$name, collapse = ','),
#         pages = x$pages,
#         volume = x$volume,
#         pmid = x$uid,
#         CATEGORY = 'article'
#     )
# }) %>%
#     df2bib('data/studies.bib')
# 
