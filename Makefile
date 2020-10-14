#!/bin/bash

# Define directory structure; for scripts
FIG_SRC=scripts/figures
TAB_SRC=scripts/tables

# Define directory structure; for output
MANUSCRIPT=manuscript
FIG_DIR=manuscript/figures
TAB_DIR=manuscript/tables

# Define directory structure; for intermediates
DATA=data
LOG=log
LOG_FIG=log/figures
LOG_TAB=log/tables

# Define commands
RFIG=R CMD BATCH --vanilla $< $(LOG_FIG)/$(<F).Rout
RTAB=R CMD BATCH --vanilla $< $(LOG_TAB)/$(<F).Rout

# All
all: ## Run all parts of the makefile
all: data figures tables clean

# Download data
data: ## Download the processed data
data: 
	wget -c -O data.zip https://ndownloader.figshare.com/articles/13088624/versions/1
	unzip -n data.zip -d $(DATA)
	
# Directories
dir_manuscript: ## Make manuscript directory tree
dir_manuscript:
	test ! -d $(MANUSCRIPT) && mkdir $(MANUSCRIPT) || exit 0
	test ! -d $(TAB_DIR) && mkdir $(TAB_DIR) || exit 0
	test ! -d $(FIG_DIR) && mkdir $(FIG_DIR) || exit 0

dir_logs: ## Make logs directory tree
dir_logs:
	test ! -d $(LOG) && mkdir $(LOG) || exit 0
	test ! -d $(LOG_FIG) && mkdir $(LOG_FIG) || exit 0
	test ! -d $(LOG_TAB) && mkdir $(LOG_TAB) || exit 0

figures: ## Generate the figures
figures: dir_manuscript \
	dir_logs \
	$(FIG_DIR)/timeline.png \
	$(FIG_DIR)/remove_batch.png \
	$(FIG_DIR)/Pparg_knockdown.png
	
tables: ## Generate the tables
tables: dir_manuscript \
	dir_logs \
	$(TAB_DIR)/datasets_pharmacological.tex \
	$(TAB_DIR)/datasets_genetic.tex

# Figures
$(FIG_DIR)/timeline.png: $(FIG_SRC)/timeline.R \
	$(DATA)/genetic_perturbation.rds \
	$(DATA)/pharmacological_perturbation.rds
	$(RFIG)
$(FIG_DIR)/remove_batch.png: $(FIG_SRC)/remove_batch.R \
	$(DATA)/genetic_perturbation.rds
	$(RFIG)
$(FIG_DIR)/Pparg_knockdown.png: $(FIG_SRC)/Pparg_knockdown.R \
	$(DATA)/genetic_perturbation.rds
	$(RFIG)
	
# Tables
$(TAB_DIR)/datasets_pharmacological.tex: $(TAB_SRC)/datasets_pharmacological.R \
	$(DATA)/pharmacological_perturbation.rds
	$(RTAB)
$(TAB_DIR)/datasets_genetic.tex: $(TAB_SRC)/datasets_genetic.R \
	$(DATA)/genetic_perturbation.rds
	$(RTAB)

# Clean Up
.PHONY: clean
clean: ## Clean up
clean:
	rm -f *.pdf
	rm -f *.RData

# Source: https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.PHONY: help
help: ## Print the current page
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'
.DEFAULT_GOAL := help
