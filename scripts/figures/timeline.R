# load library
library(tidyverse)
library(Biobase)
library(cowplot)

# load data
genetic <- read_rds('data/genetic_perturbation.rds')
pd <- as_tibble(pData(genetic))

p1 <- pd %>%
    filter(grepl('expression|down', treatment)) %>%
    group_by(treatment_target, treatment) %>%
    mutate(time = as.numeric(time),
           min = min(time),
           max = max(time),
           n = length(unique(time)),
           n2 = n()) %>%
    ungroup() %>%
    mutate(treatment_target = fct_reorder(treatment_target, n)) %>%
    ggplot() +
    geom_tile(aes(x = time, y = treatment_target, fill = series_id, width = 48), alpha = .2) +
    geom_vline(xintercept = c(0, 48, 144), color = 'gray', alpha = .5) +
    geom_segment(aes(x = min, xend = max, y = treatment_target, yend = treatment_target), size = 2.2, color = 'gray', alpha = .3) +
    geom_point(aes(x = time, y = treatment_target, shape = treatment), color = 'black') +
    geom_label(aes(x = -80, y = treatment_target, label = treatment_target), hjust = 1, label.size = NA) +
    scale_x_continuous(position = 'top',
                       limits = c(-150, 340),
                       breaks = seq(-2, 14, 2) * 24,
                       sec.axis = sec_axis(~.,
                                           breaks = c(-2, .5, 4, 10) * 24,
                                           labels = c('None', 'Early', 'Intermediate', 'Late'))) +
    scale_color_manual(values = c('darkred', 'darkgreen')) +
    scale_shape_manual(values = c(25, 24)) +
    theme_bw() +
    theme(legend.position = 'none',
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x.bottom = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank()) +
    labs(x = 'Time Point (hr)',
         y = '',
         color = 'Treatment',
         fill = 'Series ID')


pharma <- read_rds('data/pharmacological_perturbation.rds')

pd_pharma <- as_tibble(pData(pharma))

p2 <- pd_pharma %>%
    filter(grepl('drug', treatment)) %>%
    group_by(treatment_target, treatment) %>%
    mutate(time = as.numeric(time),
           min = min(time),
           max = max(time),
           n = length(unique(time)),
           n2 = n()) %>%
    ungroup() %>%
    mutate(treatment_target = fct_reorder(treatment_target, n)) %>%
    ggplot() +
    geom_tile(aes(x = time, y = treatment_target, fill = series_id, width = 48), alpha = .1) +
    geom_vline(xintercept = c(0, 48, 144), color = 'gray', alpha = .5) +
    geom_segment(aes(x = min, xend = max, y = treatment_target, yend = treatment_target), size = 2.2, color = 'gray', alpha = .3) +
    geom_point(aes(x = time, y = treatment_target), color = 'black') +
    geom_label(aes(x = -80, y = treatment_target, label = treatment_target), hjust = 1, label.size = NA) +
    scale_x_continuous(position = 'top',
                       limits = c(-300, 340),
                       breaks = seq(-2, 14, 2) * 24,
                       sec.axis = sec_axis(~.,
                                           breaks = c(-2, .5, 4, 10) * 24,
                                           labels = c('None', 'Early', 'Intermediate', 'Late'))) +
    theme_bw() +
    theme(legend.position = 'none',
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x.bottom = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank()) +
    labs(x = 'Time Point (hr)',
         y = '',
         color = 'Treatment',
         fill = 'Series ID')

plot_grid(p1, p2,
          scale = .9,
          labels = 'AUTO',
          label_fontface = 'plain',
          label_size = 10) %>%
    ggsave(plot = .,
           filename = 'manuscript/figures/timeline.png',
           width = 21, height = 15, units = 'cm')
