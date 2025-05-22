#### Jinxin Meng, 20240125, 20250318 ####
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, openxlsx, rstatix)
source('/code/R_func/calcu_diff.R')

#### SchirmerM_2018.PRJNA389280 ####
group <- read.delim('SchirmerM_2018.PRJNA389280.sample_group') %>% 
  dplyr::select(sample = run, group = group2)

data <- read.delim('SchirmerM_2018.PRJNA389280.tsv.gz') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('CT', 'CD', 'UC')))

comparisons <- dplyr::select(data, sample, value = rela_ab) %>% 
  calcu_diff(group) %>% 
  filter(pval < 0.05) %>% 
  pull(comparison) %>% 
  strsplit(split = '_vs_')

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '',
         size = 1, width = .7, palette = c('#69a7bc','#E69F00','#f08178'),
         ylab = 'Relative abundance(%)', outlier.shape = NA,
         title = 'SchirmerM_2018') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  stat_compare_means(comparisons = comparisons, label = 'p.signif', 
                     label.y = .135, tip.length = .03, step.increase = .07, 
                     vjust = .6, size = 5) +
  theme(aspect.ratio = 1.2) +
  theme(plot.title = element_text(hjust = .5, face = 'bold'))

#### LloydPriceJ_2019.PRJNA398089 ####
group <- read.delim('LloydPriceJ_2019.PRJNA398089.sample_group') %>% 
  dplyr::select(sample = run, group = group2)

data <- read.delim('LloydPriceJ_2019.PRJNA398089.tsv.gz') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('CT', 'CD', 'UC')))

comparisons <- dplyr::select(data, sample, value = rela_ab) %>% 
  calcu_diff(group) %>% 
  filter(pval < 0.05) %>% 
  pull(comparison) %>% 
  strsplit(split = '_vs_')

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '',
         size = 1, width = .7, palette = c('#69a7bc','#E69F00','#f08178'),
         ylab = 'Relative abundance(%)', outlier.shape = NA,
         title = 'LloydPriceJ_2019') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .15)) +
  stat_compare_means(comparisons = comparisons, label = 'p.signif', 
                     label.y = .14, tip.length = .012, step.increase = .07, 
                     vjust = .6, size = 5) +
  theme(aspect.ratio = 1.2) +
  theme(plot.title = element_text(hjust = .5, face = 'bold'))

#### HeQ_2017.PRJEB15371 ####
group <- read.delim('HeQ_2017.PRJEB15371.sample_group')

data <- read.delim('HeQ_2017.PRJEB15371.tsv.gz') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('Control', 'CD')))

comparisons <- dplyr::select(data, sample, value = rela_ab) %>% 
  calcu_diff(group) %>% 
  filter(pval < 0.05) %>% 
  pull(comparison) %>% 
  strsplit(split = '_vs_')

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '',
         size = 1, width = .7, palette = c('#69a7bc','#E69F00'),
         ylab = 'Relative abundance(%)', outlier.shape = NA,  
         title = 'HeQ_2017') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .2)) +
  stat_compare_means(comparisons = comparisons, label = 'p.signif', 
                     label.y = .17, tip.length = .01, step.increase = .07, 
                     vjust = .6, size = 5) +
  theme(aspect.ratio = 1.6) +
  theme(plot.title = element_text(hjust = .5, face = 'bold'))

#### YanQ_2023c.PRJEB67456 ####
group <- read.delim('YanQ_2023c.PRJEB67456.sample_group') %>% 
  dplyr::select(sample = run, group)

data <- read.delim('YanQ_2023c.PRJEB67456.tsv.gz') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('healthy', 'CD', 'UC')))

comparisons <- dplyr::select(data, sample, value = rela_ab) %>% 
  calcu_diff(group) %>% 
  filter(pval < 0.05) %>% 
  pull(comparison) %>% 
  strsplit(split = '_vs_')

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '',
         size = 1, width = .7, palette = c('#69a7bc','#E69F00','#f08178'),
         ylab = 'Relative abundance(%)', outlier.shape = NA,  
         title = 'YanQ_2023c') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .5)) +
  stat_compare_means(comparisons = comparisons, label = 'p.signif', 
                     label.y = .46, tip.length = .01, step.increase = .07, 
                     vjust = .6, size = 5) +
  theme(aspect.ratio = 1.2) +
  theme(plot.title = element_text(hjust = .5, face = 'bold'))

