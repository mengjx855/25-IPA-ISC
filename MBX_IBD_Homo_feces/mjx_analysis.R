#### Jinxin Meng, 20250217, 20250323 ####
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, openxlsx)
library(ropls)
library(clusterProfiler)
source('/code/R_func/difference_analysis.R')
source('/code/R_func/plot_PCA.R')
source('/code/R_func/calcu_diff.R')

#### info ####
group_color <- c('CD' = '#E69F00', 'HC' = '#69a7bc')

info <- read.delim('info.txt')
group <- read.delim('group.txt')
profile <- read.delim('profile.txt', row.names = 1, check.names = F)

#### class aggregate ####
profile_x <- profile %>% 
  mutate(class = info$class[match(rownames(profile), info$cpd_id)],
         class = replace(class, class == '', 'Other')) %>% 
  aggregate(. ~ class, ., sum) %>% 
  column_to_rownames('class')

diff <- difference_analysis(profile_x, group, comparison = names(group_color))

#### diff ####
plsda <- opls(data.frame(t(profile)), y = pull(group, group), orthoI = 0)
diff <- difference_analysis(profile, group, comparison = names(group_color))

diff <- data.frame(vip = plsda@vipVn) %>% 
  rownames_to_column('name') %>% 
  left_join(diff, ., by = 'name') %>% 
  mutate(enriched = ifelse((vip > 1 | pval < 0.05) & 
                             log2FC > 0, 'CD', 
                           ifelse((vip > 1 | pval < 0.05) & 
                                    log2FC < 0, 'HC', 'none'))) %>% 
  left_join(info, c('name' = 'cpd_id'))

plot_data <- diff %>% 
  dplyr::select(name, cpd_name, pval, log2FC, enriched) %>% 
  mutate(.enriched = ifelse(is.na(enriched), 'none', enriched),
         .pval = -log10(pval),
         .log2FC = ifelse(log2FC < -3.9, 3.9, log2FC),
         .label = ifelse(.pval > 4, cpd_name, ''))
  
ggscatter(plot_data, 'log2FC', '.pval', color = '.enriched',size = 1, 
          xlab = 'log2FoldChange', ylab = '-log10 P-value',
          palette = c(group_color, 'none' = 'grey')) +
  geom_text(aes(label = .label), size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_vline(xintercept = c(0), linetype = 'dashed') +
  theme(aspect.ratio = 1)
ggsave('diff.volcano.pdf', width = 5, height = 5)

#### enrich ####
bg_data <- read.delim('cpd2path_enrichment.tsv')

eKEGG <- diff %>% 
  filter(enriched != 'none') %>% 
  filter(KEGG != '') %>% 
  pull(KEGG) %>% 
  enricher(TERM2GENE = bg_data, minGSSize = 1, pvalueCutoff = 1, qvalueCutoff = 1) %>% 
  data.frame


paths <- c('map00270','map00230','map00350','map00680','map00380',
           'map00340','map00460','map00920','map00650','map00010')

plot_data <- eKEGG %>% 
  mutate(geneRatio = map_vec(strsplit(GeneRatio, '/'), \(x) as.numeric(x[1])/as.numeric(x[2])),
         KEGG = map_vec(strsplit(ID, ':'), \(x) x[1])) %>% 
  filter(KEGG %in% paths) %>% 
  arrange(geneRatio)

ggscatter(plot_data, 'ID', 'geneRatio', fill = 'pvalue', rotate = T, size = 5,
          shape = 21, legend = 'right', xlab = '', ylab = 'Gene Ratio') +
  scale_fill_viridis_c(begin = .6) +
  theme(aspect.ratio = 2)
