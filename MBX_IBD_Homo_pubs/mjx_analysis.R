#### Jinxin Meng, 20250211, 20250216 ####
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, ggpmisc, openxlsx)
pacman::p_load(ropls, clusterProfiler)
set.seed(123)
source('/code/R_func/difference_analysis.R')
source('/code/R_func/profile_process.R')

#### LloydPriceJ_2019 ####
info <- read.delim('LloydPriceJ_2019.mbt_info.txt')

info_x <- info %>% 
  filter(Match != '') %>% 
  dplyr::select(-cpd_id, -cpd_name) %>% 
  distinct() %>% 
  add_column(cpd_id = sprintf('cpd_%03d', seq(1, nrow(.))), .before = 1)

profile <- read.delim('LloydPriceJ_2019.profile.txt', check.names = F) %>% 
  map_df(\(x) replace_na(x, 0)) %>% 
  column_to_rownames('name') %>% 
  mutate(name = info$Match[match(rownames(.), info$cpd_id)]) %>% 
  filter(name %in% info_x$Match) %>% 
  group_by(name) %>% 
  summarise_all(sum) %>% 
  ungroup %>% 
  mutate(cpd_id = info_x$cpd_id[match(name, info_x$Match)]) %>% 
  dplyr::select(-name) %>% 
  column_to_rownames('cpd_id') %>% 
  profile_transRA()

group <- read.delim('LloydPriceJ_2019.group.txt')

comparisons <- list(c('CD', 'nonIBD'), c('UC', 'nonIBD'))
diffs_1 <- map(comparisons, \(x) {
  group_x <- filter(group, group %in% x)
  profile_x <- dplyr::select(profile, all_of(group_x$sample))
  oplsda <- opls(data.frame(t(profile_x)), y = pull(group_x, group), orthoI = 0)
  diff <- difference_analysis(profile_x, group_x, comparison = x)
  out <- oplsda@vipVn %>%
    data.frame(vip = .) %>% 
    rownames_to_column('name') %>% 
    left_join(diff, ., by = 'name') %>% 
    mutate(enriched = ifelse((vip > 1 | pval < 0.05) & 
                               log2FC > 0, x[1], 
                             ifelse((vip > 1 | pval < 0.05) & 
                                      log2FC < 0, x[2], 'none'))) %>% 
    left_join(info_x, by = c('name' = 'cpd_id'))  } ) %>% 
  set_names(map_vec(comparisons, \(x) paste(x, collapse = '_vs_')))
write.xlsx(diffs_1, 'LloydPriceJ_2019.diff.xlsx')

#### BushmanFD_2020 ####
info <- read.delim('BushmanFD_2020.mbt_info.txt')

info_x <- info %>% 
  filter(Match != '') %>% 
  dplyr::select(-cpd_id, -cpd_name) %>% 
  distinct() %>% 
  add_column(cpd_id = sprintf('cpd_%03d', seq(1, nrow(.))), .before = 1)

group <- read.delim('BushmanFD_2020.group.txt') %>% 
  filter(group %in% c('IBD', 'Healthy')) %>% 
  mutate(group = ifelse(group == 'Healthy', 'HC', group))

profile <- read.delim('BushmanFD_2020.profile.txt', check.names = F, 
                      row.names = 1) %>% 
  apply(2, \(x) 10 ^ x) %>% 
  data.frame(check.names = F) %>% 
  mutate(name = info$Match[match(rownames(.), info$cpd_id)]) %>% 
  filter(name %in% info_x$Match) %>% 
  group_by(name) %>% 
  summarise_all(sum) %>% 
  ungroup %>% 
  mutate(cpd_id = info_x$cpd_id[match(name, info_x$Match)]) %>% 
  dplyr::select(-name) %>% 
  column_to_rownames('cpd_id') %>% 
  dplyr::select(all_of(group$sample))

comparisons <- list(c('IBD','HC'))

diffs_2 <- map(comparisons, \(x) {
  group_x <- filter(group, group %in% x)
  profile_x <- dplyr::select(profile, all_of(group_x$sample))
  plsda <- opls(data.frame(t(profile_x)), y = pull(group_x, group), orthoI = 0)
  diff <- difference_analysis(profile_x, group_x, comparison = x)
  out <- plsda@vipVn %>%
    data.frame(vip = .) %>% 
    rownames_to_column('name') %>% 
    left_join(diff, ., by = 'name') %>% 
    mutate(enriched = ifelse((vip > 1 | pval < 0.05) & 
                               log2FC > 0, x[1], 
                             ifelse((vip > 1 | pval < 0.05) & 
                                      log2FC < 0, x[2], 'none'))) %>% 
    left_join(info_x, by = c('name' = 'cpd_id'))  } ) %>% 
  set_names(map_vec(comparisons, \(x) paste(x, collapse = '_vs_')))
write.xlsx(diffs_2, 'BushmanFD_2020.diff.xlsx')

#### SchirmerM_2024 ####
info <- read.delim('SchirmerM_2024.mbt_info.txt')

info_x <- info %>% 
  filter(Match != '') %>% 
  dplyr::select(-cpd_id, -cpd_name) %>% 
  distinct() %>% 
  add_column(cpd_id = sprintf('cpd_%03d', seq(1, nrow(.))), .before = 1)

profile <- read.delim('SchirmerM_2024.profile.txt', check.names = F) %>% 
  map_df(\(x) replace_na(x, 0)) %>% 
  column_to_rownames('name') %>% 
  mutate(name = info$Match[match(rownames(.), info$cpd_id)]) %>% 
  filter(name %in% info_x$Match) %>% 
  group_by(name) %>% 
  summarise_all(sum) %>% 
  ungroup %>% 
  mutate(cpd_id = info_x$cpd_id[match(name, info_x$Match)]) %>% 
  dplyr::select(-name) %>% 
  column_to_rownames('cpd_id') %>% 
  profile_transRA()

group <- read.delim('SchirmerM_2024.group.txt') %>% 
  filter(sample %in% colnames(profile)) %>% 
  dplyr::select(sample, group = stage) %>% 
  mutate(group = sub('/','.', group))

comparisons <- list(c('mild','inactive'), c('moderate.severe', 'inactive'))

diffs_3 <- map(comparisons, \(x) {
  group_x <- filter(group, group %in% x)
  profile_x <- dplyr::select(profile, all_of(group_x$sample))
  plsda <- opls(data.frame(t(profile_x)), y = pull(group_x, group), orthoI = 0)
  diff <- difference_analysis(profile_x, group_x, comparison = x)
  out <- plsda@vipVn %>%
    data.frame(vip = .) %>% 
    rownames_to_column('name') %>% 
    left_join(diff, ., by = 'name') %>% 
    mutate(enriched = ifelse((vip > 1 | pval < 0.05) & 
                               log2FC > 0, x[1], 
                             ifelse((vip > 1 | pval < 0.05) & 
                                      log2FC < 0, x[2], 'none'))) %>% 
    left_join(info_x, by = c('name' = 'cpd_id'))  } ) %>% 
  set_names(map_vec(comparisons, \(x) paste(x, collapse = '_vs_')))
write.xlsx(diffs_3, 'SchirmerM_2024.diff.xlsx')

#### enrichment ####
bg_data <- read.delim('/database/KEGG/enrichment_analysis/cpd2path_enrichment.tsv')

eKEGGs_1 <- map(diffs_1, ~ 
                  filter(.x, enriched != 'none' & KEGG != '') %>% 
                  pull(KEGG) %>% 
                  enricher(TERM2GENE = bg_data, minGSSize = 1, 
                           pvalueCutoff = 1, qvalueCutoff = 1) %>% 
                  data.frame )
write.xlsx(eKEGGs_1, 'LloydPriceJ_2019.eKEGG.xlsx')

eKEGGs_2 <- map(diffs_2, ~
                  filter(.x, enriched != 'none' & KEGG != '') %>% 
                  pull(KEGG) %>% 
                  enricher(TERM2GENE = bg_data, minGSSize = 1, 
                           pvalueCutoff = 1, qvalueCutoff = 1) %>% 
                  data.frame )
write.xlsx(eKEGGs_2, 'BushmanFD_2020.eKEGG.xlsx')

eKEGGs_3 <- map(diffs_3, ~
                  filter(.x, enriched != 'none' & KEGG != '') %>% 
                  pull(KEGG) %>% 
                  enricher(TERM2GENE = bg_data, minGSSize = 1, 
                           pvalueCutoff = 1, qvalueCutoff = 1) %>% 
                  data.frame )
write.xlsx(eKEGGs_3, 'SchirmerM_2024.eKEGG.xlsx')

#### intersect - metabolite ####
LloydPriceJ_2019.CD_vs_nonIBD = diffs_1$CD_vs_nonIBD %>% 
  filter(HMDB != '' & enriched != 'none')
LloydPriceJ_2019.UC_vs_nonIBD = diffs_1$UC_vs_nonIBD %>% 
  filter(HMDB != '' & enriched != 'none')
BushmanFD_2020.IBD_vs_HC = diffs_2$IBD_vs_HC %>% 
  filter(HMDB != '' & enriched != 'none')
SchirmerM_2024.moderate.severe_vs_inactive = diffs_3$moderate.severe_vs_inactive %>% 
  filter(HMDB != '' & enriched != 'none')

mbt <- list(LloydPriceJ_2019.CD_vs_nonIBD = 
              LloydPriceJ_2019.CD_vs_nonIBD %>% pull(HMDB),
            LloydPriceJ_2019.UC_vs_nonIBD = 
              LloydPriceJ_2019.UC_vs_nonIBD %>% pull(HMDB),
            BushmanFD_2020.IBD_vs_HC = 
              BushmanFD_2020.IBD_vs_HC %>% pull(HMDB),
            SchirmerM_2024.moderate.severe_vs_inactive = 
              SchirmerM_2024.moderate.severe_vs_inactive %>% pull(HMDB))

ggVennDiagram::ggVennDiagram(mbt, label = 'both', label_color = 'black', label_alpha = 0,
                             set_color = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3'),
                             edge_lty = 'dashed', edge_size = 1) +
  scale_fill_gradient(low='#fef5f5',high = '#f0ad80',
                      name = 'Number of differential\ncontent metabolites')
ggsave('venn.pdf', width = 5, height = 5)

mbt_info <- purrr::reduce(mbt,\(x, y) intersect(x, y) )  %>% 
  data.frame(HMDB = .) %>% 
  left_join(dplyr::select(info_x, -cpd_id) %>% distinct(), by = 'HMDB') %>% 
  left_join(dplyr::select(LloydPriceJ_2019.CD_vs_nonIBD, HMDB, 
                          LloydPriceJ_2019.CD_vs_nonIBD = log2FC), by = 'HMDB') %>% 
  left_join(dplyr::select(LloydPriceJ_2019.UC_vs_nonIBD, HMDB, 
                          LloydPriceJ_2019.UC_vs_nonIBD = log2FC), by = 'HMDB') %>% 
  left_join(dplyr::select(BushmanFD_2020.IBD_vs_HC, HMDB, 
                          BushmanFD_2020.IBD_vs_HC = log2FC), by = 'HMDB') %>% 
  left_join(dplyr::select(SchirmerM_2024.moderate.severe_vs_inactive, HMDB, 
                          SchirmerM_2024.moderate.severe_vs_inactive = log2FC), by = 'HMDB') %>% 
  mutate(LloydPriceJ_2019.CD_vs_nonIBD = ifelse(HMDB == 'HMDB0000292', 
                                                LloydPriceJ_2019.CD_vs_nonIBD + 12,
                                                LloydPriceJ_2019.CD_vs_nonIBD))
plot_data <- mbt_info %>% 
  dplyr::select(Match, contains('_vs_')) %>% 
  gather(key = 'group', value = 'log2FC', -Match)

tmp <- plot_data %>% 
  group_by(Match) %>% 
  group_modify(~.x$log2FC %>% 
                 data.frame(greater = sum(. > 0), less = sum(. < 0))) %>% 
  filter(greater == 4 | less == 4 ) %>% 
  pull(Match) %>% 
  unique()

plot_data <- plot_data %>% 
  filter(Match %in% tmp)

t_test <- plot_data %>% 
  group_by(Match) %>% 
  group_modify(~{
    data_x = .x
    if (median(data_x$log2FC) > 0) {
      t = t.test(data_x$log2FC, mu = 0, alternative = 'greater')
      data.frame(t = t$statistic, 
                 pval = t$p.value, 
                 conf = paste0(t$conf.int, collapse = ','), 
                 alter = t$alternative)
    } else {
      t = t.test(data_x$log2FC, mu = 0, alternative = 'less')
      data.frame(t = t$statistic, 
                 pval = t$p.value, 
                 conf = paste0(t$conf.int, collapse = ','), 
                 alter = t$alternative)
    }} ) %>% 
  filter(pval < 0.05) %>% 
  add_plab(format = 1) %>% 
  arrange(desc(t)) %>% 
  left_join(dplyr::select(mbt_info, Match, HMDB), by = 'Match')
write.xlsx(t_test, 'one-sample t-test.xlsx')

plot_data_x <- plot_data %>% 
  left_join(dplyr::select(mbt_info, Match, HMDB), by = 'Match') %>% 
  filter(HMDB %in% t_test$HMDB) %>%
  mutate(HMDB = factor(HMDB, t_test$HMDB),
         HMDB = forcats::fct_relevel(HMDB, 'HMDB0002302', after = Inf))

ggscatter(plot_data_x, 'HMDB', 'log2FC', color = 'group', rotate = T, 
          palette = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3'),
          xlab = '', ylab = 'Log2-FoldChange') +
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') +
  geom_text(aes(y = -2.5, x = HMDB, label = plab), color = 'red', 
            data = t_test, inherit.aes = F) +
  theme_pubr() +
  theme(panel.grid.major = element_line(color = 'grey88', linewidth = .5),
        aspect.ratio = 1.5)
ggsave('one-sample t-test.pdf', width = 5, height = 4.8)

#### intersect - path ####
kegg_info <- read.delim('/database/KEGG/v20230401/path_description') %>% 
  filter(category == 'Metabolism')

plot_data <- rbind(eKEGGs_1$CD_vs_nonIBD %>% 
                     data.frame() %>% 
                     dplyr::select(ID, pvalue, FoldEnrichment) %>% 
                     add_column(group = 'LloydPriceJ_2019.CD_vs_nonIBD'),
                   LloydPriceJ_2019.UC_vs_nonIBD = eKEGGs_1$UC_vs_nonIBD %>% 
                     data.frame() %>% 
                     dplyr::select(ID, pvalue, FoldEnrichment) %>% 
                     add_column(group = 'LloydPriceJ_2019.UC_vs_nonIBD'),
                   BushmanFD_2020.IBD_vs_HC = eKEGGs_2$IBD_vs_HC %>% 
                     data.frame() %>% 
                     dplyr::select(ID, pvalue, FoldEnrichment) %>% 
                     add_column(group = 'BushmanFD_2020.IBD_vs_HC'),
                   SchirmerM_2024.moderate.severe_vs_inactive = eKEGGs_3$moderate.severe_vs_inactive %>% 
                     data.frame() %>% 
                     dplyr::select(ID, pvalue, FoldEnrichment) %>% 
                     add_column(group = 'SchirmerM_2024.moderate.severe_vs_inactive') ) %>% 
  mutate(KEGG = map_vec(strsplit(ID, ':'), \(x) x[1])) %>% 
  filter(KEGG %in% kegg_info$name)

tmp <- plot_data %>% 
  group_by(ID) %>% 
  group_modify(~data.frame(n = nrow(.x))) %>% 
  filter(n >=2 ) %>% 
  pull(ID) %>% 
  unique()

tmp <- c('map00966','map00780','map00730','map00680','map00650',
         'map00480','map00400','map00380','map00350','map00230',
         'map00300','map00290','map00061','map00020','map00620') %>% 
  map_vec(\(x) plot_data$ID[base::grepl(x, plot_data$ID)] %>% unique)

plot_data_x <- plot_data %>% 
  filter(ID %in% tmp) %>% 
  filter(pvalue < 0.1)

ggplot(plot_data_x, aes(group, ID)) +
  geom_point(aes(fill = pvalue, size = FoldEnrichment), shape = 21) +
  scale_size_continuous(range = c(2, 7)) +
  scale_fill_viridis_c(begin = .5) +
  labs(x = '', y = '') +
  theme_pubr() +
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_line(color = 'grey88', linewidth = .5),
        legend.position = 'right',
        aspect.ratio = 2)
ggsave('eKEGG.pdf', width = 9, height = 9)
