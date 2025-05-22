#### Jinxin Meng, 20240801, 2025325 ####
# fastp - hisat2 - stringtie
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, openxlsx)
pacman::p_load(limma, ComplexHeatmap, clusterProfiler, enrichplot)
source('/Code/R_func/profile_process.R')

#### info ####
group <- read.delim('group.txt')

rc <- read.delim('rc_symbol.tsv', row.names = 1)
fpkm <- read.delim('fpkm_symbol.tsv', row.names = 1)

#### Heatmap ####
markers <- list(
  aISC = c('Lgr5','Ascl2','Slc12a2'),
  rISC = c('Hopx','Bmi1','Lrig1'),
  PC = c('Ang4','Lyz1','Agr2'),
  GC = c('Muc2','Ccl9'),
  TC = c('Dclk1','Cd24a'),
  EEC = c('Chga','Neurod1')
)

marker_levels <- c('aISC', 'rISC', 'PC', 'GC', 'TC', 'EEC')

# D0
group_x <- filter(group, group %in% c('PBS_D0', 'IPA_D0'))

fpkm_x <- dplyr::select(fpkm, all_of(group_x$sample)) %>%
  filter(rownames(.) %in% unlist(markers)) %>% 
  profile_transSqrt()

row_data <- markers %>% 
  map2_dfr(names(.), \(x, y) data.frame(name = x, class = y)) %>% 
  arrange(match(name, rownames(fpkm_x))) %>% 
  mutate(class = factor(class, marker_levels))

col_data <-  group_x

row_split <- row_data$class
col_split <- col_data$group

pdf('fpkm_D0.pdf', width = 5, height = 5)
pheatmap(fpkm_x, scale = 'row',
         color = colorRampPalette(c('#307cc0', 'white', '#e43589'))(100),
         split = row_split, column_split = col_split,
         cluster_row_slices = F, cluster_column_slices = F, 
         row_title_rot = 0,
         row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
         row_gap = unit(1, 'mm'), column_gap = unit(1, 'mm'),
         treeheight_col = 10, treeheight_row = 10,
         cellheight = 12, cellwidth = 12, fontsize = 10, 
         border_color = 'white', border_gp = gpar(col = 'black'), 
         heatmap_legend_param = list(title = 'Scale FPKM') )
dev.off()

# D3 
group_x <- filter(group, group %in% c('PBS_D3', 'IPA_D3'))

fpkm_x <- dplyr::select(fpkm, all_of(group_x$sample)) %>%
  filter(rownames(.) %in% unlist(markers))

row_data <- markers %>% 
  map2_dfr(names(.), \(x, y) data.frame(name = x, class = y)) %>% 
  arrange(match(name, rownames(fpkm_x))) %>% 
  mutate(class = factor(class, marker_levels))

col_data <-  group_x

row_split <- row_data$class
col_split <- col_data$group

pdf('fpkm_D3.pdf', width = 5, height = 5)
pheatmap(fpkm_x, scale = 'row',
         color = colorRampPalette(c('#307cc0', 'white', '#e43589'))(100),
         split = row_split, column_split = col_split,
         cluster_row_slices = F, cluster_column_slices = F, 
         row_title_rot = 0,
         row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
         row_gap = unit(1, 'mm'), column_gap = unit(1, 'mm'),
         treeheight_col = 10, treeheight_row = 10,
         cellheight = 12, cellwidth = 12, fontsize = 10, 
         border_color = 'white', border_gp = gpar(col = 'black'), 
         heatmap_legend_param = list(title = 'Scale FPKM') )
dev.off()

# D7 
group_x <- filter(group, group %in% c('PBS_D7', 'IPA_D7'))

fpkm_x <- dplyr::select(fpkm, all_of(group_x$sample)) %>%
  filter(rownames(.) %in% unlist(markers))

row_data <- markers %>% 
  map2_dfr(names(.), \(x, y) data.frame(name = x, class = y)) %>% 
  arrange(match(name, rownames(fpkm_x))) %>% 
  mutate(class = factor(class, marker_levels))

col_data <-  group_x

row_split <- row_data$class
col_split <- col_data$group

pdf('fpkm_D7.pdf', width = 5, height = 5)
pheatmap(fpkm_x, scale = 'row',
         color = colorRampPalette(c('#307cc0', 'white', '#e43589'))(100),
         split = row_split, column_split = col_split,
         cluster_row_slices = F, cluster_column_slices = F, 
         row_title_rot = 0,
         row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
         row_gap = unit(1, 'mm'), column_gap = unit(1, 'mm'),
         treeheight_col = 10, treeheight_row = 10,
         cellheight = 12, cellwidth = 12, fontsize = 10, 
         border_color = 'white', border_gp = gpar(col = 'black'), 
         heatmap_legend_param = list(title = 'Scale FPKM') )
dev.off()

#### diff ####
library(org.Mm.eg.db)

gene_info <- AnnotationDbi::select(org.Mm.eg.db, rownames(fpkm), 
                                   keytype = 'SYMBOL', columns = c('ENTREZID'))

fpkm_x <- profile_transLOG2(fpkm)

metadata <- data.frame(row.names = colnames(fpkm_x)) %>%
  mutate(group = group$group[match(rownames(.), group$sample)],
         group = factor(group))

grps <- metadata$group
design <- model.matrix(~ 0 + grps)
colnames(design) <- gsub('grps', '', colnames(design))

comparisons <- list(c('IPA_D0', 'PBS_D0'), 
                    c('IPA_D3', 'PBS_D3'), 
                    c('IPA_D7', 'PBS_D7'))

contrasts <- makeContrasts(contrasts = map_vec(comparisons, 
                                               ~ paste(.x, collapse = '-')), 
                           levels = design)

fit <- lmFit(fpkm_x, design)
fit_contrast <- contrasts.fit(fit, contrasts)
fit_contrast <- eBayes(fit_contrast)

diffs <- map(comparisons, \(x) 
             topTable(fit_contrast, adjust = 'BH', number = Inf, 
                      coef = paste0(x, collapse = '-') ) %>% 
               rownames_to_column('name') %>% 
               mutate(enriched = ifelse(logFC > .5 & P.Value < 0.05, x[1], 
                                        ifelse(logFC < -.5 & P.Value < 0.05, 
                                               x[2], 'none'))) %>% 
                 left_join(gene_info, by = c('name' = 'SYMBOL'))) %>% 
  set_names(map_vec(comparisons, ~ paste(.x, collapse = '_vs_')) )
write.xlsx(diffs, 'diffs.xlsx')

#### ORA KEGG ####
R.utils::setOption('clusterProfiler.download.method', 'wget')

eKEGGs <- map(diffs, \(x)
              filter(x, enriched != 'none') %>% 
               pull(ENTREZID) %>% 
               na.omit %>% 
               enrichKEGG(organism = 'mmu', pvalueCutoff = 1, 
                          qvalueCutoff = 1) %>% 
               data.frame())
write.xlsx(eKEGGs, 'eKEGGs.xlsx')

eKEGGs <- map2(diffs, map_vec(comparisons, ~ .x[1]), \(x, y)
              filter(x, enriched == y) %>% 
                pull(ENTREZID) %>% 
                na.omit %>% 
                enrichKEGG(organism = 'mmu', pvalueCutoff = 1, 
                           qvalueCutoff = 1) %>% 
                data.frame())
write.xlsx(eKEGGs, 'eKEGGs_up.xlsx')

eKEGGs <- map2(diffs, map_vec(comparisons, ~ .x[2]), \(x, y)
               filter(x, enriched == y) %>% 
                pull(ENTREZID) %>% 
                na.omit %>% 
                enrichKEGG(organism = 'mmu', pvalueCutoff = 1, 
                           qvalueCutoff = 1) %>% 
                data.frame())
write.xlsx(eKEGGs, 'eKEGGs_down.xlsx')

#### ORA GO ####
eGOs <- map(diffs, \(x)
            filter(x, enriched != 'none') %>% 
              pull(name) %>% 
              enrichGO(org.Mm.eg.db, keyType = 'SYMBOL', ont = 'ALL',
                       pvalueCutoff = 1, qvalueCutoff = 1) %>% 
              data.frame()) 
write.xlsx(eGOs, 'eGOs.xlsx')

eGOs <- map2(diffs, map_vec(comparisons, ~ .x[1]), \(x, y)
            filter(x, enriched == y) %>% 
              pull(name) %>% 
              enrichGO('org.Mm.eg.db', keyType = 'SYMBOL', ont = 'ALL',
                       pvalueCutoff = 1, qvalueCutoff = 1) %>% 
              data.frame())
write.xlsx(eGOs, 'eGOs_up.xlsx')

eGOs <- map2(diffs, map_vec(comparisons, ~ .x[2]), \(x, y)
             filter(x, enriched == y) %>% 
              pull(name) %>% 
              enrichGO('org.Mm.eg.db', keyType = 'SYMBOL', ont = 'ALL',
                       pvalueCutoff = 1, qvalueCutoff = 1) %>% 
              data.frame())
write.xlsx(eGOs, 'eGOs_down.xlsx')

#### ORA custom plot ####
data <- read.delim('up-gene enriched terms.txt')

data %>% 
  filter(group == 'D3') %>% 
  dplyr::select(name = Description, FoldEnrichment, pvalue) %>% 
  arrange(FoldEnrichment) %>% 
  mutate(name = factor(name, name)) %>% 
  ggscatter('FoldEnrichment', 'name', fill = 'pvalue', size = 'FoldEnrichment', 
            shape = 21, xlab = 'Fold enrichment', ylab = '', legend = 'right', 
            title = 'KEGG and GO enrichment analysis of up-\nregulated genes at D3 post-DSS treatment') +
  scale_fill_viridis_c(begin = .5) +
  theme(aspect.ratio = 2,
        plot.title = element_text(vjust = .5, face = 'bold'))
ggsave('enrichment_D3.pdf', width = 7, height = 5)

data %>% 
  filter(group == 'D7') %>% 
  dplyr::select(name = Description, FoldEnrichment, pvalue) %>% 
  arrange(FoldEnrichment) %>% 
  mutate(name = factor(name, name)) %>% 
  ggscatter('FoldEnrichment', 'name', fill = 'pvalue', size = 'FoldEnrichment', 
            shape = 21, xlab = 'Fold enrichment', ylab = '', legend = 'right', 
            title = 'KEGG and GO enrichment analysis of up-\nregulated genes at D7 post-DSS treatment') +
  scale_fill_viridis_c(begin = .5) +
  theme(aspect.ratio = 2,
        plot.title = element_text(vjust = .5, face = 'bold'))
ggsave('enrichment_D7.pdf', width = 7, height = 5)

