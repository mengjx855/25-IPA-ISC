#### Jinxin Meng, 20250222, 20250325 ####
pacman::p_load(dplyr, tidyr, tibble, purrr, stringr, ggplot2, ggpubr, openxlsx)
pacman::p_load(Seurat, harmony, DoubletFinder, data.table, cowplot)
pacman::p_load(clusterProfiler, enrichplot)
source('/code/R_func/scRNA_utilities.R')

#### 原始数据 ####
path <- list.files('rawdata/', '1', full.names = T)
name <- list.files('rawdata/', '1')

seurat_list <- map2(path, name, \(x, y) {
  message(paste0('Processing: ', y)) 
  count <- Read10X(x)
  CreateSeuratObject(counts = count, project = y, min.cells = 5, min.features = 100)
})

seurat <- merge(seurat_list[[1]], seurat_list[-1], add.cell.ids = name) %>% 
  JoinLayers()

saveRDS(seurat, 'seurat.raw.rds')

# table(seurat@meta.data$orig.ident)
# IPA_1 Veh_1 
# 10049  8865

#### 质控 ####
mit_genes <- rownames(seurat)[grep('^mt-', rownames(seurat), ignore.case = T)]
seurat <- PercentageFeatureSet(seurat, features = mit_genes, col.name = 'percent_mit')
rib_genes <- rownames(seurat)[grep('^Rp[sl]', rownames(seurat), ignore.case = T)]
seurat <- PercentageFeatureSet(seurat, features = rib_genes, col.name = 'percent_rib')
hb_genes <- rownames(seurat)[grep('^Hb[^(p)]', rownames(seurat), ignore.case = T)]
seurat <- PercentageFeatureSet(seurat,  features = hb_genes, col.name = 'percent_hb')

p1 <- VlnPlot(seurat, group.by = 'orig.ident', features = c('nFeature_RNA', 'nCount_RNA'), 
              pt.size = 0, ncol = 2) + 
  NoLegend()

p2 <- VlnPlot(seurat, group.by = 'orig.ident', features = c('percent_mit', 'percent_rib', 'percent_hb'), 
              pt.size = 0, ncol = 3) + 
  scale_y_continuous(breaks = seq(0, 100, 5)) + 
  NoLegend()

cowplot::plot_grid(p1, p2, ncol = 1, align = 'v')
ggsave('qc.before.pdf', width = 6, height = 10)

# 过滤, 基因表达太多了，可能就是异常点，存在双细胞混一起的可能
cells <- purrr::reduce(list(WhichCells(seurat, expression = percent_mit < 20), 
                            WhichCells(seurat, expression = nFeature_RNA > 100),
                            WhichCells(seurat, expression = nFeature_RNA < 5000),
                            WhichCells(seurat, expression = percent_hb < 1)),
                       \(x, y) base::intersect(x, y))
seurat <- subset(seurat, cells = cells)

p1 <- VlnPlot(seurat, group.by = 'orig.ident', features = c('nFeature_RNA', 'nCount_RNA'), 
              pt.size = 0, ncol = 2) + 
  NoLegend()

p2 <- VlnPlot(seurat, group.by = 'orig.ident', features = c('percent_mit', 'percent_rib', 'percent_hb'), 
              pt.size = 0, ncol = 3) + 
  NoLegend()

cowplot::plot_grid(p1, p2, ncol = 1, align = 'v')
ggsave('qc.after.pdf', width = 6, height = 10)

# dim(seurat)
# [1] 17396  8206

#### 降维聚类 #### 
seurat <- NormalizeData(seurat, normalization.method = 'LogNormalize', scale.factor = 1e4, verbose = T) # 标准化
seurat <- FindVariableFeatures(seurat, selection.method = 'vst', nfeatures = 2000) # 筛选高变基因
VariableFeaturePlot(seurat, cols = c('grey','red'))
seurat <- ScaleData(seurat) # 数据归一化，这步骤有一步回归的分析，去除噪声
seurat <- RunPCA(seurat, features = VariableFeatures(seurat)) # PCA线性降维分析
ElbowPlot(seurat, ndims = 30) # 流石图协助选择PC维度 10
seurat <- RunHarmony(seurat, 'orig.ident') # Harmony去批次
seurat <- RunUMAP(seurat, dims = 1:10, reduction = 'harmony')
seurat <- RunTSNE(seurat, dims = 1:10, reduction = 'harmony')
seurat <- FindNeighbors(seurat, reduction = 'harmony', dims = 1:10)
seurat <- FindClusters(seurat, resolution = c(seq(0, 1.9, .2)))

map2(paste0('RNA_snn_res.', seq(0, 1.9, .2)), paste0('SNN: ', seq(0, 1.9, .2)),
     \(x, y) DimPlot(seurat, reduction = 'umap', group.by = x, label = T) +
       ggtitle(y) + 
       theme(aspect.ratio = 1) + 
       guides(color = 'none')) %>% 
  plot_grid(plotlist = ., nrow = 2)
ggsave('umap.SNN.dimplot.jpg', width = 20, height = 8)

#### 注释 ####
markers <- list(
  Stem_cells = c('Lgr5','Hopx','Bmi1','Ascl2','Olfm4','Smoc2'),
  Paneth_cells = c('Ang4','Lyz1','Defa17','Defa21','Defa22','Defa24','Gm14851','Defa30'),
  Goblet_cells = c('Muc2','Ccl9','Clca1','Tff3','Agr2','Fcgbp','Zg16'),
  Enteroendocrine_cells = c('Chga','Neurog3','Chgb','Tac1','Tph1'),
  Tuft_cells = c('Dclk1','Trpm5'),
  Enterocyte = c('Aldob', 'Apoa1', 'Apoa4', 'Gsta1','Fabp1', 'Prap1'),
  T_cells = c('Cd3g','Cd4','Cd8a','Tcrg-C1','Tcrg-C2','Tcrg-C4','Gzma','Gzmb'),
  B_cells = c('Cd79a','Jchain','Igha','Igkc'))

cell_col <- c(Stem_cells = '#8dd3c7', Paneth_cells = '#ffed6f', Goblet_cells = '#bebada',
              Enteroendocrine_cells = '#fb8072', Tuft_cells = '#80b1d3',
              Enterocyte = '#fdb462', T_cells = '#b3de69', B_cells = '#bc80bd')

markers <- map(markers, \(x) x[x %in% rownames(seurat)])

# featureplot
max_len <- map_vec(markers, \(x) length(x)) %>% max()
map2(markers, names(markers), \(x, y) {
  p_text <- list(ggtext(data.frame(x = 1, y = 1), 'x', 'y', 
                        label = stringr::str_to_title(gsub('_', ' ', y)), 
                        size = 14, face = 'bold') + 
                   theme_void())
  
  p_feats <- map(x, \(m) FeaturePlot(seurat, reduction = 'umap', features = m) +
                   theme(aspect.ratio = 1))
  
  p_void <- list(ggtext(data.frame(x = 1, y = 1), 'x', 'y', size = 14, face = 'bold') + 
                   theme_void())
  
  cowplot::plot_grid(plotlist = c(p_text, p_feats, 
                                  rep(p_void, times = max_len - length(x))), nrow = 1)
  } ) %>% 
  cowplot::plot_grid(plotlist = ., align = 'v', ncol = 1)
ggsave('umap.featureplot.jpg', width = 4 * max_len + 4, height = 4 * length(markers))

seurat <- SetIdent(seurat, value = 'RNA_snn_res.1.4')
DimPlot(seurat, reduction = 'tsne', label = T, alpha = 1, raster = F) + 
  NoLegend() + 
  theme(aspect.ratio = 1)
ggsave('tsne.SNN.1.4.dimplot.pdf', width = 6, height = 6)

# 鉴定每个cluster中的marker
allmarkers <- FindAllMarkers(seurat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(allmarkers, 'cluster.allmarkers.rds')
write.xlsx(allmarkers, 'cluster.allmarkers.xlsx')

map(markers, \(x) filter(allmarkers, gene %in% x) ) %>% 
  keep(\(x) nrow(x)!=0) %>% 
  map(\(x) group_by(x, cluster) %>% 
        group_modify(~data.frame(match_n = nrow(.x), 
                                 genes = paste(.x$gene, collapse = ',')))) %>% 
  map2_df(names(.), \(x, y) add_column(x, name = y, .before = 1)) %>% 
  data.frame() %>% 
  arrange(cluster)

seurat <- RenameIdents(seurat, 
                       '0' = 'Goblet_cells','1' = 'Paneth_cells','2' = 'Enterocyte',
                       '3' = 'T_cells','4' = 'Enterocyte','5' = 'Enterocyte',
                       '6' = 'Goblet_cells','7' = 'Enterocyte', '8' = 'Goblet_cells',
                       '9' = 'Stem_cells','10' = 'Enterocyte','11' = 'Goblet_cells', 
                       '12' = 'Enterocyte','13' = 'Enterocyte','14' = 'Paneth_cells',
                       '15' = 'Enterocyte', '16' = 'Enterocyte','17' = 'Enterocyte',
                       '18' = 'Paneth_cells','19' = 'Enterocyte','20' = 'Goblet_cells',
                       '21' = 'T_cells','22' = 'Goblet_cells','23' = 'Tuft_cells',
                       '24' = 'Enteroendocrine_cells','25' = 'Goblet_cells','26' = 'T_cells',
                       '27' = 'B_cells','28' = 'Enterocyte')

# # 疑难细胞群注释
# allmarkers %>%
#   filter(p_val_adj < 0.001 & cluster == '1') %>%
#   pull() %>%
#   sc_marker_match(tissue = 'all_tissue', org = 'mouse') %>%
#   select(-species, -tissue, -cell, -gene) %>%
#   filter(grepl('intestine', item)) %>%
#   data.frame() %>%
#   head(10)

# 可视化
DimPlot(seurat, reduction = 'tsne', label = T, alpha = 1, cols = cell_col) + 
  theme(aspect.ratio = 1) + NoLegend()
ggsave('tsne.marker.dimplot.pdf', width = 6, height = 6)

DimPlot(seurat, reduction = 'tsne', label = T, alpha = 1, group.by = 'orig.ident') + 
  theme(aspect.ratio = 1)
ggsave('tsne.marker.dimplot.group.pdf', width = 6, height = 5)

DimPlot(seurat, reduction = 'tsne', label = T, alpha = 1, 
        split.by = 'orig.ident', cols = cell_col) + 
  theme(aspect.ratio = 1)
ggsave('tsne.marker.dimplot.split.pdf', width = 10, height = 5)

# marker
allmarkers <- FindAllMarkers(seurat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(allmarkers, 'cell.allmarkers.rds')
write.xlsx(allmarkers, 'cell.allmarkers.xlsx')

markers <- list(
  Stem_cells = c('Lgr5','Bmi1','Ascl2','Smoc2'),
  Paneth_cells = c('Ang4','Lyz1'),
  Goblet_cells = c('Muc2','Clca1','Tff3','Agr2'),
  Enteroendocrine_cells = c('Chga','Tac1','Tph1'),
  Tuft_cells = c('Dclk1','Trpm5'),
  Enterocyte = c('Aldob', 'Gsta1', 'Prap1'),
  T_cells = c('Cd3g','Cd8a','Gzma'),
  B_cells = c('Cd79a','Igkc'))

top_markers <- map2_df(markers, names(markers), \(x, y) 
                       rbind(filter(allmarkers, cluster == y & gene %in% x),
                             filter(allmarkers, cluster == y & !gene %in% x) %>% 
                               head(n = 4 - sum(allmarkers$cluster == y & allmarkers$gene %in% x)) ) %>% 
                         arrange(desc(avg_log2FC)) )

plot_data <- FetchData(seurat, vars = unique(top_markers$gene)) %>% 
  add_column(cell = seurat@active.ident) %>% 
  gather(key = 'gene', value = 'value', -cell) %>% 
  mutate(gene = factor(gene, unique(top_markers$gene)),
         cell = factor(cell, names(markers)))

ggplot(plot_data, aes(cell, value), color = factor(cell)) +
  geom_violin(aes(fill = cell), scale ='width') +
  facet_grid(gene ~ ., scales = 'free_y') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cell_col) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x.bottom = element_text(angle = 45, hjust = 1,
                                          vjust = NULL, color = 'black', size = 14),
        axis.text.y.left = element_blank(),
        legend.position = 'none',
        panel.spacing.y = unit(0, 'cm'),
        strip.text.y = element_text(angle = 0, size = 14, hjust = 0, face = 'italic'),
        strip.background.y = element_blank())
ggsave('cell.allmarkers.pdf', width = 5, height = 9)

saveRDS(seurat, 'seurat.out.rds')

#### 细胞组成 ####
data <- seurat@meta.data %>% 
  dplyr::select(group = orig.ident) %>% 
  add_column(cells = Idents(seurat)) %>% 
  group_by(group, cells) %>% 
  summarise(value = n()) %>% 
  group_modify(~.x %>% mutate(prec = value / sum(value) * 100)) %>% 
  ungroup

ggbarplot(data, 'group', 'prec', fill = 'cells', xlab = '', 
          ylab = 'Cell precentage (%)', legend = 'right', palette = cell_col)
ggsave('cell_prec.barplot.pdf', width = 4, height = 4)

# 干细胞柱状图
# group value
# IPA_1  4825
# Veh_1  3381
#      IPA   Veh
# Yes  234   102
# No   4591  3279
chisq.test(matrix(c(234, 102, 4591, 3279), nrow = 2))

filter(data, cells == 'Stem_cells') %>% 
  mutate(group = factor(group, c('Veh_1', 'IPA_1'))) %>% 
  ggbarplot('group', 'prec', fill ='group', palette = c('#a09fa4','#87d6f5'),
            xlab = '', ylab = 'Percentage of CBCs (%)') +
  geom_signif(annotations = '***', y_position = 5.2, xmin = 1, xmax = 2,
              tip_length = 0.1, size = .5) +
  theme(aspect.ratio = 1.5)

# 所有细胞
filter(data, !cells %in% c('T_cells','B_cells')) %>% 
  dplyr::select(-value) %>% 
  spread('group', 'prec') %>% 
  mutate(prop = IPA_1 / Veh_1,
         prop = round(prop, 2))

filter(data, !cells %in% c('T_cells','B_cells')) %>% 
  select(-prec) %>% 
  spread('group', 'value') %>% 
  mutate(IPA_other = 4825 - IPA_1,
         Veh_total = 3381 - Veh_1) %>% 
  group_by(cells) %>% 
  group_modify(~ chisq.test(matrix(unlist(.x), 2))$p.value %>% 
                 data.frame(pval = .))

filter(data, !cells %in% c('T_cells','B_cells')) %>%
  mutate(group = factor(group, c('IPA_1', 'Veh_1')),
         cells = factor(cells, c('Enterocyte','Goblet_cells','Paneth_cells',
                                 'Stem_cells','Tuft_cells','Enteroendocrine_cells'))) %>% 
  ggbarplot('cells', 'prec', fill = 'group', palette = c('#87d6f5', '#a09fa4'),
            xlab = '', ylab = 'Percentage (%)') +
  annotate('text', x = 1, y = 85, label = '1.34-Fold ***', size = 4) +
  annotate('text', x = 2, y = 60, label = '0.76-Fold ***', size = 4) +
  annotate('text', x = 3, y = 30, label = '0.92-Fold ns', size = 4) +
  annotate('text', x = 4, y = 9, label = '1.61-Fold ***', size = 4) +
  annotate('text', x = 5, y = 6, label = '0.96-Fold ns', size = 4) +
  annotate('text', x = 6, y = 3, label = '0.91-Fold ns', size = 4) +
  theme(aspect.ratio = 3/3.5)
ggsave('cell_prec.cell.barplot.pdf', width = 5, height = 4)

#### 各种类型细胞中 hmgcs2 的表达 ####
# dimplot
data <- FetchData(seurat, vars = c('Hmgcs2')) %>% 
  add_column(cells = seurat@active.ident) %>% 
  rownames_to_column('name') %>% 
  left_join(seurat@reductions$tsne@cell.embeddings %>% 
              data.frame %>% 
              rownames_to_column('name'),
            by = 'name') %>% 
  mutate(.Hmgcs2 = ifelse(Hmgcs2 > 0, 1, 0))

ggscatter(data, x = 'tSNE_1', 'tSNE_2', color = 'cells', alpha = '.Hmgcs2', 
          size = .3, palette = cell_col) +
  scale_alpha_continuous(range = c(.05, 1)) +
  theme(aspect.ratio = 1)

# accompany barplot
data <- FetchData(seurat, vars = c('Hmgcs2')) %>% 
  add_column(cells = seurat@active.ident,
             group = seurat@meta.data$orig.ident) %>% 
  rownames_to_column('name') %>% 
  mutate(value = ifelse(Hmgcs2 > 0, 'pos', 'neg')) %>% 
  count(group, value) %>% 
  group_by(group) %>% 
  group_modify(~ .x %>% mutate(total = sum(n))) %>% 
  mutate(perc = n / total * 100) %>% 
  filter(value == 'pos')

# 44.6 / 30.5
chisq.test(matrix(c(2154, 2671, 1032, 2349), nrow = 2))

ggbarplot(data, 'value', 'perc', fill = 'group', xlab = '', 
          palette = c('#87d6f5','#a09fa4')) +
  scale_y_continuous(expand = c(.01, .01)) +
  annotate('text', x = 1, y = 60, label = '44.6%', size = 4) +
  annotate('text', x = 1, y = 20, label = '30.5%', size = 4) +
  annotate('text', x = 1, y = 70, label = '1.46-Fold ***', size = 4) +
  theme(aspect.ratio = 5)

# barplot 
data <- FetchData(seurat, vars = c('Hmgcs2')) %>% 
  add_column(cells = seurat@active.ident,
             group = seurat@meta.data$orig.ident) %>% 
  rownames_to_column('name') %>% 
  mutate(value = ifelse(Hmgcs2 > 0, 'pos', 'neg')) %>% 
  count(cells, value) %>% 
  group_by(cells) %>% 
  group_modify(~ .x %>% mutate(total = sum(n))) %>% 
  mutate(perc = n / total * 100) %>% 
  filter(value == 'pos')

ggbarplot(data, 'cells', 'perc', fill = 'cells', legend = 'none', xlab = '', width = .6,
          ylab = 'Precentage (%)', sort.val = 'asc', sort.by.groups = F, rotate = T,
          title = 'The percentage of Hmgcs2+\ncells in various cell types',
          palette = cell_col, label = paste0(signif(data$perc, 3), '%'), 
          lab.vjust = .5, lab.hjust = -.1) +
  theme(aspect.ratio = 4/3)

# barplot with fold
diffs <- readRDS('cells.diff.rds')

plot_data <- map2_dfr(diffs, names(diffs), ~
                        filter(.x, name == 'Hmgcs2') %>% 
                        add_column(cell = .y)) %>% 
  filter(! cell %in% c('T_cells', 'B_cells')) %>% 
  dplyr::select(pct.1, pct.2, cell, p_val) %>% 
  gather('pct', 'value', -cell, -p_val) %>% 
  mutate(group = ifelse(pct == 'pct.1', 'IPA', 'Veh'),
         group = factor(group, c('Veh', 'IPA')))

plot_data %>% 
  dplyr::select(cell, group, value) %>% 
  spread('group', 'value') %>% 
  mutate(perc = IPA / Veh,
         perc = round(perc, 2))

ggbarplot(plot_data, 'cell', 'value', fill = 'group', xlab = '', 
          palette = c('#a09fa4','#87d6f5'),
          position = position_dodge(width = .85), 
          ylab = 'Percentage of Hmgcs2+ cell') +
  scale_y_continuous(expand = c(.01, .01)) +
  annotate('text', x = 1, y = .25, label = '1.31-Fold **', size = 4) +
  annotate('text', x = 2, y = .25, label = '1.67-Fold **', size = 4) +
  annotate('text', x = 3, y = .6, label = '1.15-Fold **', size = 4) +
  annotate('text', x = 4, y = .9, label = '0.97-Fold **', size = 4) +
  annotate('text', x = 5, y = .3, label = '1.41-Fold ns', size = 4) +
  theme(aspect.ratio = 3/5)

#### 各个细胞中的差异分析和富集分析 ####
library(org.Mm.eg.db)

seurat <- readRDS('seurat.out.rds')

cell_name <- levels(seurat@active.ident)
cell_col <- c(Stem_cells = '#8dd3c7', Paneth_cells = '#ffed6f', Goblet_cells = '#bebada',
              Enteroendocrine_cells = '#fb8072', Tuft_cells = '#80b1d3',
              Enterocyte = '#fdb462', T_cells = '#b3de69', B_cells = '#bc80bd')

gene_info <- select(org.Mm.eg.db, keys(org.Mm.eg.db), 
                    keytype = "ENTREZID", 
                    columns = "SYMBOL")

# difference
diffs <- map(cell_name, \(x) {
  seurat_x <- subset(seurat, idents = x)
  FindMarkers(seurat_x, ident.1 = 'IPA_1', ident.2 = 'Veh_1', 
              group.by = 'orig.ident') %>% 
    rownames_to_column('name') %>% 
    mutate(enriched = ifelse(avg_log2FC > .1 & p_val < .05, 'up', 
                             ifelse(avg_log2FC < -.1 & p_val < .05, 
                                    'down', 'none'))) }) %>% 
  set_names(cell_name)
write.xlsx(diffs, 'cells.diff.xlsx')
saveRDS(diffs, 'cells.diff.rds')

# enrichment
eKEGGs <- map(diffs, \(x)
              filter(x, p_val < .05 & abs(avg_log2FC) > .25) %>% 
                left_join(gene_info, by = c('name' = 'SYMBOL')) %>% 
                pull(ENTREZID) %>% 
                enrichKEGG(organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1) %>% 
                data.frame() )
write.xlsx(eKEGGs, 'cells.diff.Gene.eKEGG.xlsx')

eKEGGs <- map(diffs, \(x)
              filter(x, enriched == 'up') %>% 
                left_join(gene_info, by = c('name' = 'SYMBOL')) %>% 
                pull(ENTREZID) %>% 
                enrichKEGG(organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1) %>% 
                data.frame() )
write.xlsx(eKEGGs, 'cells.diff.upGene.eKEGG.xlsx')

eKEGGs <- map(diffs, \(x)
              filter(x, enriched == 'down') %>% 
                left_join(gene_info, by = c('name' = 'SYMBOL')) %>% 
                pull(ENTREZID) %>% 
                enrichKEGG(organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1) %>% 
                data.frame() )
write.xlsx(eKEGGs, 'cells.diff.dnGene.eKEGG.xlsx')

eGOs <- map(diffs, \(x)
            filter(x, p_val < .05 & abs(avg_log2FC) > .25) %>% 
              pull(name) %>% 
              enrichGO(org.Mm.eg.db, keyType = 'SYMBOL', pvalueCutoff = 1, 
                       qvalueCutoff = 1, ont = 'ALL') ) 
write.xlsx(eGOs, 'cells.diff.Gene.eGO.xlsx')

eGOs <- map(diffs, \(x)
            filter(x, enriched == 'up') %>% 
              pull(name) %>% 
              enrichGO(org.Mm.eg.db, keyType = 'SYMBOL', pvalueCutoff = 1, 
                       qvalueCutoff = 1, ont = 'ALL') ) 
write.xlsx(eGOs, 'cells.diff.upGene.eGO.xlsx')

eGOs <- map(diffs, \(x)
            filter(x, enriched == 'down') %>% 
              pull(name) %>% 
              enrichGO(org.Mm.eg.db, keyType = 'SYMBOL', pvalueCutoff = 1, 
                       qvalueCutoff = 1, ont = 'ALL') ) 
write.xlsx(eGOs, 'cells.diff.dnGene.eGO.xlsx')

# visualization
cell_name <- c("Goblet_cells","Paneth_cells","Enterocyte","Stem_cells",
               "Tuft_cells","Enteroendocrine_cells")

terms <- c('GO:0042180','GO:0006631','GO:0042181','GO:0010565','GO:0050678',
           'GO:0007219','GO:0019827','GO:0045927')

eGOs_data <- map_dfr(cell_name, ~ 
                       read.xlsx('cells.diff.Gene.eGO.xlsx', sheet = .x) %>% 
                       filter(ID %in% terms) %>% 
                       dplyr::select(Description, FoldEnrichment, pvalue) %>%
                       add_column(cell = .x))

terms <- c('mmu04110','mmu03030','mmu00650','mmu04152','mmu03030')

eKEGGs_data <- map_dfr(cell_name, ~ 
                         read.xlsx('cells.diff.Gene.eKEGG.xlsx', sheet = .x) %>% 
                         filter(ID %in% terms) %>% 
                         dplyr::select(Description, FoldEnrichment, pvalue) %>%
                         add_column(cell = .x)) 

plot_data <- rbind(eGOs_data, eKEGGs_data) %>% 
  filter(pvalue < 0.05)

ggscatter(plot_data, 'cell', 'Description', fill = 'pvalue', shape = 21,
          size = 'FoldEnrichment', xlab = '', ylab = '', legend = 'right',
          x.text.angle = 30) +
  scale_fill_viridis_c(begin = .5) +
  scale_size_continuous(range = c(3, 6)) +
  theme(aspect.ratio = 2,
        panel.grid.major = element_line(linewidth = .5, color = 'grey88'))
