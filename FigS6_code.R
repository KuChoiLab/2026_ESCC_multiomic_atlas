# Fig S6A ----

p <- SpatialFeaturePlot(
  merged,
  features = c('CD8 Temra'), # GC B, Memory B, Naive B, CD4 Tpex, CD4 Tfh, fDC
  images = paste0("s", 1:10),
  ncol = 1, pt.size.factor = 2
) & theme(plot.title = element_blank(), plot.margin = margin(t = 37, r = 0, b = 37, l = 0)) & NoLegend() & theme(text = element_text(family = "Arial"), legend.text = element_text(size=6))#scale_fill_viridis_c(option = "magma")
cairo_pdf("figure/Fig6/20250823_escc_CD8Temra_2_default_dl.pdf", width = 5, height = 50, family="Arial")
print(p)
dev.off()


colors = scPalette(3)
names(colors) <- c("Dense lymphoid aggregates (r/o follicle)", "Germinal Center (MALT)", "Lymphoid follicle")


cairo_pdf("figure/FigS6/20250819_escc_TLS_pathological_annotation_magma_dl.pdf", width = 5, height = 50, family="Arial")
SpatialDimPlot(merged.gcb, group.by = 'pathologist_annotation', cols = colors, ncol=1, pt.size.factor = 2) & NoLegend() & theme(plot.title=element_blank(),plot.margin = margin(t = 37, r = 0, b = 37, l = 0))
dev.off()


# Fig S6B / Fig S6C ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(CellChat)
library(showtext)
library(MAST)
library(tidyverse)
library(ComplexHeatmap)
library(colorRamp2)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

date = "20250809"
project_name = "escc"

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig6/")

seurat_obj <- readRDS("~/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")


metadat <- seurat_obj@meta.data
metadat <- metadat %>% dplyr::mutate(category = case_when(pathologist_annotation %in% c("Dense lymphoid aggregates (r/o follicle)", "Germinal Center (MALT)", "Lymphoid follicle") ~ "TLS",
                                                          T ~ "Others"))

metadat <- metadat %>% dplyr::group_by(niche, category) %>% dplyr::summarise(n=n()) %>% dplyr::ungroup() %>% dplyr::group_by(niche) %>% dplyr::mutate(sum = sum(n)) %>% dplyr::ungroup() %>% dplyr::mutate(prop = n/sum)
metadat$x <- "1"

tmp <- metadat %>% dplyr::filter(category == "TLS") %>% dplyr::arrange(desc(prop)) %>% dplyr::pull(niche)
order <- c(as.character(tmp), levels(metadat$niche)[!levels(metadat$niche) %in% tmp])

metadat$niche <- factor(metadat$niche, levels = order)

pdf("20250713_escc_visium_TLS_pathologist_annotation_bayesspace_hj.pdf", width = 4, height = 4)
metadat %>% ggplot(aes(x = niche, y = prop, fill = category)) +
  geom_bar(stat = "identity", colour = "black") +
  ylab("Proportion") +
  xlab("BayesSpace") +
  theme_classic(base_family = "Arial") +
  theme(aspect.ratio = 1, axis.text = element_text(colour = "black"), axis.text.x = element_text(hjust = 1, angle = 45)) +
  scale_fill_manual(values = c("darkgrey", "gold"))
dev.off()


# https://www.science.org/doi/full/10.1126/sciadv.abg3750
tlsgenes <- c(
  "FDCSP", "CR2", "CXCL13", "LTF", "CD52", "MS4A1", "CCL19", "LINC00926", "LTB",
  "CORO1A", "CD79B", "TXNIP", "CD19", "LIMD2", "CD37", "ARHGAP45", "BLK", "TMC8",
  "CCL21", "PTPN6", "ATP2A3", "IGHM", "SPIB", "TMSB4X", "CXCR4", "NCF1", "CD79A",
  "ARHGAP9", "DEF6", "EVL", "TBC1D10C", "RASAL3", "INPP5D", "RNASET2", "RASGRP2",
  "TNFRSF13C", "RAC2", "CD22", "ARHGEF1", "AC103591.3", "TRAF3IP3", "HLA-DQB1",
  "CD53", "ARHGAP4", "TRBC2", "POU2AF1", "TRAF5", "OGA", "FCRL3", "HLA-DQA1"
)

DefaultAssay(seurat_obj) <- "Spatial"
seurat_obj <- AddModuleScore(object = seurat_obj, features = list(tlsgenes), name = "TLS")

chemokinegenes <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18",
                    "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13")

seurat_obj <- AddModuleScore(object = seurat_obj, features = list(chemokinegenes), name = "Chemokine")

tls_totalgenes <- unique(c(tlsgenes, chemokinegenes))

seurat_obj <- AddModuleScore(object = seurat_obj, features = list(tls_totalgenes), name = "TLS_Chemokine")

fileidentity <- "umap_res0.4_sub.cluster"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4*3, height = 3, family = "Arial")
FeaturePlot(seurat_obj, features = c("TLS1", "Chemokine1", "TLS_Chemokine1"), cols = c("grey", "red"), ncol = 3, reduction = "umap.harmony") & theme_classic(base_family = "Arial") & theme(aspect.ratio = 1)
dev.off()

fileidentity <- "TLS_vlnplot"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(5), height = cm_to_inch(4), family = "Arial")
df <- FetchData(seurat_obj, vars = c("TLS1", "niche"))
p <- ggplot(df, aes(x = niche, y = TLS1, fill = niche)) +
  geom_violin(linewidth = 0.2, color = "black", trim = TRUE) + 
  theme_classic(base_family = "Arial") +
  ylab("Module score") +
  xlab("Niche") +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 5),
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 5),
    line = element_line(linewidth = 0.2)
  ) +
  scale_fill_manual(values = scPalette(n=12)) +
  guides(fill = "none", color = "none")
print(p)
dev.off()


# Fig S6D ----

escc <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250731_escc_global_final_annotation_hj.rds")
escc@meta.data$pathology %>% table()

mycelltype <- setdiff(escc$Annotation_v2 %>% unique(), c("Mito-high CD4", ".", "Mito-high CD8"))

escc <- subset(escc, subset = Annotation_v2 %in% mycelltype)
escc@meta.data <- escc@meta.data %>% dplyr::mutate(Annotation_v2 = case_when(Annotation_v2 == "GC DZ" ~ "GC B",
                                                                             Annotation_v2 == "GC LZ" ~ "GC B",
                                                                             T ~ Annotation_v2
))
test <- escc@meta.data %>% dplyr::group_by(orig.ident, Annotation_v2) %>% dplyr::summarise(n=n()) %>% dplyr::group_by(orig.ident) %>% dplyr::mutate(sum = sum(n)) %>% dplyr::mutate(proportion = n/sum)

total_celltypes <- setdiff(test$Annotation_v2 %>% unique(), "GC B")

library(tidyr)

df <- data.frame()
plot_list <- list()
for (tmp in total_celltypes) {
  
  df_wide <- test %>% dplyr::filter(Annotation_v2 %in% c(tmp, "GC B")) %>%
    select(orig.ident, Annotation_v2, proportion) %>%
    pivot_wider(names_from = Annotation_v2, values_from = proportion, values_fill = 0)
  
  cor_res <- cor.test(x = df_wide[["GC B"]], y = df_wide[[tmp]])
  df_tmp <- data.frame(celltype = tmp, pval = cor_res$p.value, correlation = cor_res$estimate)
  df <- rbind(df, df_tmp)
  
  p <- df_wide %>% ggplot(aes(x = `GC B`, y = !!sym(tmp))) +
    geom_point() & theme(aspect.ratio = 1)
  plot_list[[tmp]] <- p
}

df$Padj <- p.adjust(df$pval, method = "BH")

library(ggrepel)
df <- df %>% dplyr::mutate(category = case_when(Padj < 0.05 & correlation > 0 ~ "positive correlation",
                                                Padj < 0.05 & correlation < 0 ~ "negative correlation",
                                                T ~ "others"))
anno <- df %>% dplyr::filter(!category %in% c("others"))
ymax <- max(-log10(df$pval), na.rm = TRUE)
ymin <- min(-log10(df$pval), na.rm = TRUE)

df %>% ggplot(aes(x = correlation, y = -log10(Padj), colour = category)) +
  geom_point(alpha = 0.5, size = 0.5, stroke = 0.1) +
  theme_classic() +
  xlab("Correlation") +
  xlim(-1, 1) +
  theme_classic(base_family = "Arial") +
  theme(aspect.ratio = 1, line = element_line(linewidth = 0.2), axis.text = element_text(size = 4, colour = "black"), axis.title = element_text(size = 5), legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.2, "cm"), legend.title = element_text(size = 4), legend.text = element_text(size = 3), axis.ticks.length = unit(0.05, "cm") ) +
  scale_color_manual(values = c("blue3", "black", "red3")) +
  guides(colour = guide_legend(override.aes = list(size = 1)))

# Fig S6E ----
library(Seurat)
library(ggplot2)

escc_subset <- readRDS("20250823_escc_gcb_viper_hj.rds")
escc_subset$pathology <- factor(escc_subset$pathology, levels = c("Normal", "Dysplasia", "Microinvasive", "Macroinvasive"))
VlnPlot(escc_subset, features = c("IRF8"), group.by = "pathology", assay = "viper") & geom_boxplot()
VlnPlot(escc_subset, features = c("IRF8"), group.by = "pathology") 

# Fig S6F ----
library(Seurat)
library(ggplot2)

seurat_obj <- readRDS("~/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")
metadat <- seurat_obj@meta.data
metadat <- metadat %>% dplyr::mutate(category = case_when(pathologist_annotation %in% c("Dense lymphoid aggregates (r/o follicle)", "Germinal Center (MALT)", "Lymphoid follicle") ~ "TLS",
                                                          T ~ "Others"))

metadat <- metadat %>% dplyr::group_by(niche, category) %>% dplyr::summarise(n=n()) %>% dplyr::ungroup() %>% dplyr::group_by(niche) %>% dplyr::mutate(sum = sum(n)) %>% dplyr::ungroup() %>% dplyr::mutate(prop = n/sum)
metadat$x <- "1"

tmp <- metadat %>% dplyr::filter(category == "TLS") %>% dplyr::arrange(desc(prop)) %>% dplyr::pull(niche)
order <- c(as.character(tmp), levels(metadat$niche)[!levels(metadat$niche) %in% tmp])

metadat$niche <- factor(metadat$niche, levels = order)

mydat <- seurat_obj@assays$C2L@data # col (barcode), row (celltype)

name_df <- data.frame(orig_name = mydat %>% rownames())
name_df <- name_df %>% dplyr::mutate(new_name = case_when(orig_name == "APOD.high.pFib" ~ "APOD high CAF",
                                                          orig_name == "C1QC..TAM" ~ "C1QC+ TAM",
                                                          orig_name == "CD1C..cDC2" ~ "CD1C+ cDC2",
                                                          orig_name == "CD4.Activated.Tcm" ~ "CD4 Activated Tcm",
                                                          orig_name == "CD4.CTL" ~ "CD4 CTL",
                                                          orig_name == "CD4.Early.Tcm" ~ "CD4 Early Tcm",
                                                          orig_name == "CD4.TNFRSF9..Treg" ~ "CD4 TNFRSF9- Treg",
                                                          orig_name == "CD4.TNFRSF9..Treg.1" ~ "CD4 TNFRSF9+ Treg",
                                                          orig_name == "CD4.Tcm" ~ "CD4 Tcm",
                                                          orig_name == "CD4.Terminal.Tex" ~ "CD4 Terminal Tex",
                                                          orig_name == "CD4.Tfh" ~ "CD4 Tfh",
                                                          orig_name == "CD4.Tn" ~ "CD4 Tn",
                                                          orig_name == "CD4.Tpex" ~ "CD4 Tpex",
                                                          orig_name == "CD4.differentiating.Treg" ~ "CD4 differentiating Treg",
                                                          orig_name == "CD8.CXCR6..Trm" ~ "CD8 CXCR6+ Trm",
                                                          orig_name == "CD8.GZMK..early.Tem" ~ "CD8 GZMK+ early Tem",
                                                          orig_name == "CD8.HSP..Trm" ~ "CD8 HSP+ Trm",
                                                          orig_name == "CD8.Inflammatory.Trm" ~ "CD8 Inflammatory Trm",
                                                          orig_name == "CD8.Temra" ~ "CD8 Temra",
                                                          orig_name == "CD8.Terminal.Tex" ~ "CD8 Terminal Tex",
                                                          orig_name == "CD8.Trm" ~ "CD8 Trm",
                                                          orig_name == "CLEC9A..cDC1" ~ "CLEC9A+ cDC1",
                                                          orig_name == "Classical.Monocyte" ~ "Classical Monocyte",
                                                          orig_name == "GC.DZ" ~ "GC DZ",
                                                          orig_name == "GC.LZ" ~ "GC LZ",
                                                          orig_name == "Glandular.cell" ~ "Glandular cell",
                                                          orig_name == "Interstitial.macrophage" ~ "Interstitial Macrophage",
                                                          orig_name == "LA.TAM" ~ "LA TAM",
                                                          orig_name == "LAMP3..cDC3" ~ "LAMP3+ cDC3",
                                                          orig_name == "Mast.cell" ~ "Mast cell",
                                                          orig_name == "Memory.B" ~ "Memory B",
                                                          orig_name == "Non.classical.Monocyte" ~ "Non classical Monocyte",
                                                          orig_name == "Plasma.cell" ~ "Plasma cell",
                                                          orig_name == "Proliferating.T" ~ "Proliferating T",
                                                          orig_name == "TIMP3.high.pFib" ~ "CFD high CAF",
                                                          orig_name == "naive.B" ~ "Naive B",
                                                          orig_name == "Proliferating.macrophage" ~ "Proliferating macrophage",
                                                          orig_name == "progenitor.iCAF" ~ "Progenitor iCAF",
                                                          T ~ orig_name))


celltypes <- rownames(mydat)

c2l_df <- data.frame()
for (celltype_tmp in celltypes) {
  mydat_tmp <- mydat[celltype_tmp, ]
  df_tmp <- data.frame(barcodes = names(mydat_tmp),
                       c2l = mydat_tmp)
  df_tmp$celltype <- celltype_tmp
  c2l_df <- rbind(c2l_df, df_tmp)
}

res_df <- data.frame()
my_TLS <- c("niche7")

metadat <- seurat_obj@meta.data

c2l_df <- merge(c2l_df, name_df, by.x = "celltype", by.y = "orig_name", all.x = T)

res_df <- data.frame()
for (celltype_tmp in unique(c2l_df$new_name)) {
  
  for (slide_tmp in c(paste("s", 2:10, sep= ""))) {
    
    key1 <- metadat %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(niche %in% my_TLS) %>% rownames()
    key2 <- metadat %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(!niche %in% my_TLS) %>% rownames()
    
    mycell1 <- c2l_df %>% dplyr::filter(new_name %in% celltype_tmp) %>% dplyr::filter(barcodes %in% key1)
    mycell2 <- c2l_df %>% dplyr::filter(new_name %in% celltype_tmp) %>% dplyr::filter(!barcodes %in% key1)
    
    res <- wilcox.test(mycell1$c2l, mycell2$c2l)
    pval <- res$p.value
    
    df_tmp <- data.frame(slide = slide_tmp, celltype = celltype_tmp, pvalue = pval, log2FC = log2(mean(mycell1$c2l)/mean(mycell2$c2l)))
    res_df <- rbind(res_df, df_tmp)
    
  }
  
}

res_df$Padj <- p.adjust(res_df$pvalue, method = "bonferroni")

head(res_df)

saveRDS(res_df, "20250819_escc_TLS_enrichment_niche7_hj.rds")

res_df <- readRDS("20250819_escc_TLS_enrichment_niche7_hj.rds")

res_df <- res_df %>% dplyr::filter(Padj < 0.05) %>% dplyr::filter(log2FC > 0) 

dys_test <- res_df %>% dplyr::filter(slide %in% c("s2", "s3")) %>% dplyr::group_by(celltype) %>% dplyr::summarise(n=n()) %>% dplyr::arrange(desc(n)) 
micro_test <- res_df %>% dplyr::filter(slide %in% c("s4", "s5", "s6", "s7", "s8", "s9")) %>% dplyr::group_by(celltype) %>% dplyr::summarise(n=n()) %>% dplyr::arrange(desc(n)) 
macro_test <- res_df %>% dplyr::filter(slide %in% c("s10")) %>% dplyr::group_by(celltype) %>% dplyr::summarise(n=n()) %>% dplyr::arrange(desc(n)) 

slide_levels <- paste0("s", 2:10)

df_wide <- res_df %>%
  select(slide, celltype, log2FC) %>%
  mutate(slide = factor(slide, levels = slide_levels),
         celltype = as.character(celltype)) %>%
  complete(celltype = unique(res_df$celltype),
           slide    = slide_levels,
           fill = list(log2FC = NA_real_)) %>%
  pivot_wider(names_from = slide, values_from = log2FC) %>%
  select(celltype, all_of(slide_levels))

mat <- df_wide %>%
  column_to_rownames("celltype") %>%
  as.matrix()

min_val <- min(as.matrix(mat)[as.matrix(mat) != 0], na.rm = TRUE)
rng <- range(as.matrix(mat)[as.matrix(mat) != 0], na.rm = TRUE)
mid_val <- mean(as.matrix(mat)[as.matrix(mat) != 0], na.rm = TRUE)
col_fun <- function(x) {
  ifelse(
    x == 0,
    "darkgrey",
    colorRamp2(c(rng[1], rng[2]), c("ivory","red3"))(x)
  )
}

set.seed(1234);ht <- Heatmap(
  mat,
  name = "log2FC",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = FALSE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 5, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "Arial"),
  width = unit(ncol(mat) * 0.2, "cm"),
  height = unit(nrow(mat) * 0.2, "cm"),
  row_dend_gp = gpar(lwd = 0.4),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x, y, w, h, gp = gpar(col = "black", fill = NA, lwd = 0.4))
  }
)

lg <- Legend(
  title = "log2FC",
  col_fun = colorRamp2(c(rng[1], rng[2]), c("ivory", "red3")),
  #at = pretty(c(min_val, max_val), 4),
  labels_gp = gpar(fontsize = 3, fontfamily = "Arial"),
  title_gp  = gpar(fontsize = 3, fontfamily = "Arial"),
  legend_height = unit(0.35, "cm"),
  grid_height   = unit(0.35, "cm"),
  grid_width    = unit(0.1, "cm"),
  border = "black"
)

#htd <- draw(ht, padding = unit(c(2, 30, 2, 2), "mm"))

library(vegan)
bin_mat <- (mat > 0) * 1
bin_mat[is.na(bin_mat)] <- 0
row_dend_bin <- hclust(vegan::vegdist(bin_mat, method = "jaccard", binary = TRUE), method = "centroid")
row_groups <- cutree(row_dend_bin, k = 13) 
set.seed(1234); ht <- Heatmap(
  t(mat),
  name = "log2FC",
  col = col_fun,
  cluster_columns = as.dendrogram(row_dend_bin),
  column_split    = 13,
  cluster_rows = F,
  #cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  show_heatmap_legend = FALSE,
  column_title = NULL,
  row_names_gp = gpar(fontsize = 5, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "Arial"),
  column_dend_gp = gpar(lwd = 0.3),
  row_gap = unit(1, "mm"),
  width  = unit(ncol(t(mat)) * 0.25, "cm"),
  height = unit(nrow(t(mat)) * 0.2, "cm"),
  cell_fun = function(j,i,x,y,w,h,fill) grid.rect(x,y,w,h, gp=gpar(col="black", fill=NA, lwd=0.4))
)
draw(ht, heatmap_legend_side = "right")
dev.off()
fileidentity <- "TLS_enrichment_log2FC_celltype_niche7_heatmap_split"
cairo_pdf(paste0("20250818", "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(15) , height = cm_to_inch(8), family = "Arial")
draw(
  ht,
  heatmap_legend_list = list(lg),  
  heatmap_legend_side = "right",
  padding = unit(c(2, 2, 2, 2), "mm")
)
dev.off()


# Fig S6G, H ----
library(CellChat)

cellchat <- readRDS("results/cellchat/20250819_escc_visium_TLS_niche7_tumor_result_dl.rds")

# Fig g — interaction strength heatmap
cairo_pdf("figure/Fig6/escc_TLS_niche7_tumor_cellchat_interaction_strength_heatmap.pdf",
          width = cm_to_inch(15), height = cm_to_inch(14))
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
dev.off()

# Fig h — fDC / FRC -> CD8 subsets bubble plot
cairo_pdf("figure/Fig6/escc_TLS_niche7_tumor_fDC_FRC_CD8_bubbleplot.pdf",
          width = cm_to_inch(17), height = cm_to_inch(8))
p <- netVisual_bubble(
  cellchat,
  sources.use = c("fDC", "FRC"),
  targets.use = c("CD8 GZMK+ early Tem", "CD8 HSP+ Trm", "CD8 Temra", "CD8 Trm"),
  remove.isolate = TRUE
)
print(p + coord_flip())
dev.off()

# Fig S6I ----

library(Seurat)
library(UCell)
library(msigdbr)
library(Hmisc)
library(reshape2)
library(tidyverse)

# subset CD8 Temra
cd8.temra <- subset(wh.escc.sc, Annotation_v2 == "CD8 Temra")
DefaultAssay(cd8.temra) <- if ("SCT" %in% names(cd8.temra@assays)) "SCT" else "RNA"

# GOBP gene sets (bias 제거: receptor 유전자 제외)
cd8_receptors <- c("CD8A", "CD8B", "CD99")

m_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  select(gs_name, gene_symbol)

genes_in_data <- rownames(cd8.temra)
gs_list <- split(m_df$gene_symbol, m_df$gs_name) %>%
  lapply(function(v) setdiff(intersect(unique(v), genes_in_data), toupper(cd8_receptors))) %>%
  .[sapply(., length) >= 10 & sapply(., length) <= 300]

names(gs_list) <- str_replace_all(names(gs_list), "[^A-Za-z0-9_]", "_")

# UCell score
expr <- GetAssayData(cd8.temra, slot = "data")
scores <- ScoreSignatures_UCell(expr, features = gs_list) %>%
  as.data.frame()
cd8.temra <- AddMetaData(cd8.temra, metadata = scores[Cells(cd8.temra), , drop = FALSE])

# scatter plot function
plot_scatter <- function(obj, receptor, term_col) {
  cells <- colnames(obj)
  df <- tibble(
    expr  = as.numeric(GetAssayData(obj, slot = "data")[receptor, cells]),
    score = as.numeric(obj@meta.data[cells, term_col])
  ) %>% filter(is.finite(expr), is.finite(score))
  
  st  <- cor.test(df$expr, df$score, method = "pearson")
  lab <- sprintf("r = %.3f\np %s", st$estimate,
                 ifelse(st$p.value < 2e-16, "< 2e-16", format.pval(st$p.value, digits = 2)))
  
  ggplot(df, aes(x = expr, y = score)) +
    geom_point(alpha = 0.35, size = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    annotate("label", x = -Inf, y = Inf, label = lab,
             hjust = -0.1, vjust = 1.1, size = 3, label.size = 0.3) +
    labs(x = paste0(receptor, " expression"), y = term_col) +
    theme_classic(base_family = "Arial", base_size = 8) +
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7))
}

p1 <- plot_scatter(cd8.temra, "CD8A", "GOBP_T_CELL_MEDIATED_CYTOTOXICITY")
p2 <- plot_scatter(cd8.temra, "CD8A", "GOBP_CELL_KILLING")
p3 <- plot_scatter(cd8.temra, "CD99", "GOBP_CYTOLYSIS")

cairo_pdf("figure/FigS6/escc_cd8Temra_receptor_GOBP_scatter.pdf",
          width = cm_to_inch(18), height = cm_to_inch(6), family = "Arial")
p1 | p2 | p3
dev.off()
