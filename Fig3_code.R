# Fig 3A ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(CellChat)
library(showtext)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/step02_subclustering/TIC/")
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

date = "20250810"
project_name = "escc"

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250806_escc_squamouscells_umap_tuning_sk_hj.rds")

fileidentity <- "Fig3A"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4*4, height = 4, family = "Arial")
DimPlot(seurat_obj, group.by = "Annotation_v2", split.by = "pathology") & theme_void(base_family = "Arial") & theme(aspect.ratio = 1) 
dev.off()

# Fig 3B ----

library(Seurat)
library(Matrix)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/TIC/")

seurat_obj <- readRDS("~/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")
DefaultAssay(seurat_obj) <- "RNA"

# save metadata table:
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='20251210_escc_squamous_final_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot = "counts")
writeMM(counts_matrix, file='20251210_escc_squamous_final_metadata.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$pca@cell.embeddings, file='20251210_escc_squamous_final_harmony.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)), file='20251210_escc_squamous_final_gene_names.csv',
  quote=F,row.names=F,col.names=F
)

# paga analysis in python

# Fig 3C ----

# scVelo in python

# Fig 3D ----

# CNV cutoff:
# expr < 0.95 = loss
# expr > 1.05 = gain

myinfer <- readRDS("/Volumes/hjdrive4/escc_inferCNV/v5/infercnvR_Epi_Mono_filtered/run.final.infercnv_obj")
data.frame(exp = as.vector(myinfer@expr.data)) %>% ggplot(aes(x = exp)) +
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = c(0.95, 1.05)) +
  theme_bw()

mymeta <- seurat_obj@meta.data %>% dplyr::select(barcode, Annotation_v3)
myinfer@expr.data[(myinfer@expr.data < 1.05 & myinfer@expr.data > 0.95)] <- NA
myinfer@expr.data[(myinfer@expr.data > 1.05 | myinfer@expr.data < 0.95)] <- 1
#sum(myinfer@expr.data, na.rm = T)
myinfer@expr.data[is.na(myinfer@expr.data)] <- 0
mycount <- colSums(myinfer@expr.data)
mycount_df <- data.frame(mycount)
mycount_df$barcode <- names(mycount)
mycount_df <- merge(mycount_df, mymeta, all.x = T, by.x = "barcode", by.y = "barcode")
mycount_df <- mycount_df %>% dplyr::filter(!is.na(Annotation_v3))

pdf("20251215_escc_inferCNV_malignant_score_frequency_hj.pdf", width = 5, height = 5)
mycount_df %>% dplyr::filter(Annotation_v3 != "LGR5+ cancer stem cell") %>% ggplot(aes(x = Annotation_v3, y = mycount, fill = Annotation_v3)) +
  geom_boxplot(outliers = F) +
  ylab("Frequency of CNV") +
  xlab("Cell types") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  stat_summary(fun = mean, geom = "point", size = 1, color = "white") +       # mean point
  theme_classic() +
  scale_fill_manual(values = c("Basal" = "#F8766D", "Suprabasal" = "#0CB702", "Transitional suprabasal" = "blue4", "Tumor initiating cell" = "magenta2", Tumor = "#619CFF")) +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 20, hjust = 1), text = element_text(family = "Arial"))
dev.off()

# stat test

mycomb <- combn(unique(mycount_df$Annotation_v3), m = 2)
dim(mycomb)

wil_df <- data.frame()

for (i in 1:dim(mycomb)[2]) {
  celltype1 <- mycomb[1,i]
  celltype2 <- mycomb[2,i]
  x <- mycount_df %>% dplyr::filter(Annotation_v3 == celltype1) %>% dplyr::pull(mycount)
  y <- mycount_df %>% dplyr::filter(Annotation_v3 == celltype2) %>% dplyr::pull(mycount)
  res <- wilcox.test(x, y)
  pvalue <- res$p.value
  df_tmp <- data.frame(group1 = celltype1, group2 = celltype2, pvalue = pvalue, median1 = median(x), median2 = median(y))
  wil_df <- rbind(df_tmp, wil_df)
}

wil_df$Padj <- p.adjust(wil_df$pvalue, method = "BH")

saveRDS(wil_df, "20251219_escc_TIC_CNV_frequency_wilcoxon_hj.rds")

ttest_df <- data.frame()

for (i in 1:dim(mycomb)[2]) {
  celltype1 <- mycomb[1,i]
  celltype2 <- mycomb[2,i]
  x <- mycount_df %>% dplyr::filter(Annotation_v3 == celltype1) %>% dplyr::pull(mycount)
  y <- mycount_df %>% dplyr::filter(Annotation_v3 == celltype2) %>% dplyr::pull(mycount)
  res <- t.test(x, y)
  pvalue <- res$p.value
  df_tmp <- data.frame(group1 = celltype1, group2 = celltype2, pvalue = pvalue, median1 = median(x), median2 = median(y))
  ttest_df <- rbind(df_tmp, ttest_df)
}

ttest_df$Padj <- p.adjust(ttest_df$pvalue, method = "BH")

saveRDS(ttest_df, "20251219_escc_TIC_CNV_frequency_ttest_hj.rds")

test <- subset(seurat_obj, subset = Annotation_v3 == "LGR5+ cancer stem cell")
test <- subset(seurat_obj, subset = Annotation_v3 %in% c("Tumor initiating cell", "LGR5+ cancer stem cell"))
VlnPlot(test, features = "POLQ", group.by = "orig.ident")
VlnPlot(test, features = "", group.by = "orig.ident", split.by = "Annotation_v3")
test$orig.ident

# Fig 3E ----
seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")

old <- GeomViolin$default_aes$linewidth
update_geom_defaults("violin", list(linewidth = 0.2))
fileidentity <- "TIC_stacked_violinplot"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(6), height = cm_to_inch(4.5), family="Arial")
VlnPlot(seurat_obj, features = c("C2orf48", "CHTF18", "CREB5", "PCGF2", "AKR1C2", "PSMB3", "AKR1C1", "PIP4K2B", "STARD3"), group.by = "Annotation_v3", split.by = "Annotation_v3", stack = T, flip = F) +
  theme_classic(base_family = "Arial") +
  ylab("Cell type") +
  theme(
    line = element_line(linewidth = 0.2),
    axis.text = element_text(colour = "black", size = 4),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    #axis.ticks = element_blank(),
    axis.ticks.length = unit(0.05, "cm"),
    
    legend.position = "None",
    legend.text = element_text(size = 3),
    strip.background = element_blank(),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"),
    panel.spacing.x = unit(0, "cm"),
    strip.text.x = element_text(angle = 90, hjust = 0, size = 4),
    panel.grid = element_blank()) &
  scale_fill_manual(values = c("#F8766D", "#0CB702", "blue4", "magenta2", "gold2", "#619CFF"))
dev.off()

# Fig 3F ----
library(Seurat)
library(dplyr)
library(AUCell)
library(SCENIC)
library(ComplexHeatmap)
library(stringr)
library(showtext)
library(ComplexHeatmap)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

setwd(paste0("/Users/hojin/Dropbox/project/ESCC/submit/analysis/step02_subclustering/TIC/"))

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250806_escc_squamouscells_umap_tuning_sk_hj.rds")
metadat <- seurat_obj@meta.data

#seurat_obj <- subset(seurat_obj, subset = Annotation_v3 == "Tumor initiating cell")

seurat_obj@meta.data <- seurat_obj@meta.data %>%
  dplyr::mutate(pathology = case_when(pathology == "Microinvasive" ~ "Microinvasive carcinoma",
                                      pathology == "Macroinvasive" ~ "Macroinvasive carcinoma",
                                      T ~ pathology
  ))

seurat_obj$pathology <- factor(seurat_obj$pathology, levels = c("Normal", "Dysplasia", "Microinvasive carcinoma", "Macroinvasive carcinoma"))


library(ggalluvial)
#seurat_obj$Annotation_v3 <- as.vector(seurat_obj$Annotation_v3)
table_ESCC <- seurat_obj@meta.data %>% dplyr::select(pathology, Annotation_v3) %>% table() 
table_ESCC2 <- ( table_ESCC / table_ESCC %>% apply(1,sum) ) %>% data.frame()
table_ESCC2$pathology <- factor(table_ESCC2$pathology,  levels = c("Normal", "Dysplasia", "Microinvasive carcinoma", "Macroinvasive carcinoma"))

size_pt <- 7

mycol <- c("Basal" = "#F8766D", "Suprabasal" = "#0CB702", "Transitional suprabasal" = "blue4", "Tumor initiating cell" = "magenta2", "LGR5+ cancer stem cell" = "gold2", "Tumor" = "#619CFF")
#mycol <- c("#E41A1C", "#F781BF", "darkgreen", "#4DAF4A", "green", "#222F75","blue", "#377EB8", "#54B0E4", "cyan")
#names(mycol) <- c("Squamous cell", "Glandular cell", "Endothelial cell", "Fibroblast/fDC/FRC", "SMC/Pericyte", "B cell", "Plasma cell", "T/NK cell", "Monocyte/DC/Macrophage", "Mast")
date <- "20250812"
project_name <- "escc"

fileidentity <- "Fig3_tic_prop"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(5), height = cm_to_inch(5), family = "Arial")
ggplot(table_ESCC2, aes( x = pathology, y = Freq, alluvium = Annotation_v3)) + 
  geom_alluvium(aes(fill = Annotation_v3), colour = "black", alpha = 0.9, decreasing = FALSE, lwd = 0.3) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) + scale_x_discrete(expand = c(0,0)) + 
  theme_bw(base_family = "Arial") +
  theme(plot.title = element_blank(), legend.position = "none", 
        line = element_line(linewidth = 0.3),
        legend.text = element_text(size = size_pt-3), legend.title = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 15, size = size_pt-1, vjust = 0.9, hjust = 0.95), 
        axis.text.y = element_text(colour = "black", size = size_pt-1),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank(), panel.grid = element_blank()) +
  scale_fill_manual(values = mycol)
dev.off()

# Fig 3G ----

library(RColorBrewer)
moduleA  <- "Module.2"        
moduleB  <- "Module.0"        
assay  <- DefaultAssay(obj)   
slot   <- "scale.data"              

mat <- GetAssayData(obj, assay = assay, slot = slot)[c(moduleA, moduleB), , drop = FALSE]

diff_vec <- as.numeric(mat[moduleA, ] - mat[moduleB, ])
names(diff_vec) <- colnames(obj)

colname <- paste0("diff_", moduleA, "_minus_", moduleB)
obj[[colname]] <- diff_vec

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

p <- SpatialFeaturePlot(
  obj,
  features = colname,
  images = paste0("s", 1:10),
  ncol = 1, pt.size.factor = 1.8
) &
  theme(plot.title = element_blank(),
        plot.margin = margin(t = 37, r = 0, b = 37, l = 0)) &
  NoLegend() &
  theme(text = element_text(family = "Arial"),
        legend.text = element_text(size = 6)) &
  scale_fill_gradientn(colours=SpatialColors(n=100), na.value = 'black')

cairo_pdf("figure/Fig3/20250824_escc_Mod2_minus_Mod0_default_test_dl.pdf", width = 5, height = 50, family="Arial")
print(p)
dev.off()

# Fig 3I -----

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/NEBULA/")

seurat_obj <- readRDS("~/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig4/supple_scenic/20250806_escc_squamouscells_umap_tuning_sk_hj.rds")

set.seed(1234); sample1 <- seurat_obj@meta.data %>% dplyr::filter(Annotation_v3 == "Basal") %>% rownames() %>%
  sample(size = nrow(seurat_obj@meta.data %>% dplyr::filter(Annotation_v3 == "Basal"))*0.1, replace = F)
set.seed(1234); sample2 <- seurat_obj@meta.data %>% dplyr::filter(Annotation_v3 == "Suprabasal") %>% rownames() %>%
  sample(size = nrow(seurat_obj@meta.data %>% dplyr::filter(Annotation_v3 == "Suprabasal"))*0.1, replace = F)
set.seed(1234); sample3 <- seurat_obj@meta.data %>% dplyr::filter(Annotation_v3 == "Transitional suprabasal") %>% rownames()
set.seed(1234); sample4 <- seurat_obj@meta.data %>% dplyr::filter(Annotation_v3 == "Tumor initiating cell") %>% rownames()
set.seed(1234); sample5 <- seurat_obj@meta.data %>% dplyr::filter(Annotation_v3 == "LGR5+ cancer stem cell") %>% rownames()
set.seed(1234); sample6 <- seurat_obj@meta.data %>% dplyr::filter(sub.cluster == "4") %>% rownames() %>%
  sample(size = nrow(seurat_obj@meta.data %>% dplyr::filter(sub.cluster == "4"))*0.1, replace = F)

mycells <- c(sample1, sample2, sample3, sample4, sample5, sample6)

seurat_obj_subset <- subset(seurat_obj, cells = mycells)

# TIC vs transitional suprabasal
seurat_obj_subset_v2 <- subset(seurat_obj_subset, subset = Annotation_v3 %in% c("Tumor initiating cell", "Transitional suprabasal"))

seurat_obj_subset_v2@meta.data <- seurat_obj_subset_v2@meta.data %>% dplyr::mutate(nebula_anno = case_when(Annotation_v3 == "Tumor initiating cell" ~ "Tumor initiating cell",
                                                                                                           T ~ "Transitional suprabasal"))

seuratdata <- scToNeb(obj = seurat_obj_subset_v2, assay = "RNA", id = "donor", pred = c("nebula_anno", "donor"), offset = "nCount_RNA")
df = model.matrix(~nebula_anno+donor, data = seuratdata$pred)
## include only the first two cell types in the model to avoid separation due to too many binary variables
data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoTumor initiating cell")],offset=seuratdata$offset)
re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset)

saveRDS(re, "nebula_TIC_vs_transitional_suprabasal_hj.rds")

seurat_obj_subset@meta.data <- seurat_obj_subset@meta.data %>% dplyr::mutate(nebula_anno = case_when(Annotation_v2 == "Suprabasal" ~ "Suprabasal",
                                                                                                     T ~ "Others"))

# TIC vs Tumor
seurat_obj_subset_v2 <- subset(seurat_obj_subset, subset = sub.cluster %in% c("4", "7_0", "7_1", "7_2", "7_4"))

seurat_obj_subset_v2@meta.data <- seurat_obj_subset_v2@meta.data %>% dplyr::mutate(nebula_anno = case_when(Annotation_v3 == "Tumor initiating cell" ~ "Tumor initiating cell",
                                                                                                           T ~ "Tumor"))

seuratdata <- scToNeb(obj = seurat_obj_subset_v2, assay = "RNA", id = "donor", pred = c("nebula_anno", "donor"), offset = "nCount_RNA")
df = model.matrix(~nebula_anno+donor, data = seuratdata$pred)
## include only the first two cell types in the model to avoid separation due to too many binary variables
data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoTumor initiating cell")],offset=seuratdata$offset)
re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset)

saveRDS(re, "nebula_TIC_vs_tumor_hj.rds")

# GSEA
library(fgsea)
library(tibble)

celltype <- "Transitional suprabasal"
re <- readRDS("nebula_TIC_vs_transitional_suprabasal_hj.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoTumor initiating cell`, `se_nebula_annoTumor initiating cell`, `p_nebula_annoTumor initiating cell`)
res$Padj <- p.adjust(res$`p_nebula_annoTumor initiating cell`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "Padj")

df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)

df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("/Users/hojin/Dropbox/db/h.all.v2025.1.Hs.symbols.gmt")
set.seed(1234); fgseaRes <- fgseaMultilevel(pathways = pw, stats = ranks, minSize = 15, maxSize = 500)
fgseaRes <- as.data.frame(fgseaRes)
fgseaRes_final_df <- data.frame()
for (i in seq(1:nrow(fgseaRes))) {
  fgseaRes_final_df_tmp <- fgseaRes[i,]
  leading_edges_tmp <- fgseaRes[i,]$leadingEdge %>% unlist()
  fgseaRes_final_df_tmp$leadingEdge <- paste(leading_edges_tmp, collapse = "/")
  
  fgseaRes_final_df <- rbind(fgseaRes_final_df, fgseaRes_final_df_tmp)
}

fgseaRes_final_df <- fgseaRes_final_df %>% dplyr::filter(padj < 0.05)
test1 <- fgseaRes_final_df

celltype <- "Tumor"
re <- readRDS("nebula_TIC_vs_tumor_hj.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoTumor initiating cell`, `se_nebula_annoTumor initiating cell`, `p_nebula_annoTumor initiating cell`)
res$Padj <- p.adjust(res$`p_nebula_annoTumor initiating cell`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "Padj")

df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)

df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("/Users/hojin/Dropbox/db/h.all.v2025.1.Hs.symbols.gmt")
set.seed(1234); fgseaRes <- fgseaMultilevel(pathways = pw, stats = ranks, minSize = 15, maxSize = 500)
fgseaRes <- as.data.frame(fgseaRes)
fgseaRes_final_df <- data.frame()
for (i in seq(1:nrow(fgseaRes))) {
  fgseaRes_final_df_tmp <- fgseaRes[i,]
  leading_edges_tmp <- fgseaRes[i,]$leadingEdge %>% unlist()
  fgseaRes_final_df_tmp$leadingEdge <- paste(leading_edges_tmp, collapse = "/")
  
  fgseaRes_final_df <- rbind(fgseaRes_final_df, fgseaRes_final_df_tmp)
}

fgseaRes_final_df <- fgseaRes_final_df %>% dplyr::filter(padj < 0.05)
test2 <- fgseaRes_final_df

mypath <- c(test1$pathway, test2$pathway)

test1$celltype <- "TIC vs Transitional suprabasal"
test2$celltype <- "TIC vs Tumor"

test <- rbind(test1, test2)

max_size_val <- max(-log10(test$padj), na.rm = TRUE)
min_size_val <- min(-log10(test$padj), na.rm = TRUE)

test$pathway <- gsub(test$pathway, pattern = "HALLMARK_", replacement = "")
test$pathway <- gsub(test$pathway, pattern = "_", replacement = " ")
test$pathway <- paste0(
  toupper(substr(test$pathway, 1, 1)),
  tolower(substr(test$pathway, 2, nchar(test$pathway)))
)

library(ggplot2)
library(dplyr)

test <- test %>% dplyr::mutate(bar_color = case_when(NES > 0 ~ "TIC",
                                                     NES < 0 ~ "Others"))

test %>%
  ggplot(aes(y = pathway, x = NES, fill = bar_color
  )) +
  geom_col(width = 0.7) +
  facet_grid(. ~ celltype, scales = "free_x", space = "free_x") +
  scale_fill_manual(values= c("blue3","red3")) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text.x = element_text(size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    axis.title.y = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) + scale_x_continuous(
    limits = c(-3, 3),
    breaks = c(-3, -1.5, 0, 1.5, 3)
  )

# Fig 3J -----

library(Seurat)
library(ggplot2)
library(SeuratObject)
library(patchwork)
library(fgsea)
library(dplyr)
library(fgsea)
library(showtext)
library(reshape2)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig3/monocyte_enrichment")

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

date = "20250820"
project_name <- "escc"

old <- GeomViolin$default_aes$linewidth
update_geom_defaults("boxplot", list(linewidth = 0.2))

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")

# distribution ----

mydist <- list()

seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(niche = case_when(niche %in% c("niche8") ~ "niche8",
                                                                                 T ~ "others"))
seurat_obj$classical_monocyte <- seurat_obj@assays$C2L@data["Classical.Monocyte",]
seurat_obj$non_classical_monocyte <- seurat_obj@assays$C2L@data["Non.classical.Monocyte",]
meta_tmp <- seurat_obj@meta.data

p <- meta_tmp %>% ggplot(aes(x = classical_monocyte)) +
  geom_histogram(binwidth = 0.000001) 
#geom_vline(xintercept = 0.000005)

p2 <- meta_tmp %>% ggplot(aes(x = classical_monocyte)) +
  xlim(-0.00001, 0.0001) +
  geom_histogram(binwidth = 0.000001) +
  geom_vline(xintercept = 0.000005)

p3 <- meta_tmp %>% ggplot(aes(x = non_classical_monocyte)) +
  geom_histogram(binwidth = 0.000001) 
#geom_vline(xintercept = 0.000005)

p4 <- meta_tmp %>% ggplot(aes(x = non_classical_monocyte)) +
  xlim(-0.00001, 0.0001) +
  geom_histogram(binwidth = 0.000001) +
  geom_vline(xintercept = 0.000005)

library(patchwork)
pdf("histogram_hj.pdf", width = 30, height = 3)
p+p2+p3+p4+plot_layout(ncol = 4)
dev.off()

# classical vs non classical in niche8 ----

totalres <- list()
wilres <- list()
for (slide_tmp in paste0("s", 1:10)) {
  seurat_obj_tmp <- subset(seurat_obj, subset = orig.ident == slide_tmp)
  seurat_obj_tmp@meta.data <- seurat_obj_tmp@meta.data %>% dplyr::mutate(niche = case_when(niche %in% c("niche8") ~ "niche8",
                                                                                           T ~ "others"))
  seurat_obj_tmp$classical_monocyte <- seurat_obj_tmp@assays$C2L@data["Classical.Monocyte",]
  seurat_obj_tmp$non_classical_monocyte <- seurat_obj_tmp@assays$C2L@data["Non.classical.Monocyte",]
  
  metadat <- seurat_obj_tmp@meta.data
  metadat <- metadat %>% dplyr::select(niche, classical_monocyte, non_classical_monocyte) 
  metadat <- melt(metadat, id.vars = "niche")
  metadat$slide <- slide_tmp
  
  # test
  x <- metadat %>% dplyr::filter(niche == "niche8") %>% dplyr::filter(variable == "non_classical_monocyte") %>% dplyr::filter(value > 0.000005) %>% dplyr::pull(value)
  y <- metadat %>% dplyr::filter(niche == "niche8") %>% dplyr::filter(variable == "classical_monocyte") %>% dplyr::filter(value > 0.000005) %>% dplyr::pull(value)
  
  wilres_tmp <- wilcox.test(x, y)
  pval_tmp <- wilres_tmp$p.value
  
  myvalue <- metadat %>% dplyr::filter(value > 0.000005) %>% dplyr::group_by(niche, variable) %>% dplyr::summarise(mean_prop = mean(value)) 
  group1 <- myvalue %>% dplyr::filter(niche == "niche8") %>% dplyr::filter(variable == "non_classical_monocyte") %>% pull(mean_prop)
  group2 <- myvalue %>% dplyr::filter(niche == "niche8") %>% dplyr::filter(variable == "classical_monocyte") %>% pull(mean_prop)
  
  df_tmp <- data.frame(slide = slide_tmp, pval = pval_tmp, non_classical_monocyte_mean = group1, classical_monocyte_mean = group2)
  wilres[[slide_tmp]] <- df_tmp
  totalres[[slide_tmp]] <- metadat
}

wilres <- do.call(rbind, wilres); totalres <- do.call(rbind, totalres)

wilres$padj <- p.adjust(wilres$pval, method = "BH")
wilres <- wilres %>% dplyr::mutate(log2FC = log2(non_classical_monocyte_mean/classical_monocyte_mean))
saveRDS(totalres, "20250823_escc_classical_nonclassical_niche8_compare_total_hj.rds")
saveRDS(wilres, "20250823_escc_classical_nonclassical_niche8_compare_hj.rds")

totalres <- readRDS("20250823_escc_classical_nonclassical_niche8_compare_total_hj.rds")
wilres <- readRDS("20250823_escc_classical_nonclassical_niche8_compare_hj.rds")

wilres %>% dplyr::filter(padj < 0.05) %>% dplyr::mutate(log2FC = log2(non_classical_monocyte_mean/classical_monocyte_mean))

totalres$slide <- factor(totalres$slide, levels = paste0("s", 1:10))
totalres %>% head()

library(scales)
totalres$variable <- factor(totalres$variable, levels = c("non_classical_monocyte", "classical_monocyte"))
p <- totalres %>% dplyr::filter(niche == "niche8") %>% ggplot(aes(x = slide, y = value, fill = variable)) +
  geom_boxplot(outliers = F) +
  theme_classic(base_family = "Arial") +
  ylab("log10(Fraction for each spot)") +
  xlab("Slide") +
  theme(aspect.ratio = 0.6, axis.text = element_text(size = 5, colour = "black"), line = element_line(linewidth = 0.2), axis.title = element_text(size = 5),
        legend.key.height = unit(0.2,  "cm"), legend.key.width = unit(0.2,  "cm"), legend.title = element_text(size = 3), legend.text = element_text(size = 3),
        axis.ticks.length = unit(0.05, "cm")) +
  scale_y_log10(labels = label_number()) +
  scale_fill_manual(values = rev(c("green3", "red3")))

library(patchwork)
project_name <- "escc"
fileidentity <- "monocyte_fraction_barplot_cutoff_0.000005_main"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(9), height = cm_to_inch(6), family = "Arial")
p
dev.off()


myres_mean <- list()
for (slide_tmp in paste0("s", 1:10)) {
  seurat_obj_tmp <- subset(seurat_obj, subset = orig.ident == slide_tmp)
  seurat_obj_tmp@meta.data <- seurat_obj_tmp@meta.data %>% dplyr::mutate(sq_niche = case_when(niche %in% c("niche8") ~ "sq_niche",
                                                                                              T ~ "others"))
  seurat_obj_tmp$classical_monocyte <- seurat_obj_tmp@assays$C2L@data["Classical.Monocyte",]
  seurat_obj_tmp$non_classical_monocyte <- seurat_obj_tmp@assays$C2L@data["Non.classical.Monocyte",]
  
  metadat <- seurat_obj_tmp@meta.data
  metadat <- metadat %>% dplyr::select(sq_niche, classical_monocyte, non_classical_monocyte) 
  metadat <- melt(metadat, id.vars = "sq_niche")
  myvalue <- metadat %>% dplyr::group_by(sq_niche, variable) %>% dplyr::summarise(mean = mean(value))
  #myvalue <- metadat %>% dplyr::group_by(sq_niche, variable) %>% dplyr::summarise(sum_prop = sum(value)) 
  myvalue$slide <- slide_tmp
  myres_mean[[slide_tmp]] <- myvalue
}

myres_mean <- do.call(rbind, myres_mean)

head(myres_mean)

myres_mean$slide <- factor(myres_mean$slide, levels = paste0("s", 1:10))
myres_mean %>% dplyr::filter(variable == "non_classical_monocyte") %>% ggplot(aes(x = slide, y = mean, fill = sq_niche)) +
  geom_bar(stat = "identity", position = "dodge")
myres_mean %>% dplyr::filter(variable == "classical_monocyte") %>% ggplot(aes(x = slide, y = mean, fill = sq_niche)) +
  geom_bar(stat = "identity", position = "dodge")


# Tumor vs other ----

totalres_nonclassical <- list()
wilres_nonclassical <- list()
for (slide_tmp in paste0("s", 1:10)) {
  seurat_obj_tmp <- subset(seurat_obj, subset = orig.ident == slide_tmp)
  seurat_obj_tmp@meta.data <- seurat_obj_tmp@meta.data %>% dplyr::mutate(niche = case_when(niche %in% c("niche8") ~ "niche8",
                                                                                           T ~ "others"))
  seurat_obj_tmp$classical_monocyte <- seurat_obj_tmp@assays$C2L@data["Classical.Monocyte",]
  seurat_obj_tmp$non_classical_monocyte <- seurat_obj_tmp@assays$C2L@data["Non.classical.Monocyte",]
  
  # non calssical
  metadat <- seurat_obj_tmp@meta.data
  metadat <- metadat %>% dplyr::select(niche, non_classical_monocyte) 
  metadat <- melt(metadat, id.vars = "niche")
  metadat$slide <- slide_tmp
  
  # test
  x <- metadat %>% dplyr::filter(niche == "niche8") %>% dplyr::filter(value > 0.000005) %>%  dplyr::pull(value)
  y <- metadat %>% dplyr::filter(niche == "others")%>% dplyr::filter(value > 0.000005) %>%  dplyr::pull(value)
  
  wilres_tmp <- wilcox.test(x, y)
  pval_tmp <- wilres_tmp$p.value
  
  myvalue <- metadat %>% dplyr::group_by(niche, variable) %>% dplyr::summarise(mean_prop = mean(value)) 
  group1 <- myvalue %>% dplyr::filter(niche == "niche8") %>% pull(mean_prop)
  group2 <- myvalue %>% dplyr::filter(niche == "others") %>% pull(mean_prop)
  
  df_tmp <- data.frame(slide = slide_tmp, pval = pval_tmp, niche8_mean = group1, others_mean = group2)
  
  wilres_nonclassical[[slide_tmp]] <- df_tmp
  totalres_nonclassical[[slide_tmp]] <- metadat
}

totalres_classical <- list()
wilres_classical <- list()
for (slide_tmp in paste0("s", 1:10)) {
  seurat_obj_tmp <- subset(seurat_obj, subset = orig.ident == slide_tmp)
  seurat_obj_tmp@meta.data <- seurat_obj_tmp@meta.data %>% dplyr::mutate(niche = case_when(niche %in% c("niche8") ~ "niche8",
                                                                                           T ~ "others"))
  seurat_obj_tmp$classical_monocyte <- seurat_obj_tmp@assays$C2L@data["Classical.Monocyte",]
  seurat_obj_tmp$non_classical_monocyte <- seurat_obj_tmp@assays$C2L@data["Non.classical.Monocyte",]
  
  # non calssical
  metadat <- seurat_obj_tmp@meta.data
  metadat <- metadat %>% dplyr::select(niche, classical_monocyte) 
  metadat <- melt(metadat, id.vars = "niche")
  metadat$slide <- slide_tmp
  
  # test
  x <- metadat %>% dplyr::filter(niche == "niche8") %>% dplyr::filter(value > 0.000005) %>%  dplyr::pull(value)
  y <- metadat %>% dplyr::filter(niche == "others")%>% dplyr::filter(value > 0.000005) %>%  dplyr::pull(value)
  
  wilres_tmp <- wilcox.test(x, y)
  pval_tmp <- wilres_tmp$p.value
  
  myvalue <- metadat %>% dplyr::group_by(niche, variable) %>% dplyr::summarise(mean_prop = mean(value)) 
  group1 <- myvalue %>% dplyr::filter(niche == "niche8") %>% pull(mean_prop)
  group2 <- myvalue %>% dplyr::filter(niche == "others") %>% pull(mean_prop)
  
  df_tmp <- data.frame(slide = slide_tmp, pval = pval_tmp, niche8_mean = group1, others_mean = group2)
  
  wilres_classical[[slide_tmp]] <- df_tmp
  totalres_classical[[slide_tmp]] <- metadat
}

wilres_nonclassical <- do.call(rbind, wilres_nonclassical)
totalres_nonclassical <- do.call(rbind, totalres_nonclassical)
wilres_nonclassical$padj <- p.adjust(wilres_nonclassical$pval, method = "BH")
wilres_nonclassical <- wilres_nonclassical %>% dplyr::mutate(log2FC = log2(niche8_mean/others_mean))
saveRDS(wilres_nonclassical, "20250823_escc_nonclasscial_test_hj.rds")
saveRDS(totalres_nonclassical, "20250823_escc_nonclasscial_total_hj.rds")

wilres_classical <- do.call(rbind, wilres_classical)
totalres_classical <- do.call(rbind, totalres_classical)
wilres_classical$padj <- p.adjust(wilres_classical$pval, method = "BH")
wilres_classical <- wilres_classical %>% dplyr::mutate(log2FC = log2(niche8_mean/others_mean))
saveRDS(wilres_classical, "20250823_escc_classcial_test_hj.rds")
saveRDS(totalres_classical, "20250823_escc_classcial_total_hj.rds")

wilres_nonclassical %>% dplyr::filter(padj < 0.05)
wilres_classical %>% dplyr::filter(padj < 0.05)

library(scales)

totalres_nonclassical$slide <- factor(totalres_nonclassical$slide, levels = paste0("s", 1:10))
p2 <- totalres_nonclassical %>% ggplot(aes(x = slide, y = value, fill = niche)) +
  geom_boxplot(outliers = F) +
  theme_classic(base_family = "Arial") +
  ylab("Fraction for each spot") +
  xlab("Slide") +
  theme(aspect.ratio = 0.5, axis.text = element_text(size = 5, colour = "black"), line = element_line(linewidth = 0.2), axis.title = element_text(size = 5),
        legend.key.height = unit(0.2,  "cm"), legend.key.width = unit(0.2,  "cm"), legend.title = element_text(size = 3), legend.text = element_text(size = 3)    ,
        axis.ticks.length = unit(0.05, "cm")) +
  scale_y_log10(labels = label_number()) +
  scale_fill_manual(values = c("red3", "blue3"))

wilres_nonclassical %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FC > 0)

totalres_classical$slide <- factor(totalres_classical$slide, levels = paste0("s", 1:10))
p3 <- totalres_classical %>% ggplot(aes(x = slide, y = value, fill = niche)) +
  geom_boxplot(outliers = F) +
  theme_classic(base_family = "Arial") +
  ylab("Fraction for each spot") +
  xlab("Slide") +
  theme(aspect.ratio = 0.5, axis.text = element_text(size = 5, colour = "black"), line = element_line(linewidth = 0.2), axis.title = element_text(size = 5),
        legend.key.height = unit(0.2,  "cm"), legend.key.width = unit(0.2,  "cm"), legend.title = element_text(size = 3), legend.text = element_text(size = 3),
        axis.ticks.length = unit(0.05, "cm")) +
  scale_y_log10(labels = label_number()) +
  scale_fill_manual(values = c("red3", "blue3"))


library(patchwork)
fileidentity <- "monocyte_fraction_cutoff_0.000005"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(8), height = cm_to_inch(10), family = "Arial")
p + p2 + p3 + plot_layout(nrow = 3)
dev.off()

# Fig 3K ----

library(Seurat)
library(ggplot2)
library(SeuratObject)
library(patchwork)
library(fgsea)
library(dplyr)
library(fgsea)
library(showtext)
library(reshape2)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig3/monocyte_enrichment")

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

date = "20250820"
project_name <- "escc"

old <- GeomViolin$default_aes$linewidth
update_geom_defaults("boxplot", list(linewidth = 0.2))

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")
seurat_obj$classical_monocyte <- seurat_obj@assays$C2L@data["Classical.Monocyte",]
seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(classical_monocyte = case_when((is.na(classical_monocyte) & classical_monocyte < 0.00005) ~ 0,
                                                                                              T ~ classical_monocyte))

seurat_obj$non_classical_monocyte <- seurat_obj@assays$C2L@data["Non.classical.Monocyte",]
seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(non_classical_monocyte = case_when(
  (is.na(non_classical_monocyte) & non_classical_monocyte < 0.00005) ~ 0,
  T ~ non_classical_monocyte))
seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(classical_monocyte = case_when(niche %in% c("niche8") ~ classical_monocyte, T ~ NA),
                                                               non_classical_monocyte = case_when(niche %in% c("niche8") ~ non_classical_monocyte, T ~ NA)) 
seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(diff = (non_classical_monocyte - classical_monocyte)*100)

# slide_tmp <- "s1"
# q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
# p1 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")
# 
# slide_tmp <- "s2"
# q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
# p2 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")
# 
# slide_tmp <- "s3"
# q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
# p3 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")
# 
# slide_tmp <- "s4"
# q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
# p4 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")
# 
# slide_tmp <- "s5"
# q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
# p5 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")
# 
# slide_tmp <- "s6"
# q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
# p6 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")
# 
# slide_tmp <- "s7"
# q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
# p7 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")
# 
# slide_tmp <- "s8"
# q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
# p8 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")
# 
# slide_tmp <- "s9"
# q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
# p9 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")
# 
# slide_tmp <- "s10"
# q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
# p10 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")
# 
# library(patchwork)
# fileidentity <- "monocyte_fraction_diff_v2"
# cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4*5, height = 8, family = "Arial")
# p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+plot_layout(ncol = 5) & theme(legend.ticks = element_blank(), legend.text = element_blank())
# dev.off()


slide_tmp <- "s3"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p3 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp, pt.size.factor = 1.8) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

slide_tmp <- "s9"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p9 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp, pt.size.factor = 1.8) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

library(patchwork)
fileidentity <- "monocyte_fraction_diff_s3_s9_v2"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4*2, height = 8, family = "Arial")
p3+p9+plot_layout(ncol = 2) & theme(legend.ticks = element_blank(), legend.text = element_blank())
dev.off()

# Fig 3I ----
library(dplyr)
library(Seurat)
library(Mfuzz)
library(Biobase)
library(showtext)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
date <- "20250814"
project_name <- "escc"
cm_to_inch <- function(cm) cm / 2.54

## single cells ----

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/tumor_progress_genes/")
sq <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")
sq_subset <- subset(sq, subset = Annotation_v2 == "Tumor")
sq_subset$pathology <- factor(sq_subset$pathology, levels = c("Normal","Dysplasia","Microinvasive carcinoma","Macroinvasive carcinoma"))
VlnPlot(sq_subset, c("B2M","ANXA1"), group.by = "pathology")

Idents(sq_subset) <- "pathology"

Normal <- subset(sq_subset, pathology == c("Normal"));Dysplasia <- subset(sq_subset, pathology == c("Dysplasia"));
MIC <- subset(sq_subset, pathology == c("Microinvasive carcinoma"));MAC <- subset(sq_subset, pathology == c("Macroinvasive carcinoma"))
Normal_avg <- AverageExpression(Normal);Dysplasia_avg <- AverageExpression(Dysplasia)
MIC_avg <- AverageExpression(MIC);MAC_avg <- AverageExpression(MAC)
combined_df <- cbind(Normal_avg$RNA, Dysplasia_avg$RNA, MIC_avg$RNA, MAC_avg$RNA)
combined_df <- as.data.frame(combined_df)
colnames(combined_df) <- c("Normal","Dysplasia","MIC","MAC")
combined_df$gene <- rownames(combined_df) 
combined_df <- combined_df %>% dplyr::select("gene","Normal","Dysplasia","MIC","MAC")
write.table(combined_df, "SingleCell_AverageExpression_Tumor.txt", sep = "\t", quote = FALSE, row.names = F)

avg_exp <- table2eset(filename = './SingleCell_AverageExpression_Tumor.txt')
avg_exp <- filter.NA(avg_exp)
avg_exp <- fill.NA(avg_exp)
avg_exp.s <- standardise(avg_exp)

expr_mat <- exprs(avg_exp.s)
expr_mat <- expr_mat[complete.cases(expr_mat), ]
expr_mat <- expr_mat[!apply(expr_mat, 1, function(x) any(is.infinite(x))), ]
avg_exp.s <- ExpressionSet(assayData = expr_mat)
set.seed(1234); m1 <- mestimate(avg_exp.s)

# select clusters
set.seed(1234); dmin_values <- Dmin(avg_exp.s, m = m1, crange = seq(4, 40, 2), repeats = 3)
dmin_df <- data.frame(
  clusters = seq(4, 40, 2),
  min_dist  = dmin_values
)
dmin_df$first_diff <- c(NA, diff(dmin_df$min_dist))
dmin_df$second_diff <- c(NA, diff(dmin_df$first_diff))
elbow_idx <- which.max(abs(dmin_df$second_diff))
elbow_c   <- dmin_df$clusters[elbow_idx]

# library(ggplot2)
# 
# pdf("elbowplot1_hj.pdf", width = 5, height = 5)
# ggplot(dmin_df, aes(x = clusters, y = min_dist)) +
#   geom_line() +
#   geom_point() +
#   geom_vline(xintercept = elbow_c, linetype = "dashed", color = "red") +
#   annotate("text", x = elbow_c, y = max(dmin_df$min_dist),
#            label = paste("Elbow =", elbow_c),
#            vjust = -0.5, color = "red") +
#   theme_minimal() +
#   labs(title = "Dmin Curve and Elbow Point",
#        x = "Number of Clusters (c)",
#        y = "Minimum Centroid Distance")
# dev.off()

set.seed(1234); res <- mfuzz(avg_exp.s, c = 8, m = m1, iter = 100000)
saveRDS(res, "clustering_res1_hj.rds")

#pdf("./SingleCell_mfuzz_plot1.pdf",height=4*2, width=4*4)
#mfuzz.plot(avg_exp.s, cl = res, mfrow=c(4, 2),new.window = F)
#dev.off()

sel_clusters <- c(2, 5, 8)
selected_genes <- names(res$cluster[res$cluster %in% sel_clusters])
sub_eset <- avg_exp.s[selected_genes, ]
m_sub <- mestimate(sub_eset)

# select clusters
set.seed(1234); dmin_values <- Dmin(sub_eset, m = m1, crange = seq(2, 10, 1), repeats = 3)
dmin_df <- data.frame(
  clusters = seq(2, 10, 1),
  min_dist  = dmin_values
)
dmin_df$first_diff <- c(NA, diff(dmin_df$min_dist))
dmin_df$second_diff <- c(NA, diff(dmin_df$first_diff))
elbow_idx <- which.max(abs(dmin_df$second_diff))
elbow_c   <- dmin_df$clusters[elbow_idx]

library(ggplot2)

# pdf("elbowplot2_hj.pdf", width = 5, height = 5)
# ggplot(dmin_df, aes(x = clusters, y = min_dist)) +
#   geom_line() +
#   geom_point() +
#   geom_vline(xintercept = elbow_c, linetype = "dashed", color = "red") +
#   annotate("text", x = elbow_c, y = max(dmin_df$min_dist),
#            label = paste("Elbow =", elbow_c),
#            vjust = -0.5, color = "red") +
#   theme_minimal() +
#   labs(title = "Dmin Curve and Elbow Point",
#        x = "Number of Clusters (c)",
#        y = "Minimum Centroid Distance")
# dev.off()

set.seed(1234); res <- mfuzz(sub_eset, c = 20, m = m1, iter = 100000)
saveRDS(res, "clustering_res2_hj.rds")

# pdf("./SingleCell_mfuzz_plot2.pdf",height=4*5, width=4*5)
# mfuzz.plot(sub_eset, cl = res, mfrow=c(4, 5),new.window = F)
# dev.off()

sel_clusters <- c(10)
selected_genes <- names(res$cluster[res$cluster %in% sel_clusters])
saveRDS(selected_genes, "decreasing_genes_hj.rds")

sel_clusters <- c(2, 6, 9, 11, 13, 14, 19)
selected_genes <- names(res$cluster[res$cluster %in% sel_clusters])
saveRDS(selected_genes, "increasing_genes_hj.rds")

## visium ----

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/tumor_progress_genes/")
sq <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")
DefaultAssay(sq) <- "Spatial"
sq_subset <- subset(sq, subset = niche == "niche8")
sq_subset@meta.data <- sq_subset@meta.data %>% dplyr::mutate(pathology = case_when(orig.ident %in% c("s1") ~ "Normal",
                                                                                   orig.ident %in% c("s2", "s3") ~ "Dysplasia",
                                                                                   orig.ident %in% c("s4", "s5", "s6", "s7", "s8", "s9") ~ "Microinvasive carcinoma",
                                                                                   orig.ident %in% c("s10") ~ "Macroinvasive carcinoma",
))

sq_subset$pathology <- factor(sq_subset$pathology, levels = c("Normal","Dysplasia","Microinvasive carcinoma","Macroinvasive carcinoma"))

VlnPlot(sq_subset, c("B2M","ANXA1"), group.by = "pathology")

Idents(sq_subset) <- "pathology"

Normal <- subset(sq_subset, pathology == c("Normal"));Dysplasia <- subset(sq_subset, pathology == c("Dysplasia"));
MIC <- subset(sq_subset, pathology == c("Microinvasive carcinoma"));MAC <- subset(sq_subset, pathology == c("Macroinvasive carcinoma"))
Normal_avg <- AverageExpression(Normal);Dysplasia_avg <- AverageExpression(Dysplasia)
MIC_avg <- AverageExpression(MIC);MAC_avg <- AverageExpression(MAC)
combined_df <- cbind(Normal_avg$Spatial, Dysplasia_avg$Spatial, MIC_avg$Spatial, MAC_avg$Spatial)
combined_df <- as.data.frame(combined_df)
colnames(combined_df) <- c("Normal","Dysplasia","MIC","MAC")
combined_df$gene <- rownames(combined_df) 
combined_df <- combined_df %>% dplyr::select("gene","Normal","Dysplasia","MIC","MAC")
write.table(combined_df, "Visium_AverageExpression_Tumor.txt", sep = "\t", quote = FALSE, row.names = F)

avg_exp <- table2eset(filename = './Visium_AverageExpression_Tumor.txt')
avg_exp <- filter.NA(avg_exp)
avg_exp <- fill.NA(avg_exp)
avg_exp.s <- standardise(avg_exp)

expr_mat <- exprs(avg_exp.s)
expr_mat <- expr_mat[complete.cases(expr_mat), ]
expr_mat <- expr_mat[!apply(expr_mat, 1, function(x) any(is.infinite(x))), ]
avg_exp.s <- ExpressionSet(assayData = expr_mat)
set.seed(1234); m1 <- mestimate(avg_exp.s)

# select clusters
set.seed(1234); dmin_values <- Dmin(avg_exp.s, m = m1, crange = seq(4, 40, 2), repeats = 3)
dmin_df <- data.frame(
  clusters = seq(4, 40, 2),
  min_dist  = dmin_values
)
dmin_df$first_diff <- c(NA, diff(dmin_df$min_dist))
dmin_df$second_diff <- c(NA, diff(dmin_df$first_diff))
elbow_idx <- which.max(abs(dmin_df$second_diff))
elbow_c   <- dmin_df$clusters[elbow_idx]

library(ggplot2)

# pdf("visium_elbowplot1_hj.pdf", width = 5, height = 5)
# ggplot(dmin_df, aes(x = clusters, y = min_dist)) +
#   geom_line() +
#   geom_point() +
#   geom_vline(xintercept = elbow_c, linetype = "dashed", color = "red") +
#   annotate("text", x = elbow_c, y = max(dmin_df$min_dist),
#            label = paste("Elbow =", elbow_c),
#            vjust = -0.5, color = "red") +
#   theme_minimal() +
#   labs(title = "Dmin Curve and Elbow Point",
#        x = "Number of Clusters (c)",
#        y = "Minimum Centroid Distance")
# dev.off()

set.seed(1234); res <- mfuzz(avg_exp.s, c = 8, m = m1, iter = 100000)
saveRDS(res, "visium_clustering_res1_hj.rds")

# pdf("./visium_mfuzz_plot1.pdf",height=4*2, width=4*4)
# mfuzz.plot(avg_exp.s, cl = res, mfrow=c(4, 2),new.window = F)
# dev.off()

sel_clusters <- c(1, 5, 7, 8)
selected_genes <- names(res$cluster[res$cluster %in% sel_clusters])
sub_eset <- avg_exp.s[selected_genes, ]
m_sub <- mestimate(sub_eset)

# select clusters
set.seed(1234); dmin_values <- Dmin(sub_eset, m = m1, crange = seq(2, 10, 1), repeats = 3)
dmin_df <- data.frame(
  clusters = seq(2, 10, 1),
  min_dist  = dmin_values
)
dmin_df$first_diff <- c(NA, diff(dmin_df$min_dist))
dmin_df$second_diff <- c(NA, diff(dmin_df$first_diff))
elbow_idx <- which.max(abs(dmin_df$second_diff))
elbow_c   <- dmin_df$clusters[elbow_idx]

library(ggplot2)

# pdf("visium_elbowplot2_hj.pdf", width = 5, height = 5)
# ggplot(dmin_df, aes(x = clusters, y = min_dist)) +
#   geom_line() +
#   geom_point() +
#   geom_vline(xintercept = elbow_c, linetype = "dashed", color = "red") +
#   annotate("text", x = elbow_c, y = max(dmin_df$min_dist),
#            label = paste("Elbow =", elbow_c),
#            vjust = -0.5, color = "red") +
#   theme_minimal() +
#   labs(title = "Dmin Curve and Elbow Point",
#        x = "Number of Clusters (c)",
#        y = "Minimum Centroid Distance")
# dev.off()

set.seed(1234); res <- mfuzz(sub_eset, c = 20, m = m1, iter = 100000)
saveRDS(res, "visium_clustering_res2_hj.rds")

# pdf("visium_mfuzz_plot2.pdf",height=4*5, width=4*5)
# mfuzz.plot(sub_eset, cl = res, mfrow=c(4, 5),new.window = F)
# dev.off()

sel_clusters <- c(8, 10, 15)
selected_genes <- names(res$cluster[res$cluster %in% sel_clusters])
saveRDS(selected_genes, "visium_decreasing_genes_hj.rds")

sel_clusters <- c(4, 6, 12, 14)
selected_genes <- names(res$cluster[res$cluster %in% sel_clusters])
saveRDS(selected_genes, "visium_increasing_genes_hj.rds")


# overlapping genes
d1 <- readRDS("decreasing_genes_hj.rds")
d2 <- readRDS("visium_decreasing_genes_hj.rds")

d <- intersect(d1, d2)

saveRDS(d, "sc_visium_decreasing_genes_hj.rds")

i1 <- readRDS("increasing_genes_hj.rds")
i2 <- readRDS("visium_increasing_genes_hj.rds")

i <- intersect(i1, i2)

saveRDS(i, "sc_visium_increasing_genes_hj.rds")


## stacked vlnplot ----
ggplot2::update_geom_defaults("violin", list(lwd = 0.3, linewidth = 0.3))
mygenes1 <- c("SUSD4", "ARPC2", "RGS12", "SH3RF2", "SLC22A23", "NIPAL2", "ST5")
mygenes2 <- c("ADAR", "KHDC4", "COPA", "FXR1", "PCMTD1", "GRB2", "FTL")

sq <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")
sq_subset <- subset(sq, subset = Annotation_v2 == "Tumor")
sq_subset$pathology <- factor(sq_subset$pathology, levels = c("Normal","Dysplasia","Microinvasive carcinoma","Macroinvasive carcinoma"))
Idents(sq_subset) <- "pathology"

# pdf("20250815_singlecell_decreasing_vlnplot_hj.pdf", width = 4*5, height = 4*28)
# VlnPlot(sq_subset, d, group.by = "pathology", pt.size = 0, cols = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"), ncol = 5) & theme(aspect.ratio = 1) & NoLegend()
# dev.off()
# 
# pdf("20250815_singlecell_increasing_vlnplot_hj.pdf", width = 4*5, height = 4*62)
# VlnPlot(sq_subset, i, group.by = "pathology", pt.size = 0, cols = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"), ncol = 5) & theme(aspect.ratio = 1) & NoLegend()
# dev.off()
# 
# fileidentity <- "sc_decreasing"
# cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(5*4), height = cm_to_inch(5*2), family="Arial")
# VlnPlot(sq_subset, group.by = "pathology", features = mygenes1, pt.size = 0, cols = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"), ncol = 4) & ylab("Protein activity") & theme_classic(base_family = "Arial") & theme(line = element_line(linewidth = 0.3), aspect.ratio = 0.8, axis.text = element_text(size = 5), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 5, hjust = 0.5), axis.title = element_text(size = 5)) & NoLegend()
# dev.off()
# 
# fileidentity <- "sc_increasing"
# cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(5*4), height = cm_to_inch(5*2), family="Arial")
# VlnPlot(sq_subset, group.by = "pathology", features = mygenes2, pt.size = 0, cols = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"), ncol = 4) & ylab("Protein activity") & theme_classic(base_family = "Arial") & theme(line = element_line(linewidth = 0.3), aspect.ratio = 0.8, axis.text = element_text(size = 5), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 5, hjust = 0.5), axis.title = element_text(size = 5)) & NoLegend()
# dev.off()

old <- GeomViolin$default_aes$linewidth
update_geom_defaults("violin", list(linewidth = 0.2))
fileidentity <- "sc_progression_increasing_stacked_violinplot"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(4), height = cm_to_inch(6), family="Arial")
VlnPlot(sq_subset, features = mygenes2, group.by = "pathology", split.by = "pathology", stack = T, flip = T) +
  theme_classic(base_family = "Arial") +
  #ylab("Cell type") +
  xlab("") +
  theme(
    line = element_line(linewidth = 0.2),
    axis.text = element_text(colour = "black", size = 3),
    axis.text.x = element_text(size = 4, angle = 30, hjust = 1),
    axis.title = element_text(size = 4),
    panel.background = element_blank(),
    #axis.ticks = element_blank(),
    axis.ticks.length = unit(0.05, "cm"),
    #axis.ticks.y = element_blank(),
    legend.position = "None",
    legend.text = element_text(size = 3),
    strip.background = element_blank(),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 4),
    panel.grid = element_blank()) &
  scale_fill_manual(values = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"))
dev.off()


old <- GeomViolin$default_aes$linewidth
update_geom_defaults("violin", list(linewidth = 0.2))
fileidentity <- "sc_progression_decreasing_stacked_violinplot"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(4), height = cm_to_inch(6), family="Arial")
VlnPlot(sq_subset, features = mygenes1, group.by = "pathology", split.by = "pathology", stack = T, flip = T) +
  theme_classic(base_family = "Arial") +
  #ylab("Cell type") +
  xlab("") +
  theme(
    line = element_line(linewidth = 0.2),
    axis.text = element_text(colour = "black", size = 3),
    axis.text.x = element_text(size = 4, angle = 30, hjust = 1),
    axis.title = element_text(size = 4),
    panel.background = element_blank(),
    #axis.ticks = element_blank(),
    axis.ticks.length = unit(0.05, "cm"),
    #axis.ticks.y = element_blank(),
    legend.position = "None",
    legend.text = element_text(size = 3),
    strip.background = element_blank(),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 4),
    panel.grid = element_blank()) &
  scale_fill_manual(values = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"))
dev.off()

sq <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")
DefaultAssay(sq) <- "Spatial"
sq_subset <- subset(sq, subset = niche == "niche8")
sq_subset@meta.data <- sq_subset@meta.data %>% dplyr::mutate(pathology = case_when(orig.ident %in% c("s1") ~ "Normal",
                                                                                   orig.ident %in% c("s2", "s3") ~ "Dysplasia",
                                                                                   orig.ident %in% c("s4", "s5", "s6", "s7", "s8", "s9") ~ "Microinvasive carcinoma",
                                                                                   orig.ident %in% c("s10") ~ "Macroinvasive carcinoma",
))

sq_subset$pathology <- factor(sq_subset$pathology, levels = c("Normal","Dysplasia","Microinvasive carcinoma","Macroinvasive carcinoma"))

pdf("20250815_visium_decreasing_vlnplot_hj.pdf", width = 4*5, height = 4*28)
VlnPlot(sq_subset, d, group.by = "pathology", pt.size = 0, cols = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"), ncol = 5) & theme(aspect.ratio = 1) & NoLegend()
dev.off()

pdf("20250815_visium_increasing_vlnplot_hj.pdf", width = 4*5, height = 4*62)
VlnPlot(sq_subset, i, group.by = "pathology", pt.size = 0, cols = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"), ncol = 5) & theme(aspect.ratio = 1) & NoLegend()
dev.off()

VlnPlot(sq_subset, c("B2M","ANXA1"), group.by = "pathology")

fileidentity <- "visium_decreasing"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(5*4), height = cm_to_inch(5*2), family="Arial")
VlnPlot(sq_subset, group.by = "pathology", features = mygenes1, pt.size = 0, cols = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"), ncol = 4) & ylab("Protein activity") & theme_classic(base_family = "Arial") & theme(line = element_line(linewidth = 0.3), aspect.ratio = 0.8, axis.text = element_text(size = 5), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 5, hjust = 0.5), axis.title = element_text(size = 5)) & NoLegend()
dev.off()

fileidentity <- "visium_increasing"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(5*4), height = cm_to_inch(5*2), family="Arial")
VlnPlot(sq_subset, group.by = "pathology", features = mygenes2, pt.size = 0, cols = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"), ncol = 4) & ylab("Protein activity") & theme_classic(base_family = "Arial") & theme(line = element_line(linewidth = 0.3), aspect.ratio = 0.8, axis.text = element_text(size = 5), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 5, hjust = 0.5), axis.title = element_text(size = 5)) & NoLegend()
dev.off()

old <- GeomViolin$default_aes$linewidth
update_geom_defaults("violin", list(linewidth = 0.2))
fileidentity <- "spatial_progression_increasing_stacked_violinplot"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(4), height = cm_to_inch(6), family="Arial")
VlnPlot(sq_subset, features = mygenes2, group.by = "pathology", split.by = "pathology", stack = T, flip = T) +
  theme_classic(base_family = "Arial") +
  #ylab("Cell type") +
  xlab("") +
  theme(
    line = element_line(linewidth = 0.2),
    axis.text = element_text(colour = "black", size = 3),
    axis.text.x = element_text(size = 4, angle = 30, hjust = 1),
    axis.title = element_text(size = 4),
    panel.background = element_blank(),
    #axis.ticks = element_blank(),
    axis.ticks.length = unit(0.05, "cm"),
    #axis.ticks.y = element_blank(),
    legend.position = "None",
    legend.text = element_text(size = 3),
    strip.background = element_blank(),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 4),
    panel.grid = element_blank()) &
  scale_fill_manual(values = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"))
dev.off()


old <- GeomViolin$default_aes$linewidth
update_geom_defaults("violin", list(linewidth = 0.2))
fileidentity <- "spatial_progression_decreasing_stacked_violinplot"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(4), height = cm_to_inch(6), family="Arial")
VlnPlot(sq_subset, features = mygenes1, group.by = "pathology", split.by = "pathology", stack = T, flip = T) +
  theme_classic(base_family = "Arial") +
  #ylab("Cell type") +
  xlab("") +
  theme(
    line = element_line(linewidth = 0.2),
    axis.text = element_text(colour = "black", size = 3),
    axis.text.x = element_text(size = 4, angle = 30, hjust = 1),
    axis.title = element_text(size = 4),
    panel.background = element_blank(),
    #axis.ticks = element_blank(),
    axis.ticks.length = unit(0.05, "cm"),
    #axis.ticks.y = element_blank(),
    legend.position = "None",
    legend.text = element_text(size = 3),
    strip.background = element_blank(),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 4),
    panel.grid = element_blank()) &
  scale_fill_manual(values = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"))
dev.off()
