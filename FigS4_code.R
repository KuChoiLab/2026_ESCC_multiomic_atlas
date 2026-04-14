# FigS4A ----

library(dplyr)
library(Seurat)
library(CellChat)
library(ggplot2)
library(showtext)
library(fgsea)
library(stringr)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/cellchat/")

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

date <- "20250806"
project_name <- "escc"

cellchat_normal <- readRDS(file = "/Users/hojin/Dropbox/project/ESCC/submit/analysis/cellchat/20250806_escc_normal_hj.rds")
cellchat_dysp <- readRDS(file = "/Users/hojin/Dropbox/project/ESCC/submit/analysis/cellchat/20250806_escc_dysplasia_hj.rds")
cellchat_micro <- readRDS(file = "/Users/hojin/Dropbox/project/ESCC/submit/analysis/cellchat/20250806_escc_microinvasive_carcinoma_hj.rds")
cellchat_macro <- readRDS(file = "/Users/hojin/Dropbox/project/ESCC/submit/analysis/cellchat/20250806_escc_macroinvasive_carcinoma_hj.rds")
object.list <- list(Normal = cellchat_normal,
                    Dysplasia = cellchat_dysp,
                    Micro_invasive = cellchat_micro,
                    Macro_invasive = cellchat_macro)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4), measure = "weight")
change_p1 <- netVisual_heatmap(cellchat, comparison = c(1,2))
fileidentity <- "normal_dysplasia_cellchat_change"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 9, height = 9, family="Arial")
change_p1
dev.off()

# Fig S4B ----

cellchat_carcionma <- readRDS("~/Dropbox/projects/ESCC/submit/analysis/cellchat/20250823_escc_micro_macro_merged_hj.rds")
object.list <- list(Dysplasia = cellchat_dysp,
                    Carcionma = cellchat_carcionma)

change_p1 <- netVisual_heatmap(cellchat, comparison = c(1,2))
mat <- change_p1@matrix

groupSize <- mat[, colnames(mat) == "Tumor"]
groupSize[groupSize < 0 | is.na(groupSize) | is.infinite(groupSize)] <- 0
mat[, colnames(mat) != "Tumor"] <- 0
mat[is.na(mat)] <- 0

par(family = "Arial", mar = c(0, 0, 0, 0))
netVisual_circle(mat, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, vertex.label.cex = 0.2) 

# Fig S4C ----

## * part1 (in myCAF - COL, LAM family) ----
library(scCustomize)
library(Seurat)
library(ggplot2)

my_colors <- c(
  "Normal"= "#83f52c",
  "Dysplasia" = "#f3f315",
  "Microinvasive carcinoma" = "#ff6600",
  "Macroinvasive carcinoma" = "#ff0099"
)

fibroblast_final_sk <- readRDS("~/Dropbox/20250724_ESCC_scRNAseq_annotation_T_cells/RData/20250812_fibroblast_final_sk.rds")
fibroblast_final_sk$pathology <- factor(fibroblast_final_sk$pathology, levels = c('Normal','Dysplasia','Microinvasive','Macroinvasive'))
levels(fibroblast_final_sk$pathology) <- c('Normal','Dysplasia','Microinvasive carcinoma','Macroinvasive carcinoma')
table(fibroblast_final_sk$Annotation_v2)
# APOD high pFib            iCAF           myCAF             NAF             NMF progenitor iCAF TIMP3 high pFib 
# 538            1374             340             829            4096             874             756
myCAF_seurat <- subset(fibroblast_final_sk, Annotation_v2 %in% "myCAF")

genelist_S2B_COL <- c('COL1A1',
                      'COL1A2',
                      'COL4A1',
                      'COL4A2',
                      'COL4A5',
                      'COL6A1',
                      'COL6A2',
                      'COL6A3')

genelist_S2B_LAM <- c('LAMA2',
                      'LAMA4')

genelist_S2B_ITG <- c('ITGA2',
                      'ITGA6',
                      'ITGB1',
                      'ITGB4')

part1_genes <- c(genelist_S2B_COL, genelist_S2B_LAM)

Stacked_VlnPlot(myCAF_seurat, features = part1_genes, colors_use = my_colors, group.by = 'pathology')

Idents(myCAF_seurat) <- myCAF_seurat$pathology
markers_results_part1 <- FindMarkers(
  myCAF_seurat,
  ident.1 = "Dysplasia", 
  ident.2 = "Normal",    
  features = part1_genes,
  logfc.threshold = 0,   
  min.pct = 0
)

markers_results_part1[part1_genes,]

## * part 2 (in tumor - ITG family) ----
squamous_final_sk <- readRDS("~/Dropbox/20250724_ESCC_scRNAseq_annotation_T_cells/RData/20250812_squamous_final_sk.rds")

tumor_rna_seurat <- subset(squamous_final_sk, Annotation_v3 %in% "Tumor")
Idents(tumor_rna_seurat) <- tumor_rna_seurat$pathology

table(tumor_rna_seurat$pathology)

fileidentity <- "S4B_cellchat_genes_VlnPlot_combined_part2_updated" 
cairo_pdf(paste0("Figure4_sk/Supplementary_4B/","20250819" ,"_", project_name, "_", fileidentity, "_", myname, ".pdf"), width = cm_to_inch(10), height = cm_to_inch(40), family = "Arial") # $)C3m9. GG0\ Fd@LAv4B 0!7N 18cm, <<7N 22.5cm 0! CV4k@T4O4Y. 
Stacked_VlnPlot(tumor_rna_seurat, features = genelist_S2B_ITG, colors_use = my_colors)
dev.off()

markers_results_part2 <- FindMarkers(
  tumor_rna_seurat,
  ident.1 = "Dysplasia",
  ident.2 = "Normal",   
  features = part2_genes,
  logfc.threshold = 0,  
  min.pct = 0
)

markers_results_part2[genelist_S2B_ITG,]

# Fig S4D ----
library(dplyr)
library(Seurat)
library(CellChat)
library(ggplot2)
library(showtext)
library(stringr)
library(GLAD)
library(RColorBrewer)
library(viridis)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

date <- "20250823"
project_name <- "escc"

cellchat_normal <- readRDS(file = "Figure4_sk/Supplementary_cellchat//20250806_escc_normal_hj.rds")
cellchat_dysp <- readRDS(file = "Figure4_sk/Supplementary_cellchat//20250806_escc_dysplasia_hj.rds")
cellchat_micro <- readRDS(file = "Figure4_sk/Supplementary_cellchat//20250806_escc_microinvasive_carcinoma_hj.rds")
cellchat_macro <- readRDS(file = "Figure4_sk/Supplementary_cellchat//20250806_escc_macroinvasive_carcinoma_hj.rds")
cellchat_carcinoma <- readRDS(file = "Figure4_sk/Supplementary_cellchat/20250823_escc_micro_macro_merged_hj.rds")

object.list <- list(Normal = cellchat_normal,
                    Dysplasia = cellchat_dysp,
                    Carcinoma = cellchat_carcinoma)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

df.net_carc <- subsetCommunication(object.list$Carcinoma, slot.name = "net", 
                                   sources.use = "myCAF", targets.use = "Tumor")
df.net_carc$stage <- "Carcinoma"

df.net_carc = df.net_carc %>% dplyr::select(ligand, prob, receptor,stage,interaction_name_2)
min_prob = min(df.net_carc$prob)
max_prob = max(df.net_carc$prob)
df.net_carc = df.net_carc %>% dplyr::mutate(prob_norm = (prob - min_prob)/(max_prob - min_prob))

df.net_dys = df.net_dys %>% dplyr::select(ligand, prob, receptor,stage,interaction_name_2)
min_prob = min(df.net_dys$prob)
max_prob = max(df.net_dys$prob)
df.net_dys = df.net_dys %>% dplyr::mutate(prob_norm = (prob - min_prob)/(max_prob - min_prob))

mydf = rbind(df.net_carc, df.net_dys)

common_cci = mydf %>% dplyr::group_by(interaction_name_2) %>% dplyr::summarise(n=n()) %>% dplyr::filter(n == 2) %>% dplyr::pull(interaction_name_2)
private_cci = mydf %>% dplyr::group_by(interaction_name_2) %>% dplyr::summarise(n=n()) %>% dplyr::filter(n == 1) %>% dplyr::pull(interaction_name_2)

mydf = mydf %>% dplyr::mutate(Pair_info = case_when(
  interaction_name_2 %in% common_cci ~ "Common pair",
  interaction_name_2 %in% private_cci ~ "Private pair"))
mydf$stage <- factor(mydf$stage, levels = c("Dysplasia", "Carcinoma"))

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

sc <- scale_fill_gradientn(colours = myPalette(100))


ligand_level <- unique(sort(mydf$ligand))
mydf$ligand <- factor(mydf$ligand, levels = rev(ligand_level))

fileidentity <- "dysplasia_carcinoma_cellchat_dotplot"
cairo_pdf(paste0("Figure4_sk/Supplementary_cellchat/",date, "_", project_name, "_", fileidentity, "_sk.pdf"), width = cm_to_inch(10), height = cm_to_inch(10), family="Arial")
mydf %>%
  ggplot(aes(y = ligand, x = receptor, fill = prob_norm), color = "black") +
  geom_point(size = 1.5, aes(shape = Pair_info), stroke = 0.2) +
  facet_wrap(~ stage, ncol = 2) + 
  labs(y = "Ligand", x = "Receptor") +
  scale_shape_manual(values = c(21,24)) + 
  theme_bw(base_family = "Arial")+
  scale_y_discrete(position = "right") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(colour = "black", size = 4),
        aspect.ratio = 2,
        #line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.3),
        legend.ticks = element_blank(),
        legend.text = element_text(size = 3),
        #legend.title = element_text(size = size_pt-2),
        legend.key.height = unit(0.2, "cm"),
        panel.grid = element_line(linewidth = 0.3),
        legend.key.width  = unit(0.2, "cm"),
        legend.spacing.x  = unit(0.02, "cm"),
        strip.background = element_rect(color = "black", fill = "gray95", linewidth = 0.5),
        panel.border = element_rect(color = "black", linewidth = 0.5),
        strip.text = element_text(size = 5),
        axis.title = element_text(colour = "black", size = 4),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0.05, unit = "cm"),
        legend.title = element_blank()) + sc
dev.off()

# Fig S4E ----

library(Seurat)

# pull C2L abundances for blend
features_c2l <- c(
  "Tumor", "myCAF", "iCAF", "Progenitor iCAF",
  "Basal", "Suprabasal", "Glandular cell",
  "GC DZ", "CD4 Tpex", "SMC", "NMF",
  "CD8 HSP+ Trm", "CD8 CXCR6+ Trm",
  "Proliferating T", "Proliferating macrophage",
  "Classical Monocyte", "Non classical Monocyte"
)

features_spatial <- c(
  "NOTCH1", "CD4",
  "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"
)

blend_mat <- rbind(
  merged@assays$C2L[features_c2l, ],
  merged@assays$Spatial[features_spatial, ]
)

# clean rownames (no spaces or special chars)
rownames(blend_mat) <- c(
  "Tumor", "myCAF", "iCAF", "ProgenitoriCAF",
  "Basal", "Suprabasal", "GlandularCell",
  "GCDZ", "CD4Tpex", "SMC", "NMF",
  "CD8HSPTrm", "CD8CXCR6Trm",
  "ProliferatingT", "ProliferatingMac",
  "ClassicalMono", "NonClassicalMono",
  "NOTCH1", "CD4",
  "HLADMA", "HLADMB", "HLADPA1", "HLADPB1", "HLADQA1", "HLADQB1"
)

merged[["blend"]] <- CreateAssayObject(blend_mat)
merged@assays$blend@data <- blend_mat
merged@active.assay <- "blend"

# blend spatial plot — Tumor / myCAF / iCAF / Progenitor iCAF
cairo_pdf("figure/Fig4/escc_tumor_mycaf_icaf_progenitoriCAF_blend.pdf",
          width = cm_to_inch(30), height = cm_to_inch(60), family = "Arial")

SpatialFeaturePlotBlend3(
  merged,
  features = c("Tumor", "myCAF", "iCAF", "ProgenitoriCAF"),
  alpha = c(1, 1),
  ncol = 2,
  pt.size = 1.9
) &
  guides(
    fill   = guide_colourbar(barheight = unit(3, "mm"), barwidth = unit(20, "mm")),
    colour = guide_colourbar(barheight = unit(3, "mm"), barwidth = unit(20, "mm"))
  ) &
  theme(
    aspect.ratio  = 1,
    legend.text   = element_text(size = 6),
    legend.title  = element_text(size = 7),
    plot.title    = element_text(size = 9)
  )

dev.off()

# Fig S4F ----

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(AUCell)
library(SCENIC)
library(ComplexHeatmap)
library(ggplot2)
library(stringr)
library(showtext)

seurat_obj  <- readRDS(file.path(base_dir, "20250812_squamous_final_sk.rds"))
regulonAUC  <- readRDS(file.path(base_dir, "3.4_regulonAUC.Rds"))
res_df      <- readRDS(file.path(base_dir, "tumor_stage_test_hj.rds"))

res_df <- res_df %>%
  dplyr::filter(group1 == "Dysplasia") %>%
  dplyr::mutate(
    padj    = p.adjust(pval, method = "bonferroni"),
    log2FC  = log2(group1_mean / group2_mean),
    TF_name = str_extract(TF, "^[^ ]+")
  )

res_df_filtered <- res_df %>%
  dplyr::filter(
    group2_mean < group1_mean,   # higher in Dysplasia
    padj        < 0.05,
    log2FC      > 0
  )

pathology_levels <- c("Normal", "Dysplasia", "Microinvasive carcinoma", "Macroinvasive carcinoma")

meta <- seurat_obj@meta.data %>%
  dplyr::filter(Annotation_v2 == "Tumor") %>%
  dplyr::mutate(Pathology = factor(pathology, levels = pathology_levels)) %>%
  dplyr::select(Pathology)

regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
regulonAUC <- regulonAUC[, rownames(meta)]

auc_mat <- getAUC(regulonAUC)
auc_scaled <- t(scale(t(auc_mat), center = TRUE, scale = TRUE))   # z-score per regulon

cells_by_path <- split(rownames(meta), meta$Pathology)
cell_order    <- unlist(cells_by_path[pathology_levels])

auc_scaled_ordered <- auc_scaled[res_df_filtered$TF, cell_order, drop = FALSE]

pth_df <- data.frame(
  cell      = rownames(meta),
  Pathology = meta$Pathology,
  row.names = NULL
)

summ <- as.data.frame(auc_scaled_ordered) %>%
  rownames_to_column("TF") %>%
  pivot_longer(-TF, names_to = "cell", values_to = "score") %>%
  left_join(pth_df, by = "cell") %>%
  dplyr::group_by(TF, Pathology) %>%
  dplyr::summarise(
    mean_score = mean(score, na.rm = TRUE),
    n          = dplyr::n(),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(
    TF        = factor(TF, levels = rev(rownames(auc_scaled_ordered))),  # match heatmap row order
    Pathology = factor(Pathology, levels = pathology_levels)
  )

p <- ggplot(summ, aes(x = Pathology, y = TF)) +
  geom_point(aes(size = n, colour = mean_score)) +
  scale_size_continuous(range = c(2, 9), name = "Cells (n)") +
  scale_colour_gradient2(
    low      = "#1f3a93",
    mid      = "white",
    high     = "#ff6b6b",
    midpoint = 0,
    name     = "Mean score"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y  = element_text(size = 7),
    axis.title   = element_blank(),
    legend.title = element_text(size = 7),
    legend.text  = element_text(size = 6),
    panel.grid.major.x = element_line(colour = "grey88"),
    panel.grid.major.y = element_line(colour = "grey88")
  )


# Fig S4G ----
intersect_genes <- readRDS("~/Dropbox/20250724_ESCC_scRNAseq_annotation_T_cells/RData/intersect_genes.rds")

## * RNA ----
squamous_final_sk <- readRDS("~/Dropbox/20250724_ESCC_scRNAseq_annotation_T_cells/RData/20250812_squamous_final_sk.rds")
tumor_rna_seurat <- subset(squamous_final_sk, Annotation_v3 %in% "Tumor")
Idents(tumor_rna_seurat) <- tumor_rna_seurat$pathology
table(tumor_rna_seurat$Annotation_v2)
# Basal Suprabasal      Tumor 
# 14         20       8846
tumor_only_rna_seurat <- subset(tumor_rna_seurat, Annotation_v2 %in% 'Tumor')
Stacked_VlnPlot(tumor_only_rna_seurat, features = intersect_genes, colors_use = my_colors)

SGE_markers_RNA <- FindMarkers(
  tumor_only_rna_seurat,
  ident.1 = "Dysplasia",
  ident.2 = "Normal",   
  features = intersect_genes,
  logfc.threshold = 0,  
  min.pct = 0
)

SGE_markers_RNA[intersect_genes, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    avg_log2FC = round(avg_log2FC, 2),                           
    p_val_adj  = formatC(p_val_adj, format = "E", digits = 2)    
  )

## * protein ----
escc <- readRDS("RData/escc_tumor_viper_hj.rds")
Idents(escc) <- escc$pathology

DefaultAssay(escc) <- "viper"
Stacked_VlnPlot(escc, features = intersect_genes, colors_use = my_colors) + theme(text = element_text(family = "Arial"))

SGE_markers_protein <- FindMarkers(
  escc,
  ident.1 = "Dysplasia",  
  ident.2 = "Normal",     
  features = intersect_genes,
  logfc.threshold = 0,    
  min.pct = 0
)

SGE_markers_protein[intersect_genes, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  mutate(
    p_val_adj  = formatC(p_val_adj, format = "E", digits = 2) 
  )

## viper protein activity mean values

library(dplyr)
library(Seurat)
library(viper)
library(ggplot2)
library(RColorBrewer)
library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
date <- "20250810"
project_name <- "escc"
cm_to_inch <- function(cm) cm / 2.54

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/viper/tumor/")

## make viper objects 
escc <- readRDS("~/Dropbox/project/ESCC/submit/object_share/20250806_escc_squamouscells_umap_tuning_sk_hj.rds")
escc <- subset(escc, subset = Annotation_v2 == "Tumor")
escc$donor %>% table()
cell1 <- escc@meta.data %>% dplyr::filter(pathology == "Normal") %>% rownames()

set.seed(1234)
escc$cellname <- colnames(escc)
cell2 <- escc@meta.data %>%
  filter(pathology != "Normal") %>%
  group_by(pathology) %>%
  slice_sample(prop = 0.1, replace = FALSE) %>% dplyr::pull(cellname)

cell3 <- c(cell1, cell2)
saveRDS(cell3, "random_cells_hj.rds")

cell3 <- readRDS("random_cells_hj.rds")

escc <- subset(escc, cells = cell3)

myscobj <- escc 
mymat <- myscobj@assays$RNA@data %>% as.matrix()
set.seed(1234)
mymat <- as.data.frame(mymat)
mymat$gene <- rownames(mymat)
mycol = c("gene", colnames(mymat))
mycol = mycol[1:length(mycol)-1]
mymat <- mymat[,mycol]
write.table(mymat, "tumor_total_exp_mat.txt", col.names = T, sep = "\t", row.names = F, quote = F)

## Regulon analysis was conducted in mccleary 

exp.mat = read.table("tumor_total_exp_mat.txt", header = T, check.names = F)
rownames(exp.mat) <- exp.mat$gene
exp.mat <- exp.mat %>% dplyr::select(-c("gene"))
exp.mat <- as.matrix(exp.mat)
reg = readRDS("./tumor_v2/Pruned.rds")
viperes = viper(eset = exp.mat, regulon = reg, verbose = T)
viperes %>% dim()

dim(escc)
assay_tmp <- CreateAssayObject(data = viperes)
escc[['viper']] <- assay_tmp

saveRDS(escc, file = "escc_tumor_viper_hj.rds")

rm(escc)
rm(exp.mat)

escc_subset <- readRDS(file = "escc_tumor_viper_hj.rds")

## normal vs dys
pathology_target = c("Dysplasia")
viper_mat <- as.matrix(escc_subset[["viper"]]@data)
colnames(viper_mat) <- colnames(escc_subset)
metadat <- escc_subset@meta.data
metadat$cell_name <- colnames(escc_subset)

wilres_df <- data.frame()
for (i in 1:nrow(viper_mat)) {
  reg_tmp <- rownames(viper_mat)[i]
  df_tmp <- as.data.frame(viper_mat[reg_tmp,])
  colnames(df_tmp) <- "protein_activity"
  
  df_tmp <- merge(df_tmp, metadat, by.x = "row.names", by.y = "cell_name")
  x <- df_tmp %>% dplyr::filter(pathology %in% pathology_target) %>% dplyr::pull(protein_activity)
  y <- df_tmp %>% dplyr::filter(pathology %in% "Normal") %>% dplyr::pull(protein_activity)
  wilres <- wilcox.test(x = x, y = y)
  
  wilres_df_tmp <- data.frame(TF = reg_tmp,
                              mean_group1 = mean(x),
                              mean_group2 = mean(y),
                              diff = mean(x) - mean(y),
                              pvalue = wilres$p.value)
  wilres_df <- rbind(wilres_df, wilres_df_tmp)
}

wilres_df$adj_pvalue <- p.adjust(wilres_df$pvalue, method = "bonferroni")

write.table(wilres_df, "20250807_escc_viper_dys_normal_wilcoxon_hj.txt", col.names = T, sep = "\t", row.names = F, quote = F)


# Fig S4H ----
library(Seurat)
library(dplyr)
library(AUCell)
library(SCENIC)
library(ComplexHeatmap)
library(stringr)
library(showtext)
library(ComplexHeatmap)
library(CellChat)
library(colorRamp2)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

setwd(paste0("/Users/hojin/Dropbox/project/ESCC/submit/analysis/scenic/fib/"))

## 1. wilcoxon ----

wilcoxon_res <- list()
regulonAUC <- readRDS('/Users/hojin/Dropbox/project/ESCC/submit/analysis/scenic/fib/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

## 2. scaling ----
## sc ------
seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/analysis/step02_subclustering/fibroblast_fDC_SMC_pericyte/20250728_escc_Fib_fDC_FRC_hj.rds")

seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(Annotation_v2 = case_when(Annotation_v2 == "APOD high pFib" ~ "APOD high CAF",
                                                                                         Annotation_v2 == "TIMP3 high pFib" ~ "CFD high CAF",
                                                                                         T ~ Annotation_v2))

seurat_obj <- subset(seurat_obj, subset = Annotation_v2 %in% c("NMF", "progenitor iCAF", "iCAF", "myCAF"))

seurat_obj <- subset(seurat_obj, subset = platform == "scRNA")

metadat <- seurat_obj@meta.data
regulonmat <- regulonAUC@assays@data$AUC
regulonmat <- regulonmat[, colnames(seurat_obj)]
TFs <- rownames(regulonmat) %>% unique()
res_df <- data.frame()
TFs <- rownames(regulonmat) %>% unique()
metadat <- seurat_obj@meta.data
metadat$cellname <- rownames(metadat)

res_df <- data.frame()
for (i in 1:length(TFs)) {
  tf_tmp <- TFs[i]
  react <- regulonmat[tf_tmp, ]
  df_tmp <- data.frame(cellname = names(react), regulon_activity = react)
  df_tmp <- merge(df_tmp, metadat, by.x = "cellname", by.y = "cellname", all.x = T)
  
  for (celltype_tmp in c(c(seurat_obj$Annotation_v2 %>% unique()))) {
    df1 <- df_tmp %>% dplyr::filter(Annotation_v2 == celltype_tmp)
    df2 <- df_tmp %>% dplyr::filter(Annotation_v2 != celltype_tmp)
    res <- wilcox.test(df1$regulon_activity, df2$regulon_activity)
    pval_tmp <- res$p.value
    res_tmp <- data.frame(TF = tf_tmp,
                          group1 = celltype_tmp,
                          group1_mean = mean(df1$regulon_activity),
                          group2_mean = mean(df2$regulon_activity),
                          pval = pval_tmp)
    res_tmp$padj <- p.adjust(res_tmp$pval, method = "bonferroni")
    res_df <- rbind(res_df, res_tmp)
  }
}

saveRDS(res_df, "fib_sc_test_hj.rds")
res_df <- readRDS("fib_sc_test_hj.rds")


## sn ------
seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/analysis/step02_subclustering/fibroblast_fDC_SMC_pericyte/20250728_escc_Fib_fDC_FRC_hj.rds")

seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(Annotation_v2 = case_when(Annotation_v2 == "APOD high pFib" ~ "APOD high CAF",
                                                                                         Annotation_v2 == "TIMP3 high pFib" ~ "CFD high CAF",
                                                                                         T ~ Annotation_v2))

seurat_obj <- subset(seurat_obj, subset = Annotation_v2 %in% c("NMF", "progenitor iCAF", "iCAF", "myCAF"))

seurat_obj <- subset(seurat_obj, subset = platform == "snRNA")

metadat <- seurat_obj@meta.data
regulonmat <- regulonAUC@assays@data$AUC
regulonmat <- regulonmat[, colnames(seurat_obj)]
TFs <- rownames(regulonmat) %>% unique()
res_df <- data.frame()
TFs <- rownames(regulonmat) %>% unique()
metadat <- seurat_obj@meta.data
metadat$cellname <- rownames(metadat)

res_df <- data.frame()
for (i in 1:length(TFs)) {
  tf_tmp <- TFs[i]
  react <- regulonmat[tf_tmp, ]
  df_tmp <- data.frame(cellname = names(react), regulon_activity = react)
  df_tmp <- merge(df_tmp, metadat, by.x = "cellname", by.y = "cellname", all.x = T)
  
  for (celltype_tmp in c(c(seurat_obj$Annotation_v2 %>% unique()))) {
    df1 <- df_tmp %>% dplyr::filter(Annotation_v2 == celltype_tmp)
    df2 <- df_tmp %>% dplyr::filter(Annotation_v2 != celltype_tmp)
    res <- wilcox.test(df1$regulon_activity, df2$regulon_activity)
    pval_tmp <- res$p.value
    res_tmp <- data.frame(TF = tf_tmp,
                          group1 = celltype_tmp,
                          group1_mean = mean(df1$regulon_activity),
                          group2_mean = mean(df2$regulon_activity),
                          pval = pval_tmp)
    res_tmp$padj <- p.adjust(res_tmp$pval, method = "bonferroni")
    res_df <- rbind(res_df, res_tmp)
  }
}

saveRDS(res_df, "fib_sn_test_hj.rds")
res_df <- readRDS("fib_sn_test_hj.rds")

res_sc_df <- readRDS("fib_sc_test_hj.rds")
res_sc_df$platform <- "sc"
res_sn_df <- readRDS("fib_sn_test_hj.rds")
res_sn_df$platform <- "sn"

res_df <- rbind(res_sc_df, res_sn_df)
quant <- quantile(res_df$padj, 0.1)

filtered_df <- res_df %>%
  dplyr::group_by(group1, platform) %>%
  dplyr::filter(padj < quantile(padj, 0.15, na.rm = TRUE)) %>%
  ungroup()

res_df <- filtered_df %>% dplyr::mutate(log2FC = log2(group1_mean/group2_mean)) %>% dplyr::filter(group2_mean < group1_mean) %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FC > 0.8)

res_counts <- res_df %>% dplyr::group_by(TF, group1) %>% dplyr::summarise(n=n()) %>% dplyr::arrange(desc(n))
mytfs <- res_counts %>% dplyr::filter(n == 2) %>% dplyr::pull(TF) %>% unique()

mycount_df <- res_counts %>% dplyr::filter(n == 2)

## 3. merge sc, sn ----

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/analysis/step02_subclustering/fibroblast_fDC_SMC_pericyte/20250728_escc_Fib_fDC_FRC_hj.rds")

seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(Annotation_v2 = case_when(Annotation_v2 == "APOD high pFib" ~ "APOD high CAF",
                                                                                         Annotation_v2 == "TIMP3 high pFib" ~ "CFD high CAF",
                                                                                         T ~ Annotation_v2))

seurat_obj <- subset(seurat_obj, subset = Annotation_v2 %in% c("NMF", "progenitor iCAF", "iCAF", "myCAF"))

seurat_obj <- subset(seurat_obj, subset = platform == "scRNA")
regulonAUC_scaled_sc <- getAUC(regulonAUC)[mytfs, colnames(seurat_obj)]

regulonAUC_scaled_sc <- t(scale(t(regulonAUC_scaled_sc), center = T, scale=T))
# regulonAUC_scaled_sc <- t(apply(regulonAUC_scaled_sc, 1, function(x) {
#   (x - min(x)) / (max(x) - min(x))
# }))
#regulonAUC_scaled_sc <- t(scale(t(getAUC(regulonAUC)[res_df_filtered$TF, colnames(seurat_obj)]), center = T, scale=T))

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/analysis/step02_subclustering/fibroblast_fDC_SMC_pericyte/20250728_escc_Fib_fDC_FRC_hj.rds")

seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(Annotation_v2 = case_when(Annotation_v2 == "APOD high pFib" ~ "APOD high CAF",
                                                                                         Annotation_v2 == "TIMP3 high pFib" ~ "CFD high CAF",
                                                                                         T ~ Annotation_v2))

seurat_obj <- subset(seurat_obj, subset = Annotation_v2 %in% c("NMF", "progenitor iCAF", "iCAF", "myCAF"))

seurat_obj <- subset(seurat_obj, subset = platform == "snRNA")
regulonAUC_scaled_sn <- getAUC(regulonAUC)[mytfs, colnames(seurat_obj)]
regulonAUC_scaled_sn <- t(scale(t(regulonAUC_scaled_sn), center = T, scale=T))
regulonAUC_scaled <- cbind(regulonAUC_scaled_sn, regulonAUC_scaled_sc)

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/analysis/step02_subclustering/fibroblast_fDC_SMC_pericyte/20250728_escc_Fib_fDC_FRC_hj.rds")
seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(Annotation_v2 = case_when(Annotation_v2 == "APOD high pFib" ~ "APOD high CAF",
                                                                                         Annotation_v2 == "TIMP3 high pFib" ~ "CFD high CAF",
                                                                                         T ~ Annotation_v2))

seurat_obj <- subset(seurat_obj, subset = Annotation_v2 %in% c("NMF", "progenitor iCAF", "iCAF", "myCAF"))

celltypes <- seurat_obj@meta.data %>% select(Annotation_v2)
celltypes$New_Identity <- factor(celltypes$Annotation_v2, levels = c("NMF", "progenitor iCAF", "iCAF", "myCAF"))
selectedResolution <- "New_Identity"
cellsPerTypes <- split(rownames(celltypes), celltypes[,selectedResolution]) 

## 4. scenic heatmap -----
cellorder <- unlist(cellsPerTypes)

cell_type <- data.frame(
  "cell type" = unlist(mapply(rep, names(cellsPerTypes), lapply(cellsPerTypes, length))))

regulonAUC_scaled_order <- regulonAUC_scaled[mytfs, cellorder]
row.names(cell_type) <- colnames(regulonAUC_scaled_order)
palette_length = 100

date <- "20250806"
project_name <- "escc"
cm_to_inch <- function(cm) cm / 2.54
fileidentity <- "fib_scenic"

library(CellChat)
mycol <- scPalette(4)
names(mycol) <- names(cellsPerTypes)

platform_col <- c("scRNA" = "grey", "snRNA" = "orange3")
platform_vec <- seurat_obj@meta.data[colnames(regulonAUC_scaled_order), "platform"] |> as.character()
platform_vec <- factor(platform_vec, levels = names(platform_col))

ha = HeatmapAnnotation(Platform = platform_vec,
                       Celltype = unlist(mapply(rep, names(cellsPerTypes), lapply(cellsPerTypes, length))),
                       col = list(
                         Platform = platform_col,
                         Celltype  = mycol
                       ),
                       simple_anno_size = unit(0.3, "cm"),
                       annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param = list(
                         title_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                         labels_gp = gpar(fontsize = 5, fontfamily = "Arial")
                       ))

col_fun = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "red"))

png(filename = paste0(date, "_", project_name, "_", fileidentity, "_hj.png"),
    width = 8, height = 5, units = "in", res = 600, family = "Arial")
set.seed(1234);ht <- Heatmap(matrix = regulonAUC_scaled_order, name = "mat",
                             col = col_fun,
                             cluster_columns = F,
                             cluster_column_slices = F,
                             clustering_method_rows = "ward.D2",
                             show_column_names = F,
                             use_raster = TRUE,
                             raster_device = "png",
                             raster_quality = 5,
                             raster_device_param = list(res = 600),
                             row_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
                             column_names_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                             heatmap_legend_param = list(
                               title = "Regulon activity",
                               title_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                               labels_gp = gpar(fontsize = 5, fontfamily = "Arial")
                             ),
                             top_annotation = ha)
draw(ht, use_raster = TRUE)
dev.off()

cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"),
          width = cm_to_inch(10), height = cm_to_inch(7), family = "Arial")

col_fun <- colorRamp2(c(-1, 0, 1), c("darkblue", "white", "red"))

set.seed(1234); ht <- Heatmap(
  matrix = regulonAUC_scaled_order,
  name = "mat",
  col  = col_fun,
  use_raster = TRUE,    
  raster_device = "png",
  raster_quality = 5,
  raster_device_param = list(res = 600),
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  clustering_method_rows = "ward.D2",
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 5, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "Arial"),
  heatmap_legend_param = list(
    title = "Score",
    title_gp = gpar(fontsize = 4, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 4, fontfamily = "Arial")
  ),
  top_annotation = ha
)

draw(ht)
dev.off()


## 5. ORA ----
setwd(paste0("/Users/hojin/Dropbox/project/ESCC/submit/analysis/scenic/fib/"))

library(org.Hs.eg.db)
library(msigdbr)
library(clusterProfiler)

regulons <- readRDS('./2.6_regulons_asGeneSet.Rds')
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, entrez_gene)
hs <- org.Hs.eg.db
row_order_vec <- row_order(ht)
ordered_row_names <- rownames(regulonAUC_scaled_order)[row_order_vec]

tfs <- gsub(ordered_row_names, pattern = " \\([0-9]*g\\)", replacement = "")

for (i in 1:length(tfs)) {
  if ("COL1A1" %in% regulons[[tfs[i]]]) {
    print(tfs[i])
  }
}

enricher_res <- data.frame()

for (i in 1:length(tfs)) {
  
  print(i)
  target_genes <- regulons[[tfs[i]]]
  print(i)
  ENTREZID <- AnnotationDbi::select(hs, keys = target_genes, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>% dplyr::select(ENTREZID) %>% na.omit() %>% pull()
  
  em <- tryCatch({
    res <- clusterProfiler::enricher(ENTREZID,
                                     TERM2GENE = m_t2g,
                                     pvalueCutoff = 1,
                                     qvalueCutoff = 1,
                                     minGSSize = 1)
  }, error = function(e) {
    message("enricher error (", tfs[i], "): ", e$message)
    return(NULL)
  })
  
  if (!is.null(em) && nrow(em) > 0) {
    em <- as.data.frame(em)
    em$TFs <- tfs[i]
    enricher_res <- rbind(enricher_res, em)
  } else {
    "complete"
  }
  
}

library(WriteXLS)
WriteXLS(enricher_res, "20250807_escc_fib_targetgenes_gobp_hj.xlsx")

enricher_res <- readxl::read_excel("20250807_escc_fib_targetgenes_gobp_hj.xlsx")
enricher_res <- enricher_res %>% dplyr::filter(p.adjust < 0.05)
enricher_res <- enricher_res %>% dplyr::group_by(TFs) %>% dplyr::arrange(p.adjust, .by_group = T) %>% dplyr::slice(1)
rownames(enricher_res) <- enricher_res$TFs

tf_filtered <- intersect(tfs, enricher_res$TFs)
enricher_res <- enricher_res[tf_filtered,] # sorting top1 results


WriteXLS(enricher_res, "20250807_escc_fib_targetgenes_gobp_top1_hj.xlsx")

## 6. GOBP plot ----

library(readxl)

enricher_res <- read_xlsx("20250807_escc_fib_targetgenes_gobp_top1_hj.xlsx")

library(plyr)
myplot_df <- data.frame()

for (tf_tmp in tfs) {
  
  nrow_res <- enricher_res[enricher_res$TFs == tf_tmp,] %>% nrow()
  
  if (nrow_res == 0) {
    
    myplot_df_tmp <- data.frame(TFs = tf_tmp)
    myplot_df <- rbind.fill(myplot_df, myplot_df_tmp)
    
  }
  
  myplot_df_tmp <- enricher_res[enricher_res$TFs == tf_tmp,]
  myplot_df <- rbind(myplot_df, myplot_df_tmp)
  
}

myplot_df$X = 1
myplot_df$TFs <- factor(myplot_df$TFs, levels = myplot_df$TFs)

myplot_df <- myplot_df %>%
  dplyr::mutate(color = case_when(is.na(ID) ~ "white",
                                  T ~ "black")) %>%
  dplyr::mutate(ID = case_when(is.na(ID) ~ TFs,
                               T ~ ID))

myplot_df <- myplot_df %>% dplyr::select(ID, p.adjust)

df <- myplot_df %>%
  transmute(
    ID    = as.character(ID),
    padj  = suppressWarnings(as.numeric(p.adjust))
  )

use_log <- TRUE
val <- if (use_log) -log10(df$padj) else df$padj
val_name <- if (use_log) "-log10(padj)" else "padj"

mat <- as.matrix(val)
colnames(mat) <- val_name
rownames(mat) <- make.unique(df$ID, sep = " ") 
row_lbl <- df$ID                               

rng <- range(val, na.rm = TRUE)
if (use_log) {
  col_fun <- colorRamp2(c(0, rng[2]), c("white", "red"))
} else {
  col_fun <- colorRamp2(c(rng[1], rng[2]), c("white", "red"))
}

rownames(mat) <- gsub(rownames(mat), pattern = "GOBP_", replacement = "") %>% gsub(., pattern = "_", replacement = " ")
row_lbl <- gsub(row_lbl, pattern = "GOBP_", replacement = "") %>% gsub(., pattern = "_", replacement = " ")
row_lbl <- paste0(
  toupper(substr(row_lbl, 1, 1)),
  tolower(substr(row_lbl, 2, nchar(row_lbl)))
)

ht <- Heatmap(
  mat,
  name = val_name,
  col  = col_fun,
  na_col = "darkgrey",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_labels = row_lbl,
  row_names_gp = gpar(fontsize = 5, fontfamily = "Arial"),
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(col = "black", fill = NA, lwd = 0.5))
  },
  width = unit(ncol(mat)/5.1, "cm"), 
  height = unit(nrow(mat)/5.1, "cm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 3, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 3, fontfamily = "Arial"),
    legend_height = unit(1, "cm"),                        
    grid_height = unit(0.2, "cm"),                        
    grid_width = unit(0.1, "cm")                          
  )
)

# save plots
fileidentity <- "scenic_top1_gobp"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"),
          width = cm_to_inch(10), height = cm_to_inch(8), family = "Arial")
draw(
  ht,
  heatmap_legend_side = "bottom"
)
dev.off()


# Fig S4I ----

library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(rstatix)
library(ggh4x)

# 1. Sample group annotation
visium$sample_group <- NA
visium$sample_group[visium$orig.ident == "s1"]                          <- "Normal"
visium$sample_group[visium$orig.ident %in% c("s2", "s3")]              <- "Dysplasia"
visium$sample_group[visium$orig.ident %in% c("s4","s5","s6","s7","s8","s9")] <- "Microinvasive carcinoma"
visium$sample_group[visium$orig.ident == "s10"]                         <- "Macroinvasive carcinoma"

visium$sample_group <- factor(
  visium$sample_group,
  levels = c("Normal", "Dysplasia", "Microinvasive carcinoma", "Macroinvasive carcinoma")
)

# 2. Extract feature vectors
col1a1_vec <- visium@assays$Spatial@data["COL1A1", , drop = TRUE]
mycaf_vec  <- visium@assays$C2L@data["myCAF",   , drop = TRUE]

common_bc  <- intersect(names(col1a1_vec), names(mycaf_vec))
col1a1_vec <- col1a1_vec[common_bc]
mycaf_vec  <- mycaf_vec[common_bc]

# 3. Build data frame
meta <- visium@meta.data[common_bc, , drop = FALSE] %>%
  tibble::rownames_to_column("barcode")

df <- tibble::tibble(
  barcode = common_bc,
  COL1A1  = as.vector(col1a1_vec),
  myCAF   = as.vector(mycaf_vec)
) %>%
  left_join(meta %>% dplyr::select(barcode, sample_group), by = "barcode") %>%
  filter(!is.na(sample_group)) %>%
  mutate(
    sample_group = factor(
      sample_group,
      levels = c("Normal", "Dysplasia", "Microinvasive carcinoma", "Macroinvasive carcinoma")
    )
  )

# 4. Correlation per group (Pearson, BH-adjusted)
cor_by_group <- df %>%
  dplyr::group_by(sample_group) %>%
  rstatix::cor_test(myCAF, COL1A1, method = "pearson") %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  dplyr::mutate(label = sprintf("r = %.2f\np = %.2g", cor, p.adj)) %>%
  dplyr::ungroup()

# 5. Color palette
col_conditions <- c(
  "Normal"                  = "#83f52c",
  "Dysplasia"               = "#f3f315",
  "Microinvasive carcinoma" = "#ff6600",
  "Macroinvasive carcinoma" = "#ff0099"
)

# 6. Plot
p_scatter_group <- ggplot(df, aes(x = myCAF, y = COL1A1, color = sample_group)) +
  geom_point(alpha = 0.6, size = 0.8, stroke = 0) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "red") +
  ggh4x::facet_wrap2(
    ~ sample_group, nrow = 1, scales = "free_y",
    strip = ggh4x::strip_themed(
      text_x = list(element_text(face = "bold", size = size_pt, color = "black"))
    )
  ) +
  scale_color_manual(values = col_conditions) +
  labs(
    x     = "myCAF abundance (Cell2Location)",
    y     = "COL1A1 expression",
    color = "Condition"
  ) +
  geom_text(
    data = cor_by_group %>% dplyr::mutate(x = Inf, y = -Inf),
    aes(x = x, y = y, label = label),
    hjust = 1.05, vjust = -0.5,
    inherit.aes = FALSE,
    size   = 4,
    family = "Arial"
  ) +
  theme_classic(base_size = size_pt, base_family = "Arial") +
  theme(
    aspect.ratio = 1,
    text         = element_text(family = "Arial"),
    line         = element_line(linewidth = 0.3),
    axis.text    = element_text(size = size_pt + 2, color = "black"),
    strip.text   = element_text(face = "bold", size = size_pt + 2),
    legend.position = "none"
  )

print(p_scatter_group)

# Fig S4J ----

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

escc_global_final_annotation_hj <- readRDS("RData/20250731_escc_global_final_annotation_hj.rds")
escc_global_final_annotation_hj$pathology <- factor(escc_global_final_annotation_hj$pathology, levels = c('Normal','Dysplasia','Microinvasive','Macroinvasive'))
levels(escc_global_final_annotation_hj$pathology) <- c('Normal','Dysplasia','Microinvasive carcinoma','Macroinvasive carcinoma')

table(escc_global_final_annotation_hj$pathology)
# Normal               Dysplasia Microinvasive carcinoma Macroinvasive carcinoma 
# 48500                   43948                   33892                   24683

unique(escc_global_final_annotation_hj$Annotation_v2)

# 1. cell type order
cell_order <- c(
  # lymphocyte
  "naive B", "Memory B", "GC LZ", "GC DZ", "Plasma cell",
  "CD4 Tn", "CD4 Tfh", "CD4 Tcm", "CD4 Early Tcm", "CD4 Activated Tcm",
  "CD4 CTL", "CD4 differentiating Treg", "CD4 TNFRSF9- Treg", "CD4 TNFRSF9+ Treg",
  "CD4 Terminal Tex", "CD4 Tpex", "Mito-high CD4", "Proliferating T",
  "CD8 Trm", "CD8 CXCR6+ Trm", "CD8 HSP+ Trm", "CD8 Inflammatory Trm",
  "CD8 GZMK+ early Tem", "CD8 Temra", "CD8 Terminal Tex", "Mito-high CD8",
  "NK1", "NK2", "NK3",
  # myeloid
  "Classical Monocyte", "Non-classical Monocyte", 
  "C1QC+ TAM", "LA TAM", "Interstitial macrophage", "Proliferating macrophage",
  "CLEC9A+ cDC1", "CD1C+ cDC2", "LAMP3+ cDC3", "pDC", "Mast cell",
  # fibroblast
  "NMF", "NAF", "APOD+ CAF", "CCBE1+ CAF", "progenitor iCAF", "iCAF", "myCAF",
  # squamous
  "Basal", "Suprabasal", "Tumor"
)

# 2. immune related genes order
immune_related_genes <- c(
  "NKG7", "GZMA", "GZMB", "GNLY", "PRF1", "GZMK", "IFNG", # Cytotoxic/effector
  "CD28", "ICOS", "CD27", "CD40LG", "TNFRSF9", "TNFRSF4", "TNFRSF18", "TNFRSF14", # Co-stimulatory
  "CD274", "PDCD1LG2", "IDO1", "LILRB4", "LILRB2", "TDO2", "CD276", "VSIR", "LGALS9", "LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4", "BTLA", "KLRC1", "ENTPD1", "LAYN" # Co-inhibitory/exhuastion
)

# 3. subset data
obj_sub <- subset(
  escc_global_final_annotation_hj, 
  Annotation_v2 %in% cell_order
)

obj_sub$Annotation_v2 <- factor(obj_sub$Annotation_v2, levels = cell_order)

# 4. transform to dataframe
df <- FetchData(
  obj_sub,
  vars = c("Annotation_v2", immune_related_genes, "pathology")
)

# 5. transform to longer form
df_long <- df %>%
  pivot_longer(
    cols = all_of(immune_related_genes),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  mutate(gene = factor(gene, levels = immune_related_genes))

# 6. plotting
p <- ggplot(df_long, aes(x = Annotation_v2, y = expression, fill = pathology)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_grid(gene ~ ., scales = "free_y", space = "free_y", switch = "y" ) +
  scale_fill_manual(values = c(
    "Normal"= "#83f52c",
    "Dysplasia" = "#f3f315",
    "Microinvasive carcinoma" = "#ff6600",
    "Macroinvasive carcinoma" = "#ff0099"
  )) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y.left = element_text(angle = 0),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside"
  ) +
  labs(x = NULL, y = NULL)
p
