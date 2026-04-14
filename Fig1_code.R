library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(EnhancedVolcano)
library(patchwork)
library(reshape2)
library(purrr)
library(SoupX)
library(CellChat)
library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

pt_to_pixels <- function(pt_size, dpi = 600) {
  px <- (pt_size * dpi / 72)
  return(round(px, 2))
}

pt_to_geom_size <- function(pt) round(pt / 2.85, 2)

date = "20250731"
project_name = "escc"

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/v2/Fig1/20250731/")

escc <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250722_escc_total_no_cellcycle_regressout_final_hj.rds")

DimPlot(escc, label = T) & theme(aspect.ratio = 1)

escc@meta.data <- escc@meta.data %>% dplyr::mutate(Annotation = case_when(
  RNA_snn_res.0.8 %in% c("0", "12", "16", "29") ~ "B cell",
  RNA_snn_res.0.8 %in% c("13") ~ "Plasma cell",
  RNA_snn_res.0.8 %in% c("8", "15", "18", "31") ~ "Monocyte/DC/Macrophage",
  RNA_snn_res.0.8 %in% c("19") ~ "Mast",
  RNA_snn_res.0.8 %in% c("2", "7", "28", "33") ~ "T/NK cell",
  RNA_snn_res.0.8 %in% c("1", "3", "5", "17", "21", "23", "24", "26",
                         "30", "32", "34") ~ "Squamous cell",
  RNA_snn_res.0.8 %in% c("4", "11") ~ "Glandular cell",
  RNA_snn_res.0.8 %in% c("9", "20", "27") ~ "SMC/Pericyte",
  RNA_snn_res.0.8 %in% c("14", "22") ~ "Endothelial cell",
  RNA_snn_res.0.8 %in% c("6", "10", "25") ~ "Fibroblast/fDC/FRC",
  T ~ "Others"))


escc$Annotation %>% table() 
mylevel <- c("Squamous cell", "Glandular cell", "Endothelial cell", "Fibroblast/fDC/FRC", "SMC/Pericyte", "B cell", "Plasma cell", "T/NK cell", "Monocyte/DC/Macrophage", "Mast")
escc$Annotation <- factor(escc$Annotation, levels = mylevel)

# Fig1B p1 ----

size_pt <- 6
label_size <- pt_to_pixels(size_pt, dpi = 300)
ggplot_font_size <- pt_to_geom_size(pt = 6)

mycol <- c("#E41A1C", "#F781BF", "darkgreen", "#4DAF4A", "green", "#222F75","blue", "#377EB8", "#54B0E4", "cyan")
names(mycol) <- mylevel

set.seed(1234); p1 <- DimPlot(escc, group.by = "Annotation", cols = mycol, raster = T, pt.size = 1, label = F) & theme_void(base_family = "Arial") & theme_void(base_family = "Arial") & theme(plot.title = element_text(hjust = 0.5, size = 6), aspect.ratio = 1) & NoLegend()

fileidentity <- "Fig1B"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(6), height = cm_to_inch(6), family="Arial")
print(p1)
dev.off()

anno <- unique(escc@meta.data$Annotation) %>% sort()
legend_df <- data.frame(celltype = factor(anno, levels = anno), x = 1, y = 1)
p2 <- ggplot(legend_df, aes(x, y, color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = mycol) +
  theme_void(base_family = "Arial") +
  coord_cartesian(xlim = c(2,3)) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.size = unit(0.01, "cm"),
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.1, "cm"),
    legend.text = element_text(size = size_pt-3, family = "Arial")
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

fileidentity <- "Fig1B_legend"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(6), height = cm_to_inch(6), family="Arial")
p2
dev.off()

# Fig1C ----

mygenes <- c("TP63", "EGFR", "DSG3", "KRT7", "MUC5B", "AGR2", "PECAM1", "VWF", "DCN", "LUM", "FCAMR", "MYH11", "ACTA2", "LMOD1", "RGS5", "CD19", "MS4A1", "PRDM1", "XBP1", "CD2", "CD3E", "CD3G", "KLRD1", "NKG7", "CD68", "FCER1G", "ITGAX", "CLEC10A", "KIT")

escc$Annotation <- factor(escc$Annotation, levels = rev(mylevel))

fileidentity <- "Fig1C"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(11), height = cm_to_inch(5), family = "Arial")
p <- DotPlot(escc, features = mygenes, group.by = "Annotation", scale = T) +
  xlab("Genes") +
  theme_set(theme_classic(base_family = "Arial")) +
  theme(line = element_line(linewidth = 0.3),
        axis.text.y = element_text(colour = "black", size = size_pt-1),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 45, size = size_pt-1,
                                   vjust = 0.9, hjust = 0.95),
        legend.text = element_text(size = size_pt-2),
        legend.position = "bottom",
        legend.title = element_text(size = size_pt-2),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width  = unit(0.2, "cm"),
        legend.spacing.x  = unit(0.02, "cm"),
        legend.ticks = element_blank(),
        axis.title = element_text(colour = "black", size = size_pt-1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  ylab("") + scale_color_gradientn(colours=c("blue","white", "red"), name = 'Average\nExpression',limits= c(-3, 3)) +
  scale_size_area(max_size = 2) 
#guides(size = guide_legend(override.aes = list(size = 2)))
print(p)
dev.off()

# Fig1B 2 ----

df <- escc@meta.data %>% dplyr::select(Annotation, platform)
df <- df %>% dplyr::group_by(Annotation, platform) %>% dplyr::summarise(n=n())
df <- df %>% dplyr::group_by(platform) %>% dplyr::mutate(total = sum(n)) %>% dplyr::ungroup() %>%
  dplyr::mutate(prop = n/total) %>% dplyr::mutate(direction_prop = case_when(platform == "snRNA" ~ -prop,
                                                                             T ~ prop))
max_prop <- max(abs(df$direction_prop))

df$Annotation <- factor(df$Annotation, levels = rev(mylevel))

fileidentity <- "Fig1B2"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(6), height = cm_to_inch(4.5), family = "Arial")
df %>% ggplot(aes(x = Annotation, y = direction_prop,  fill = Annotation)) +
  geom_bar(stat = "identity", position = "identity", width = 0.7) +
  xlab("Celltype") +
  ylab("Proportion") +
  coord_flip() +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(labels = function(x) (abs(x)), limits = c(-max_prop, max_prop)) +
  theme_set(theme_classic(base_family = "Arial")) +
  theme(line = element_line(linewidth = 0.3),
        axis.text = element_text(size = size_pt-1, colour = "black"),
        axis.title = element_text(size = size_pt-1),
        legend.text = element_text(size = size_pt-3),
        legend.position = "None",
        legend.title = element_text(size = size_pt),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width  = unit(0.1, "cm"),
        legend.spacing.x  = unit(0.02, "cm"),
        plot.margin = unit(c(0, 0.1, 0, 0), "cm")) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))
dev.off()


# Fig1C ----

mygenes <- c("TP63", "EGFR", "DSG3", "KRT7", "MUC5B", "AGR2", "PECAM1", "VWF", "DCN", "LUM", "FCAMR", "MYH11", "ACTA2", "LMOD1", "RGS5", "CD19", "MS4A1", "PRDM1", "XBP1", "CD2", "CD3E", "CD3G", "KLRD1", "NKG7", "CD68", "FCER1G", "ITGAX", "CLEC10A", "KIT")

escc$Annotation <- factor(escc$Annotation, levels = rev(mylevel))

fileidentity <- "Fig1C"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(11), height = cm_to_inch(5), family = "Arial")
p <- DotPlot(escc, features = mygenes, group.by = "Annotation", scale = T) +
  xlab("Genes") +
  theme_set(theme_classic(base_family = "Arial")) +
  theme(line = element_line(linewidth = 0.3),
        axis.text.y = element_text(colour = "black", size = size_pt-1),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 45, size = size_pt-1,
                                   vjust = 0.9, hjust = 0.95),
        legend.text = element_text(size = size_pt-2),
        legend.position = "bottom",
        legend.title = element_text(size = size_pt-2),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width  = unit(0.2, "cm"),
        legend.spacing.x  = unit(0.02, "cm"),
        legend.ticks = element_blank(),
        axis.title = element_text(colour = "black", size = size_pt-1),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  ylab("") + scale_color_gradientn(colours=c("blue","white", "red"), name = 'Average\nExpression',limits= c(-3, 3)) +
  scale_size_area(max_size = 2) 
#guides(size = guide_legend(override.aes = list(size = 2)))
print(p)
dev.off()

# Fig1D ----

library(dplyr)
library(Seurat)
library(ggalluvial)
library(showtext)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig1/20250810/")

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

date = "20250811"
project_name <- "escc"

escc <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250731_escc_global_final_annotation_hj.rds")

escc@meta.data <- escc@meta.data %>% dplyr::mutate(Annotation = case_when(
  RNA_snn_res.0.8 %in% c("0", "12", "16", "29") ~ "B cell",
  RNA_snn_res.0.8 %in% c("13") ~ "Plasma cell",
  RNA_snn_res.0.8 %in% c("8", "15", "18", "31") ~ "Monocyte/DC/Macrophage",
  RNA_snn_res.0.8 %in% c("19") ~ "Mast",
  RNA_snn_res.0.8 %in% c("2", "7", "28", "33") ~ "T/NK cell",
  RNA_snn_res.0.8 %in% c("1", "3", "5", "17", "21", "23", "24", "26",
                         "30", "32", "34") ~ "Squamous cell",
  RNA_snn_res.0.8 %in% c("4", "11") ~ "Glandular cell",
  RNA_snn_res.0.8 %in% c("9", "20", "27") ~ "SMC/Pericyte",
  RNA_snn_res.0.8 %in% c("14", "22") ~ "Endothelial cell",
  RNA_snn_res.0.8 %in% c("6", "10", "25") ~ "Fibroblast/fDC/FRC",
  T ~ "Others"))

escc@meta.data <- escc@meta.data %>%
  dplyr::mutate(pathology = case_when(pathology == "Microinvasive" ~ "Microinvasive carcinoma",
                                      pathology == "Macroinvasive" ~ "Macroinvasive carcinoma",
                                      T ~ pathology
  ))

escc$pathology <- factor(escc$pathology, levels = c("Normal", "Dysplasia", "Microinvasive carcinoma", "Macroinvasive carcinoma"))

mylevel <- c("Squamous cell", "Glandular cell", "Endothelial cell", "Fibroblast/fDC/FRC", "SMC/Pericyte", "B cell", "Plasma cell", "T/NK cell", "Monocyte/DC/Macrophage", "Mast")
escc$Annotation <- factor(escc$Annotation, levels = mylevel)

library(ggalluvial)
table_ESCC <- escc@meta.data %>% dplyr::select(pathology, Annotation) %>% table() 
table_ESCC2 <- ( table_ESCC / table_ESCC %>% apply(1,sum) ) %>% data.frame()
table_ESCC2$pathology <- factor(table_ESCC2$pathology,  levels = c("Normal", "Dysplasia", "Microinvasive carcinoma", "Macroinvasive carcinoma"))

size_pt <- 7

mycol <- c("#E41A1C", "#F781BF", "darkgreen", "#4DAF4A", "green", "#222F75","blue", "#377EB8", "#54B0E4", "cyan")
names(mycol) <- c("Squamous cell", "Glandular cell", "Endothelial cell", "Fibroblast/fDC/FRC", "SMC/Pericyte", "B cell", "Plasma cell", "T/NK cell", "Monocyte/DC/Macrophage", "Mast")

fileidentity <- "Fig1D"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(5), height = cm_to_inch(5), family = "Arial")
ggplot(table_ESCC2, aes( x = pathology, y = Freq, alluvium = Annotation)) + 
  geom_alluvium(aes(fill = Annotation), colour = "black", alpha = 0.9, decreasing = FALSE, lwd = 0.3) + 
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

# Fig 1E ----

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")

seurat_obj$pathologist_annotation %>% unique()
dim(seurat_obj)
seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(pathologist_annotation = case_when(pathologist_annotation == "" ~ "No information",
                                                                                                  T ~ pathologist_annotation))

mycol <- scPalette(n = 21)
mylevel <- c(unique(seurat_obj$pathologist_annotation)[unique(seurat_obj$pathologist_annotation) != "No information"], "No information")
seurat_obj$pathologist_annotation <- factor(seurat_obj$pathologist_annotation, levels = mylevel)
names(mycol) <- mylevel

cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4*4, height = 5, family="Arial")
p1 <- SpatialDimPlot(seurat_obj, cols = mycol, pt.size.factor = 0, images = c("s1")) & NoLegend() & theme(aspect.ratio = 1)
p2 <- SpatialDimPlot(seurat_obj, cols = mycol, pt.size.factor = 0, images = c("s2")) & NoLegend() & theme(aspect.ratio = 1)
p3 <- SpatialDimPlot(seurat_obj, cols = mycol, pt.size.factor = 0, images = c("s4")) & NoLegend() & theme(aspect.ratio = 1)
p4 <- SpatialDimPlot(seurat_obj, cols = mycol, pt.size.factor = 0, images = c("s10")) & NoLegend() & theme(aspect.ratio = 1)
p1+p2+p3+p4+plot_layout(ncol = 4)

mycol <- c("Basal (Hyperplasia)" = "#4DAF4A", "Basal (HGD)" = "#4DAF4A", "Basal (LGD to HGD)" = "#4DAF4A", "Basal (normal)" = "#4DAF4A",
           "Dense lymphoid aggregates (r/o follicle)" =  "#E3BE00", "Germinal Center (MALT)" = "#E3BE00", "Lymphoid follicle" = "#E3BE00",
           "Invasive" = "#E41A1C",
           "Lamina propria" = "#00CDD1", "Lamina Propria (below HGD) and loose immune environments (below invasive)" = "#00CDD1",
           "Mucin ducts" = "#BC9DCC",
           "Muscularis mucosa" = "#984EA3", "Stroma (musculous)" = "#984EA3",
           "Proper muscle" = "#222F75",
           "Stroma (connective tissue)" = "#A6CEE3", "Submucosal connective tissues" = "#A6CEE3",
           "Submucosal glands" = "#A65628",
           "Suprabasal" = "#F29403",
           "Submucosal vessels" = "#E7298A", "Submucosal vessels/lymphatics" = "#E7298A",
           "No information" = "darkgrey")

cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4*4, height = 5, family="Arial")
p1 <- SpatialDimPlot(seurat_obj, group.by = "pathologist_annotation", cols = mycol, pt.size.factor = 1.8, images = c("s1")) & NoLegend() & theme(aspect.ratio = 1)
p2 <- SpatialDimPlot(seurat_obj, group.by = "pathologist_annotation", cols = mycol, pt.size.factor = 1.6, images = c("s2")) & NoLegend() & theme(aspect.ratio = 1)
p3 <- SpatialDimPlot(seurat_obj, group.by = "pathologist_annotation", cols = mycol, pt.size.factor = 1.6, images = c("s4")) & NoLegend() & theme(aspect.ratio = 1)
p4 <- SpatialDimPlot(seurat_obj, group.by = "pathologist_annotation", cols = mycol, pt.size.factor = 1.8, images = c("s10")) & NoLegend() & theme(aspect.ratio = 1)
p1+p2+p3+p4+plot_layout(ncol = 4)

# niche
mycol <- scPalette(n=12)
names(mycol) <- levels(seurat_obj$niche)

cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4*4, height = 5, family="Arial")
p1 <- SpatialDimPlot(seurat_obj, group.by = "niche", cols = mycol, pt.size.factor = 1.8, images = c("s1")) & NoLegend() & theme(aspect.ratio = 1)
p2 <- SpatialDimPlot(seurat_obj, group.by = "niche", cols = mycol, pt.size.factor = 1.6, images = c("s2")) & NoLegend() & theme(aspect.ratio = 1)
p3 <- SpatialDimPlot(seurat_obj, group.by = "niche", cols = mycol, pt.size.factor = 1.6, images = c("s4")) & NoLegend() & theme(aspect.ratio = 1)
p4 <- SpatialDimPlot(seurat_obj, group.by = "niche", cols = mycol, pt.size.factor = 1.8, images = c("s10")) & NoLegend() & theme(aspect.ratio = 1)
p1+p2+p3+p4+plot_layout(ncol = 4)

# Fig1F ----
library(dplyr)
library(tidyr)
library(purrr)
library(rstatix)
library(ggplot2)
library(scales)

seurat_obj <- readRDS('data/visium/20250821_escc_visium_c2f_scaled_latest_dl.rds')

# prep cell type composition per niche
c2l_data <- GetAssayData(seurat_obj, assay = "C2L", slot = "data")
cluster_info <- seurat_obj$niche

keep_ct <- c(
  "Tumor", "Basal", "Suprabasal", "Glandular cell",
  "Memory B", "GC DZ", "GC LZ", "Plasma cell",
  "CD8 inflammatory Trm", "CD8 CXCR6+ Trm", "CD4 Tn",
  "SMC", "NMF", "Pericyte", "BEC"
)

niche_order <- c("niche8", "niche1", "niche5", "niche4",
                 "niche12", "niche7", "niche10", "niche11",
                 "niche3", "niche9", "niche2", "niche6")

# per-niche mean abundance
df <- t(c2l_data) %>%
  as.data.frame() %>%
  mutate(cluster = cluster_info) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(-cluster, names_to = "CellType", values_to = "Expression") %>%
  filter(CellType %in% keep_ct) %>%
  mutate(
    cluster  = factor(as.character(cluster),  levels = niche_order),
    CellType = factor(as.character(CellType), levels = keep_ct)
  ) %>%
  complete(cluster, CellType, fill = list(Expression = 0))

# spot-level data for stats
df_raw <- t(c2l_data) %>%
  as.data.frame() %>%
  mutate(cluster = cluster_info) %>%
  pivot_longer(-cluster, names_to = "CellType", values_to = "value") %>%
  filter(CellType %in% keep_ct) %>%
  mutate(
    cluster  = factor(as.character(cluster),  levels = niche_order),
    CellType = factor(as.character(CellType), levels = keep_ct)
  ) %>%
  filter(!is.na(cluster), is.finite(value))

# wilcoxon: in-niche vs out-of-niche, bonferroni correction
pvals <- df_raw %>%
  group_by(CellType) %>%
  group_modify(~ {
    map_dfr(levels(.x$cluster), function(cl) {
      a <- .x %>% filter(cluster == cl)  %>% pull(value)
      b <- .x %>% filter(cluster != cl) %>% pull(value)
      
      if (length(a) >= 5 && length(b) >= 5 && (sd(a) > 0 || sd(b) > 0)) {
        wt <- wilcox_test(
          data.frame(v = c(a, b),
                     g = factor(c(rep("in", length(a)), rep("out", length(b))),
                                levels = c("in", "out"))),
          v ~ g, alternative = "greater", detailed = TRUE, exact = FALSE
        )
        tibble(cluster = cl, p = wt$p,
               median_in = median(a), median_out = median(b),
               log2FC = log2((median(a) + 1e-6) / (median(b) + 1e-6)))
      } else {
        tibble(cluster = cl, p = NA_real_,
               median_in = median(a), median_out = median(b), log2FC = NA_real_)
      }
    })
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "bonferroni")) %>%
  mutate(
    enriched = !is.na(p_adj) & p_adj < 0.05 & median_in > median_out & log2FC >= log2(1.5),
    star = case_when(
      enriched & p_adj < 0.001 ~ "**",
      enriched & p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# row z-score
df_scaled <- df %>%
  group_by(CellType) %>%
  mutate(Expression_z = if (sd(Expression, na.rm = TRUE) == 0) 0
         else as.numeric(scale(Expression))) %>%
  ungroup() %>%
  mutate(Expression_z = ifelse(is.finite(Expression_z), Expression_z, 0))

max_abs <- max(abs(df_scaled$Expression_z), na.rm = TRUE)

df_anno <- df_scaled %>%
  left_join(pvals %>% select(CellType, cluster, p_adj, star),
            by = c("CellType", "cluster")) %>%
  mutate(
    cluster  = factor(cluster,  levels = niche_order),
    CellType = factor(CellType, levels = keep_ct)
  )

# heatmap
pretty_colors <- c("#2166AC", "#4393C3", "#F7F7F7", "#FDDBC7", "#B2182B")

p <- ggplot(df_anno, aes(x = cluster, y = CellType, fill = Expression_z)) +
  geom_tile() +
  geom_tile(fill = NA, colour = "black", linewidth = 0.1) +
  scale_fill_gradientn(
    colours = pretty_colors,
    values  = scales::rescale(c(-max_abs, -max_abs/2, 0, max_abs/2, max_abs)),
    limits  = c(-max_abs, max_abs),
    name    = "Z-score"
  ) +
  geom_text(
    data = ~ filter(.x, Expression_z >= 0),
    aes(label = star), size = 2.1, color = "black"
  ) +
  scale_x_discrete(limits = niche_order) +
  scale_y_discrete(limits = keep_ct) +
  theme_minimal(base_size = 3, base_family = "Arial") +
  theme(
    axis.text.x  = element_text(size = 4, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 4),
    axis.title.x = element_text(size = 5, face = "bold"),
    axis.title.y = element_text(size = 5, face = "bold"),
    legend.title = element_text(size = 4, face = "bold"),
    legend.text  = element_text(size = 4),
    panel.grid   = element_blank(),
    plot.margin  = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  ) +
  labs(x = "Niche", y = "Cell Type")

cairo_pdf("figure/Fig3/escc_c2l_celltype_niche_heatmap.pdf",
          width = cm_to_inch(6), height = cm_to_inch(5), family = "Arial")
print(p)
dev.off()


