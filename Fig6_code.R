
# Fig 6A / 6D ----

library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

# GC B fraction per sample
gc_b_values <- GetAssayData(merged.gcb, assay = "C2L", slot = "data")["GC B", ]

df_sample <- data.frame(
  sample_id          = merged.gcb$orig.ident,
  pathological_stage = merged.gcb$pathological_stage,
  GC_B_fraction      = as.numeric(gc_b_values)
) %>%
  group_by(sample_id, pathological_stage) %>%
  summarise(GC_B_fraction = mean(GC_B_fraction), .groups = "drop") %>%
  mutate(sample_id = trimws(as.character(sample_id)))

# merge with colocalization result
df_plot <- df_sample %>%
  left_join(summary_df %>% mutate(sample = trimws(as.character(sample))),
            by = c("sample_id" = "sample")) %>%
  rename(Colocalization_fraction = 6) %>%
  mutate(pathological_stage = factor(pathological_stage,
                                     levels = c("Normal", "Dysplasia", "Microinvasive carcinoma", "Macroinvasive carcinoma")))

stage_pal <- c(
  "Normal"                  = "#83f52c",
  "Dysplasia"               = "#f3f315",
  "Microinvasive carcinoma" = "#ff6600",
  "Macroinvasive carcinoma" = "#ff0099"
)

base_theme <- theme_classic(base_family = "Arial") +
  theme(
    axis.text.x  = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
    axis.text.y  = element_text(size = 8),
    axis.title.y = element_text(size = 9, face = "bold"),
    plot.margin  = unit(c(0.1, 0.2, 0.1, 0.1), "cm")
  )

# GC B fraction boxplot (log10)
p_box <- ggplot(df_plot, aes(x = pathological_stage, y = GC_B_fraction,
                             fill = pathological_stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, size = 1, shape = 21, stroke = 0.5) +
  scale_fill_manual(values = stage_pal) +
  scale_y_continuous(trans = "log10", labels = percent_format(accuracy = 0.1)) +
  labs(x = "", y = "GC B fraction") +
  base_theme +
  theme(legend.position = "none")

# TLS colocalization line plot
prop_sum <- df_plot %>%
  group_by(pathological_stage) %>%
  summarise(
    mean = mean(Colocalization_fraction, na.rm = TRUE),
    se   = sd(Colocalization_fraction, na.rm = TRUE) / sqrt(sum(!is.na(Colocalization_fraction))),
    .groups = "drop"
  )

p_line <- ggplot(prop_sum, aes(x = pathological_stage, y = mean, color = pathological_stage)) +
  geom_line(aes(group = 1), linewidth = 0.8, color = "black") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15, na.rm = TRUE) +
  scale_color_manual(values = stage_pal, name = "Pathological stage") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), limits = c(0, NA)) +
  labs(x = "", y = "TLS component\ncolocalization fraction") +
  base_theme +
  theme(legend.position = "right",
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.2), "cm"))

cairo_pdf("figure/Fig6/escc_Fig6A_GCB_fraction_colocalization.pdf",
          width = cm_to_inch(16), height = cm_to_inch(8), family = "Arial")
p_box | p_line
dev.off()

# Fig 6B ----

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(CellChat)
library(showtext)
library(tidyr)
library(tibble)
library(circlize)
library(tidyverse)
library(ComplexHeatmap)
library(colorRamp2)
library(vegan)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

date = "20250809"
project_name = "escc"

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig6/")

seurat_obj <- readRDS("~/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")

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
my_TLS <- c("niche10")

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

saveRDS(res_df, "20250810_escc_TLS_enrichment_hj.rds")
res_df <- readRDS("20250810_escc_TLS_enrichment_hj.rds")



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

bin_mat <- (mat > 0) * 1
bin_mat[is.na(bin_mat)] <- 0
row_dend_bin <- hclust(vegan::vegdist(bin_mat, method = "jaccard", binary = TRUE), method = "centroid")
row_groups <- cutree(row_dend_bin, k = 14) 
set.seed(1234); ht <- Heatmap(
  t(mat),
  name = "log2FC",
  col = col_fun,
  cluster_columns = as.dendrogram(row_dend_bin),
  column_split    = 14,
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


fileidentity <- "TLS_enrichment_log2FC_celltype_heatmap_split"
cairo_pdf(paste0("20250818", "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(15) , height = cm_to_inch(8), family = "Arial")
draw(
  ht,
  heatmap_legend_list = list(lg),  
  heatmap_legend_side = "right",
  padding = unit(c(2, 2, 2, 2), "mm")
)
dev.off()

# Fig 6C ----

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/final_figure_codes/Fig6/Fig6A/")

seurat_obj <- readRDS("~/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250821_escc_visium_c2f_scaled_latest_dl.rds")

seurat_obj$Tumor <- as.vector(seurat_obj@assays$C2L["Tumor",])
seurat_obj$GCB <- as.vector(seurat_obj@assays$C2L["GC DZ",]) + as.vector(seurat_obj@assays$C2L["GC LZ",])

SpatialFeaturePlotBlend <- function(object, features, combine = TRUE,
                                    feature_1_alt_name = NULL,
                                    feature_2_alt_name = NULL, assay = NULL,
                                    bottom_left = "#000000",
                                    bottom_right = "green1", # #ff0000
                                    top_left = "blue1", # green1, blue1
                                    top_right = "yellow1", ...)  {
  
  # Generate a grid of RGB color values given the requested corner colours.
  gen_color_grid <- function(side_length, bottom_left, bottom_right,
                             top_left, top_right) {
    
    grad_gen <- function(start, end, n = side_length) {
      colfunc <- colorRampPalette(c(start, end))
      return(colfunc(n))
    }
    
    # x_y = "x to y"; "bl" = "bottom left", etc
    bl_tl <- grad_gen(bottom_left, bottom_right)
    br_tr <- grad_gen(top_left, top_right)
    
    l <- lapply(seq_len(length(bl_tl)),
                function(i) {
                  start <- bl_tl[i]
                  end <- br_tr[i]
                  new_grad <- grad_gen(start, end)
                })
    
    return(t(matrix(unlist(l), ncol = side_length, nrow = side_length)))
  }
  
  if (length(features) != 2) {
    stop(paste(c("Incorrect number of features. ",
                 "Requires two features, received ",
                 length(features))))
  }
  
  if (!is.null(assay)) {
    DefaultAssay(object) <- assay
  }
  
  blend_plot_theme <- theme(legend.position = "none",
                            plot.title = element_text(hjust = 0.5))
  
  plot_list_outer <- list()
  
  for (i in Images(object)) {
    plot_list <- lapply(features,
                        function(feature) {
                          max_color <- ifelse(feature == features[1],
                                              bottom_right, top_left)
                          SpatialFeaturePlot(object, feature,
                                             images = i, ...) +
                            scale_fill_gradient(low = bottom_left,
                                                high = max_color) +
                            ggtitle(feature) +
                            blend_plot_theme
                        })
    
    cell_barcodes <- Seurat:::CellsByImage(object, images = i,
                                           unlist = TRUE)
    cells_obj_sub <- object[, cell_barcodes]
    images_sub_list <- list(object[[i]])
    names(images_sub_list) <- i
    cells_obj_sub@images <- images_sub_list
    dat <- FetchData(cells_obj_sub, features)
    side_length <- 100
    col_grid <- gen_color_grid(side_length, bottom_left, bottom_right,
                               top_left, top_right)
    dat_norm <- apply(dat, 2,
                      function(x) {
                        round((side_length - 1) * x / max(x)) + 1
                      })
    colors <- sapply(seq_len(nrow(dat_norm)), function(x) {
      col_grid[dat_norm[x, 1], dat_norm[x, 2]]
    })
    new_md_column <- paste0(features[1], "_vs_", features[2])
    cells_obj_sub[[new_md_column]] <- colors
    names(colors) <- as.character(colors)
    
    plot_list[[3]] <- SpatialDimPlot(cells_obj_sub, new_md_column,
                                     cols = colors, images = i, ...) +
      ggtitle(paste0(features[1], "_", features[2])) +
      blend_plot_theme
    legend_grid <- expand.grid(
      seq(from = min(dat[, features[1]], na.rm = TRUE),
          to   = max(dat[, features[1]], na.rm = TRUE),
          length.out = side_length),
      seq(from = min(dat[, features[2]], na.rm = TRUE),
          to   = max(dat[, features[2]], na.rm = TRUE),
          length.out = side_length)
    )
    colnames(legend_grid) <- features
    legend_colors <- c(col_grid)
    legend_grid$color <- legend_colors
    names(legend_colors) <- legend_colors
    
    legend <- ggplot(legend_grid,
                     aes(x = .data[[features[1]]], y = .data[[features[2]]],
                         color = color)) +
      geom_point(shape = 15, size = 1.9) +
      scale_color_manual(values = legend_colors) +
      coord_cartesian(expand = FALSE) +
      theme(legend.position = "none", aspect.ratio = 1,
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      xlab(ifelse(is.null(feature_1_alt_name),
                  features[1], feature_1_alt_name)) +
      ylab(ifelse(is.null(feature_2_alt_name),
                  features[2], feature_2_alt_name))
    
    plot_list[[4]] <- wrap_plots(ggplot() + theme_void(), legend,
                                 ggplot() + theme_void(), ncol = 1,
                                 heights = c(0.2, 0.6, 0.2))
    
    plot_list_outer[[i]] <- plot_list
  }
  
  if (combine == FALSE) {
    return(plot_list_outer)
  } else {
    plot_list_outer <- lapply(plot_list_outer,
                              function(p) {
                                wrap_plots(p, nrow = 1,
                                           widths = c(0.28, 0.28,
                                                      0.28, 0.16))
                              })
    p <- wrap_plots(plot_list_outer, ncol = 1)
    
    return(p)
  }
}

SpatialFeaturePlotBlend(seurat_obj, features = c("Tumor", "GCB"), bottom_right = "magenta", top_left = "cyan", top_right = scales::colour_ramp(c("magenta", "cyan"))(0.5), pt.size.factor = 1.8)


# Fig 6E ----

library(dplyr)
library(Seurat)
library(viper)
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(showtext)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

project_name <- "escc"
setwd("~/Dropbox/project/ESCC/submit/analysis/viper/Immune_cells")

## 1. make viper objects -------
escc <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250731_escc_global_final_annotation_hj.rds")
escc$Annotation_v2 %>% table()

escc$Annotation %>% table()
#escc <- subset(escc, subset = Annotation %in% c("B cell", "Myeloid cell", "Plasma cell", "T/NK cell"))
#escc <- subset(escc, subset = Annotation_v2 != ".")
#escc <- subset(escc, subset = Annotation_v2 != "Mito-high CD4")
#escc <- subset(escc, subset = Annotation_v2 != "Mito-high CD8")

escc_v2 <- subset(escc, subset = Annotation_v2 %in% c("naive B"))
set.seed(1234); samples <- sample(colnames(escc_v2), size = 3000)
escc_v2 <- subset(escc, subset = Annotation_v2 %in% c("GC DZ", "GC LZ", "Plasma cell", "Classical Monocyte"))
escc$barcode <- colnames(escc)
escc <- subset(escc, subset = barcode %in% c(colnames(escc_v2), samples))
escc$Annotation_v2 %>% table()

myscobj <- escc 
mymat <- myscobj@assays$RNA@data %>% as.matrix()
set.seed(1234)
mymat <- as.data.frame(mymat)
mymat$gene <- rownames(mymat)
mycol = c("gene", colnames(mymat))
mycol = mycol[1:length(mycol)-1]
mymat <- mymat[,mycol]
write.table(mymat, "GCB_positive_control_total_exp_mat.txt", col.names = T, sep = "\t", row.names = F, quote = F)


## Regulon analysis was conducted in mccleary ----

exp.mat = read.table("GCB_positive_control_total_exp_mat.txt", header = T, check.names = F)
rownames(exp.mat) <- exp.mat$gene
exp.mat <- exp.mat %>% dplyr::select(-c("gene"))
exp.mat <- as.matrix(exp.mat)
#reg = readRDS("~/Dropbox/project/ESCC/submit/analysis/viper/Immune_cells/Immune_cells/Pruned.rds")
reg = readRDS("~/Dropbox/project/ESCC/submit/analysis/viper/Immune_cells/Immune_cells/unPruned.rds")
viperes = viper(eset = exp.mat, regulon = reg, verbose = T)
viperes %>% dim()

dim(escc)
assay_tmp <- CreateAssayObject(data = viperes)
escc[['viper']] <- assay_tmp

escc@assays$viper@data["IRF8",]

saveRDS(escc, file = "escc_gcb_viper_hj.rds")

rm(escc)
rm(exp.mat)

## 2. viper wilcoxon (GCB) -------
escc_subset <- readRDS(file = "escc_gcb_viper_hj.rds")
escc_subset <- subset(escc_subset, subset = Annotation_v2 %in% c("GC DZ", "GC LZ"))

# normal vs dys
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

wilres_df$adj_pvalue <- p.adjust(wilres_df$pvalue, method = "BH")
wilres_df %>% dplyr::filter(TF == "IRF8")
write.table(wilres_df, "20251222_escc_gcb_viper_dys_normal_wilcoxon_hj.txt", col.names = T, sep = "\t", row.names = F, quote = F)

# micro vs macro
pathology_target = c("Microinvasive")
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
  y <- df_tmp %>% dplyr::filter(pathology %in% "Macroinvasive") %>% dplyr::pull(protein_activity)
  wilres <- wilcox.test(x = x, y = y)
  
  wilres_df_tmp <- data.frame(TF = reg_tmp,
                              mean_group1 = mean(x),
                              mean_group2 = mean(y),
                              diff = mean(x) - mean(y),
                              pvalue = wilres$p.value)
  wilres_df <- rbind(wilres_df, wilres_df_tmp)
}

wilres_df$adj_pvalue <- p.adjust(wilres_df$pvalue, method = "BH")

write.table(wilres_df, "20251222_escc_gcb_viper_micro_macro_wilcoxon_hj.txt", col.names = T, sep = "\t", row.names = F, quote = F)

tb5 = read.table("20251222_escc_gcb_viper_dys_normal_wilcoxon_hj.txt", header = T, sep = "\t")
tb6 = read.table("20251222_escc_gcb_viper_micro_macro_wilcoxon_hj.txt", header = T, sep = "\t")
tb5 <- tb5 %>% dplyr::filter(adj_pvalue < 0.05 & diff > 0.5)
tb6 <- tb6 %>% dplyr::filter(adj_pvalue < 0.05 & diff > 0.5) 

mygenes <- intersect(tb5$TF, tb6$TF)
mygenes # "CDCA7L" "ELF1"   "IRF8"   "LMO2"   "RBBP7"  "RBM39"  "REL"

mygenes <- tb5 %>% dplyr::filter(TF %in% mygenes) %>% dplyr::arrange(desc(diff)) %>% dplyr::pull(TF)

escc_subset <- readRDS(file = "escc_gcb_viper_hj.rds")

escc_subset@meta.data <- escc_subset@meta.data %>% dplyr::mutate(pathology = case_when(pathology == "Macroinvasive" ~ "Macroinvasive carcinoma",
                                                                                       pathology == "Microinvasive" ~ "Microinvasive carcinoma",
                                                                                       T ~ pathology))
escc_subset@meta.data$pathology <- factor(escc_subset@meta.data$pathology, levels = c("Normal", "Dysplasia", "Microinvasive carcinoma", "Macroinvasive carcinoma"))

escc_subset@meta.data <- escc_subset@meta.data %>% dplyr::mutate(Annotation_v2 = case_when(Annotation_v2 == "GC LZ" ~ "GC B",
                                                                                           Annotation_v2 == "GC DZ" ~ "GC B",
                                                                                           Annotation_v2 == "naive B" ~ "Naive B",
                                                                                           T ~ Annotation_v2
))

Idents(escc_subset) <- "Annotation_v2"
DefaultAssay(escc_subset) <- "viper"


escc_subset$Annotation_v2 <- factor(escc_subset$Annotation_v2, levels = c("GC B", "Naive B", "Plasma cell", "Classical Monocyte"))
p <- FetchData(
  escc_subset,
  vars = c(mygenes, "Annotation_v2", "pathology")
) %>%
  tidyr::pivot_longer(
    cols = all_of(mygenes),
    names_to = "gene",
    values_to = "Protein activity"
  )

p$gene <- factor(p$gene, levels = mygenes)

pdf("20251222_escc_viper_GCB_positive_controls_hj.pdf", width = 5, height = 6)
ggplot(p, aes(x = Annotation_v2, y = `Protein activity`, fill = pathology)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7,
    linewidth = 0.3) +
  xlab("Cell type") +
  facet_grid(gene ~ ., scales = "free_y", switch = "y") +
  theme_classic(base_family = "Arial") +
  scale_fill_manual(values = c("#83f52c", "#f3f315", "#ff6600", "#ff0099")) +
  theme(
    text = element_text(family = "Arial"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top")
dev.off()

# Fig 3F ----

library(readxl)
library(survival)
library(ggplot2)
library(survminer)
library(ggsci)
library(GSVA)
library(dplyr)
library(showtext)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

axis_theme = theme(axis.text.x = element_text(size=10, color = "black"),
                   axis.title.x = element_text(size=10, color = "black"),
                   axis.text.y = element_text(size=10, color = "black"),
                   axis.title.y = element_text(size=10, color = "black"), panel.border = element_rect(size = 1))

### input check
Original_TLS_df = read_xlsx("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/final_figure_codes/Fig6/Fig6F/TLS_data.xlsx", sheet = "Sheet1")
rownames(Original_TLS_df) = Original_TLS_df$Samples
Original_TLS_df = Original_TLS_df[,-1]

### output check
work_dir = "."


### TLS area
TLS_area_df = Original_TLS_df[,grepl("TLS.area", colnames(Original_TLS_df))]
colnames(TLS_area_df) = gsub(".TLS.area", "", colnames(TLS_area_df))
TLS_area_df = apply(TLS_area_df, 2, function(x) {x[x=="no TLS"] = 0;return(as.numeric(x))}) %>% data.frame()
rownames(TLS_area_df) = rownames(Original_TLS_df)
colnames(TLS_area_df) = c("Normal", "Low Grade\nDysplasia", "High Grade\nDysplasia", "Microinvasive\nCarcinoma", "Macroinvasive\nCarcinoma")
TLS_area_df$sample_name = rownames(TLS_area_df)


TLS_area_df_melt = reshape2::melt(TLS_area_df)

### total area
total_area_df = Original_TLS_df[,grepl("total.area", colnames(Original_TLS_df))]
colnames(total_area_df) = gsub(".total.area", "", colnames(total_area_df))
colnames(total_area_df) = c("Normal", "Low Grade\nDysplasia", "High Grade\nDysplasia", "Microinvasive\nCarcinoma", "Macroinvasive\nCarcinoma")
total_area_df$sample_name = rownames(total_area_df)
total_area_df_melt = reshape2::melt(total_area_df)

identical(paste0(total_area_df_melt$sample_name,"_", total_area_df_melt$variable), paste0(TLS_area_df_melt$sample_name,"_",TLS_area_df_melt$variable))

TLS_density = TLS_area_df_melt$value/total_area_df_melt$value

TLS_density = data.frame("variable" = TLS_area_df_melt$variable, "value" = TLS_density)

low_wilcox = wilcox.test(TLS_density[TLS_density$variable %in% c("Normal", "Low Grade\nDysplasia"),"value"],
                         TLS_density[TLS_density$variable %in% c("High Grade\nDysplasia"),"value"])
high_wilcox = wilcox.test(TLS_density[TLS_density$variable %in% c("High Grade\nDysplasia"),"value"],
                          TLS_density[TLS_density$variable %in% c("Microinvasive\nCarcinoma", "Macroinvasive\nCarcinoma"),"value"])


p = ggplot(TLS_density, aes(x = variable, y = value, fill=variable)) + geom_boxplot(outlier.shape = NA, width=0.7) + geom_point(position = position_jitter(seed = 42), alpha = 0.6)
p = p + theme_bw() + axis_theme + ggtitle(paste0("low: ", low_wilcox$p.value,"\nhigh: ", high_wilcox$p.value)) + guides(fill="none") + xlab("") + ylab("TLS density") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), text = element_text(family = "Arial")) +
  scale_fill_manual(values = c("Normal"= "#83f52c", "Low Grade\nDysplasia" = "ivory", "High Grade\nDysplasia" = "#f3f315", "Microinvasive\nCarcinoma" = "#ff6600", "Macroinvasive\nCarcinoma" = "#ff0099"))

pdf("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/final_figure_codes/Fig6/Fig6F/20251225_escc_FigF_hj.pdf", width = 5, height = 5)
print(p)
dev.off()

# Fig 6G ----

library(Seurat)
library(openxlsx)
library(survival)
library(ggplot2)
library(survminer)
library(ggsci)
library(GSVA)
library(dplyr)
library(SoupX)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/survival_analysis")

tumor <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250806_escc_squamouscells_umap_tuning_sk_hj.rds")
escc <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250731_escc_global_final_annotation_hj.rds")
escc$key <- colnames(escc)

##################################
thickness = theme(legend.title = element_text(size = 15, color = "black"),
                  legend.text = element_text(size = 15, color = "black"),
                  axis.text.x = element_text(size = 15, color = "black"),
                  axis.text.y = element_text(size = 15, color = "black"),
                  axis.title.x = element_text(size = 15, color = "black"),
                  axis.title.y = element_text(size = 15, color = "black"),
                  plot.title = element_text(size = 20, color = "black"))
ESCA = read.table("/Users/hojin/Downloads/TCGA.ESCA.sampleMap_HiSeq.gz", header = TRUE)
colnames(ESCA) = gsub("\\.", "-", colnames(ESCA))
rownames(ESCA) = ESCA[, 1]
ESCA = ESCA[,seq(2, ncol(ESCA))]
clinical = read.csv("/Users/hojin/Downloads/clinical.tsv", sep = '\t', header = TRUE)
clinical = clinical[!duplicated(clinical[,2]),]
clinical = clinical[grepl("quamous",clinical$primary_diagnosis),]
clinical$'pathologic stage' = clinical$ajcc_pathologic_stage
clinical$'pathologic stage'[clinical$'pathologic stage' %in% paste0("Stage ", c("IA", "IB", "II", "IIA", "IIB"))] = "I & II"
clinical$'pathologic stage'[clinical$'pathologic stage' %in% paste0("Stage ", c("III", "IIIA", "IIIB", "IIIC", "IV", "IVA"))] = "III & IV"
clinical$'pathologic stage'[clinical$'pathologic stage' == "'--"] = NA
clinical$ajcc_pathologic_stage3 = clinical$ajcc_pathologic_stage
clinical$ajcc_pathologic_stage3[clinical$ajcc_pathologic_stage3 %in% paste0("Stage ", c("IA", "IB"))] = "I"
clinical$ajcc_pathologic_stage3[clinical$ajcc_pathologic_stage3 %in% paste0("Stage ", c("II", "IIA", "IIB"))] = "II ~ IV"
clinical$ajcc_pathologic_stage3[clinical$ajcc_pathologic_stage3 %in% paste0("Stage ", c("III", "IIIA", "IIIB", "IIIC"))] = "II ~ IV"
clinical$ajcc_pathologic_stage3[clinical$ajcc_pathologic_stage3 %in% paste0("Stage ", c("IV", "IVA"))] = "II ~ IV"
clinical$ajcc_pathologic_stage3[clinical$ajcc_pathologic_stage3 == "'--"] = NA
clinical$Age = rep("test", nrow(clinical))
clinical$Age[clinical$age_at_index<60] = "<60"
clinical$Age[clinical$age_at_index>=60] = ">=60"
clinical$Age = factor(clinical$Age, levels = c("<60", ">=60"))
##########
cell_type_interest <- c('myCAF', 'GC B cells')

escc <- SetIdent(escc, value="Annotation_v2")
escc$Annotation_v2 %>% table()
escc <- RenameIdents(escc, "GC DZ" = "GC B Cell", "GC LZ" = "GC B Cell")
escc$Annotation <- escc@active.ident
escc <- subset(escc, Annotation != c("."))
escc <- subset(escc, Annotation != c("Mito-high CD4"))
escc <- subset(escc, Annotation != c("Mito-high CD8"))

escc_1000 <- quickMarkers(escc@assays$RNA@counts, escc$Annotation, N=1000)
soup_markers <- escc_1000 %>% filter(qval < 0.05) %>% filter(geneFrequency >= 0.1) %>% filter(geneFrequencyOutsideCluster < 0.03)
myCAF <- soup_markers %>% filter(cluster == c("myCAF")) %>% pull(gene)

GC_B <- soup_markers %>% filter(cluster == c("GC B Cell")) %>% pull(gene)

Gene_list <- c(list(myCAF), list(GC_B))
names(Gene_list) <- c("myCAF", "GC_B")
Hallmark_MAT_ssgsea <- gsva(as.matrix(ESCA), Gene_list, method="ssgsea", mx.diff=TRUE, tau=1, verbose=TRUE)

ESCA = ESCA[,paste0(clinical$case_submitter_id,"-01")]
Hallmark_MAT_ssgsea = Hallmark_MAT_ssgsea[,paste0(clinical$case_submitter_id,"-01")]
surv_plotting = function(GSVA_results, clinical_table, markers, cutoff_option, output_path){
  ESCA_tmp = GSVA_results
  Exp_standard = seq(1, ncol(ESCA_tmp))
  ESCA_tmp = ESCA_tmp[,paste0(clinical[,2],"-01")]
  if(cutoff_option == "Q3"){
    Exp_standard[ESCA_tmp[markers,]>=quantile(as.numeric(ESCA_tmp[markers,]))[4]]="High"
    Exp_standard[ESCA_tmp[markers,]<quantile(as.numeric(ESCA_tmp[markers,]))[4]]="Low"
  }else if(cutoff_option == "Q1"){
    Exp_standard[ESCA_tmp[markers,]>=quantile(as.numeric(ESCA_tmp[markers,]))[2]]="High"
    Exp_standard[ESCA_tmp[markers,]<quantile(as.numeric(ESCA_tmp[markers,]))[2]]="Low"
  }else if(cutoff_option == "mean"){
    Exp_standard[ESCA_tmp[markers,]>=mean(as.numeric(ESCA_tmp[markers,]), na.rm = TRUE)]="High"
    Exp_standard[ESCA_tmp[markers,]<mean(as.numeric(ESCA_tmp[markers,]), na.rm = TRUE)]="Low"
  }else if(cutoff_option == "median"){
    Exp_standard[ESCA_tmp[markers,]>=median(as.numeric(ESCA_tmp[markers,]), na.rm = TRUE)]="High"
    Exp_standard[ESCA_tmp[markers,]<median(as.numeric(ESCA_tmp[markers,]), na.rm = TRUE)]="Low"
  }
  clinical_tmp = data.frame(clinical_table, Exp_standard)
  clinical_tmp$vital_status[clinical_tmp$vital_status == "Alive"] = 0
  clinical_tmp$vital_status[clinical_tmp$vital_status == "Dead"] = 1
  clinical_tmp$cont_score = as.numeric(ESCA_tmp[markers,])
  clinical_tmp$vital_status = as.numeric(clinical_tmp$vital_status)
  clinical_tmp$Exp_standard <-factor(clinical_tmp$Exp_standard)
  followup = as.numeric(apply(clinical_tmp, 1, function(x) {if(x[10]>=0){x[10]}else{x[48]}}))
  clinical_tmp = data.frame(clinical_tmp, followup)
  clinical_tmp = clinical_tmp[grepl("Low", clinical_tmp$Exp_standard) | grepl("High", clinical_tmp$Exp_standard),]
  fit<-survfit(Surv((followup/30.5),vital_status)~Exp_standard, data=clinical_tmp)
  test=survdiff(Surv((followup/30.5),vital_status)~Exp_standard,data=clinical_tmp)
  pval = 1 - pchisq(test$chisq, length(test$n) - 1)
  assign("clinical_tmp", clinical_tmp, envir = .GlobalEnv)
  NGS_OS_plot = ggsurvplot(
    fit = fit,size = 1.5,censor.size=2,xlab = "Time in Months",ylab = "OS",surv.median.line = "hv",
    palette = c("Exp_standard=Low"="#0072B5FF", "Exp_standard=High"="#BC3C29FF"),
    xlim=c(0,36),ylim=c(0,1),break.time.by=6,pval = T, pval.coord = c(1, 0.25),surv.scale = "percent", title=markers,
    risk.table = TRUE, risk.table.height = 0.3, conf.int = F,
    legend.title=paste0("Cutoff Option is ", cutoff_option))
  NGS_OS_plot$plot <- NGS_OS_plot$plot + thickness+ scale_y_continuous(breaks=seq(0, 1, by=0.2),labels=seq(0, 100, by=20),lim=c(0,1))
  pdf(file = paste0(output_path, "/", cutoff_option, "_", round(pval,3),"_", markers, "_survival_plot.pdf"), width=7, height=6)
  print(NGS_OS_plot)
  dev.off()
}

for(i_geneset in rownames(Hallmark_MAT_ssgsea)){surv_plotting(Hallmark_MAT_ssgsea, clinical, i_geneset, "median", './')}

clinical_tmp$Exp_standard <- relevel(factor(clinical_tmp$Exp_standard), ref = "Low")

# HR 
cox_bin <- coxph(Surv(followup/30.5, vital_status) ~ Exp_standard, data = clinical_tmp)
s_bin   <- summary(cox_bin)

HR  <- s_bin$coef["Exp_standardHigh", "exp(coef)"]
LCL <- s_bin$conf.int["Exp_standardHigh", "lower .95"]
UCL <- s_bin$conf.int["Exp_standardHigh", "upper .95"]

cat(sprintf("HR (High vs Low) = %.3f (95%% CI %.3f–%.3f)\n", HR, LCL, UCL))

NGS_OS_plot$plot <- NGS_OS_plot$plot +
  annotate("text", x = 1, y = 0.15,
           label = sprintf("HR=%.2f (95%% CI %.2f–%.2f)\nLog-rank p=%.3g", HR, LCL, UCL, pval),
           hjust = 0, vjust = 0)
print(NGS_OS_plot)




