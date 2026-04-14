# Fig S1B ----

library(dplyr)
library(ggplot2)
library(ggforce)
library(tibble)
library(gridExtra)
library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")  # 또는 /opt/X11/... 가능
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

date <- "20250806"
project_name <- "escc"
fileidentity <- "sunburst_combined"
myname <- "jk"

svg_filename <- paste0(date, "_", project_name, "_", fileidentity, "_", myname, ".svg")
pdf_filename <- paste0(date, "_", project_name, "_", fileidentity, "_", myname, ".pdf")

donor_color <- c(
  "19" = "darkgreen", "20" = "green2", "21" = "#26A63A", "22" = "#B4B61A",
  "24" = "#704D9E", "26" = "#CF63A6", "27" = "#F7A086", "28" = "#F3E79A",
  "31" = "#273871", "33" = "blue", "36" = "#5087C1", "37" = "#ACCCE4"
)

pathology_color <- c(
  "Normal" = "#83f52c",
  "Dysplasia" = "#f3f315",
  "Microinvasive" = "#ff6600",
  "Macroinvasive" = "#ff0099"
)

generate_sunburst <- function(data, platform_label, center_label, center_donor = NULL) {
  plot_df <- data %>%
    group_by(donor, orig.ident, pathology) %>%
    summarise(count = n(), .groups = "drop")
  
  total_count <- sum(plot_df$count)
  donor_df <- plot_df %>%
    group_by(donor) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    mutate(
      fraction = count / sum(count),
      start = cumsum(lag(fraction, default = 0)) * 2 * pi
    )
  
  rotation_angle <- 0
  if (!is.null(center_donor)) {
    donor_start_angle <- donor_df %>%
      filter(donor == center_donor) %>%
      pull(start)
    if (length(donor_start_angle) > 0) {
      rotation_angle <- -donor_start_angle
    }
  }
  
  inner_df <- donor_df %>%
    mutate(
      start = start + rotation_angle,
      end = start + fraction * 2 * pi,
      mid = (start + start + fraction * 2 * pi) / 2,
      label = paste0("Pt. ", donor, "\n", format(count, big.mark = ","))
    )
  
  outer_df <- plot_df %>%
    left_join(inner_df %>% select(donor, donor_start = start, donor_end = end), by = "donor") %>%
    group_by(donor) %>%
    mutate(
      donor_range = donor_end[1] - donor_start[1],
      sub_fraction = count / sum(count),
      sub_start = cumsum(lag(sub_fraction, default = 0)) * donor_range + donor_start,
      sub_end = sub_start + sub_fraction * donor_range,
      mid = (sub_start + sub_end) / 2,
      percent_label = paste0(orig.ident, "\n(", round(sub_fraction * 100, 1), "%)")
    ) %>%
    ungroup()
  
  outer_df$pathology <- factor(outer_df$pathology, levels = names(pathology_color))
  
  ggplot() +
    geom_arc_bar(
      data = inner_df,
      aes(x0 = 0, y0 = 0, r0 = 0.4, r = 0.7,
          start = start, end = end, fill = donor),
      color = "white"
    ) +
    geom_arc_bar(
      data = outer_df,
      aes(x0 = 0, y0 = 0, r0 = 0.7, r = 1,
          start = sub_start, end = sub_end, fill = pathology),
      color = "white"
    ) +
    geom_text(
      data = inner_df,
      aes(x = 0.55 * cos(mid), y = 0.55 * sin(mid), label = label),
      size = 2, family = "Arial"
    ) +
    geom_text(
      data = outer_df,
      aes(x = 1.05 * cos(mid), y = 1.05 * sin(mid), label = percent_label),
      size = 1.5, family = "Arial"
    ) +
    annotate("text", x = 0, y = 0, label = center_label, size = 3, fontface = "plain", family = "Arial") +
    coord_fixed() +
    theme_void(base_family = "Arial") +
    ggtitle(paste0(platform_label, " Cell Composition by Donor and Sample")) +
    theme(plot.title = element_text(hjust = 0.5, size = 6, family = "Arial")) +
    scale_fill_manual(values = c(donor_color, pathology_color), guide = "none")
}

object <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Latest_Figures/Figures_20250806/object/20250731_escc_global_final_annotation_hj.rds")

metadata <- object@meta.data[, c("platform", "donor", "orig.ident", "pathology")] %>%
  mutate(pathology = as.character(pathology)) %>%
  filter(!is.na(pathology), pathology != "")

scRNA_df <- metadata %>% filter(platform == "scRNA")
snRNA_df <- metadata %>% filter(platform == "snRNA")

p1 <- generate_sunburst(scRNA_df, "scRNA", "Cells\n45,744", center_donor = "22")
p2 <- generate_sunburst(snRNA_df, "snRNA", "Nuclei\n105,279", center_donor = "20")

svg(svg_filename, width = cm_to_inch(18), height = cm_to_inch(9), family = "Arial")
grid.arrange(p1, p2, nrow = 1)
dev.off()

saveRDS(list(p1 = p1, p2 = p2), paste0(date, "_", project_name, "_", fileidentity, "_plots.rds"))

# Fig S1C ----

# Fig S1D ----

library(dplyr)
library(Seurat)

setwd("~/Dropbox/ESCC_snrna/202507/")

escc <- readRDS("R_objects/20250731_escc_global_final_annotation_hj.rds")

pt_to_pixels <- function(pt_size, dpi = 600) {
  px <- (pt_size * dpi / 72)
  return(round(px, 2))
}
pt_to_geom_size <- function(pt) round(pt / 2.85, 2)

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
mylevel <- c("Squamous cell", "Glandular cell", "Endothelial cell", "Fibroblast/fDC/FRC", "SMC/Pericyte", "B cell", "Plasma cell", "T/NK cell", "Monocyte/DC/Macrophage", "Mast")
escc$Annotation <- factor(escc$Annotation, levels = mylevel)

mycol <- c("#E41A1C", "#F781BF", "darkgreen", "#4DAF4A", "green", "#222F75","blue", "#377EB8", "#54B0E4", "cyan")
names(mycol) <- mylevel

set.seed(1234); p1 <- DimPlot(escc, group.by = "Annotation", split.by = 'platform', cols = mycol, raster = T, pt.size = 1, label = F) & theme_void(base_family = "Arial") & theme_void(base_family = "Arial") & theme(plot.title = element_text(hjust = 0.5, size = 6), aspect.ratio = 1) & NoLegend()

cairo_pdf("~/Dropbox/escc_figure_share/FigS1/20250817_escc_FigS1D_km.pdf", width = cm_to_inch(12), height = cm_to_inch(6), family="Arial")
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

# Fig S1E ----

library(Seurat)
library(dplyr)
library(ggplot2)
library(MAST)
library(showtext)
library(stringr)

project_name <- "escc"
date <- "20250828"

cm_to_inch <- function(cm) cm / 2.54

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

setwd("~/Dropbox/project/ESCC/submit/analysis/platform_deg")

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250731_escc_global_final_annotation_hj.rds")

seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(Annotation = case_when(
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


unique(seurat_obj$Annotation)

## Plasma cell ----
celltype_tmp <- "Plasma cell"
seurat_obj_tmp <- subset(seurat_obj, subset = Annotation == celltype_tmp)
seurat_obj_tmp$platform %>% table()

set.seed(1234); seurat_obj_tmp_subset <- seurat_obj_tmp@meta.data %>% group_by(platform) %>% slice_sample(prop = 0.1, replace = F) %>% pull(cellname); seurat_obj_tmp <- subset(seurat_obj_tmp, cells = seurat_obj_tmp_subset)

mydat <- seurat_obj_tmp@assays$RNA@data[rownames(seurat_obj_tmp)[!str_detect(rownames(seurat_obj_tmp), pattern = "MT\\-")],]
sn_cells <- seurat_obj_tmp@meta.data %>% dplyr::filter(platform == "snRNA") %>% rownames()
sc_cells <- seurat_obj_tmp@meta.data %>% dplyr::filter(platform == "scRNA") %>% rownames()
mydat_tmp <- mydat[, sn_cells]
mydat_tmp <- mydat[, sc_cells]

res <- FindMarkers(seurat_obj_tmp, ident.1 = "snRNA", ident.2 = "scRNA", group.by = "platform", test.use = "MAST", latent.vars = "donor", logfc.threshold = 0, min.pct = 0, verbose = T)
saveRDS(res, paste0("res", "_plasma_cell", "_MAST_fixed_donor_hj.rds"))

## Mast cell ----
celltype_tmp <- "Mast"
seurat_obj_tmp <- subset(seurat_obj, subset = Annotation == celltype_tmp)
seurat_obj_tmp$platform %>% table()

set.seed(1234);sn_cell <-  seurat_obj_tmp@meta.data %>% dplyr::filter(platform == "snRNA") %>% dplyr::pull(cellname) %>% sample(., size = 244, replace = F)
set.seed(1234);sc_cell <-  seurat_obj_tmp@meta.data %>% dplyr::filter(platform == "scRNA") %>% dplyr::pull(cellname)
seurat_obj_tmp <- subset(seurat_obj_tmp, cells = c(sc_cell, sn_cell))

res <- FindMarkers(seurat_obj_tmp, ident.1 = "snRNA", ident.2 = "scRNA", group.by = "platform", test.use = "MAST", latent.vars = "donor", logfc.threshold = 0, min.pct = 0, verbose = T)
saveRDS(res, paste0("res", "_mastcell", "_MAST_fixed_donor_hj.rds"))

## B cell ----
celltype_tmp <- "B cell"
seurat_obj_tmp <- subset(seurat_obj, subset = Annotation == celltype_tmp)
seurat_obj_tmp$platform %>% table()

set.seed(1234); seurat_obj_tmp_subset <- seurat_obj_tmp@meta.data %>% group_by(platform) %>% slice_sample(prop = 0.1, replace = F) %>% pull(cellname); seurat_obj_tmp <- subset(seurat_obj_tmp, cells = seurat_obj_tmp_subset)

res <- FindMarkers(seurat_obj_tmp, ident.1 = "snRNA", ident.2 = "scRNA", group.by = "platform", test.use = "MAST", latent.vars = "donor", logfc.threshold = 0, min.pct = 0, verbose = T)
saveRDS(res, paste0("res", celltype_tmp, "_MAST_fixed_donor_hj.rds"))

## "Monocyte/DC/Macrophage" ----
celltype_tmp <- "Monocyte/DC/Macrophage"
seurat_obj_tmp <- subset(seurat_obj, subset = Annotation == celltype_tmp)
seurat_obj_tmp$platform %>% table()

set.seed(1234); seurat_obj_tmp_subset <- seurat_obj_tmp@meta.data %>% group_by(platform) %>% slice_sample(prop = 0.1, replace = F) %>% pull(cellname); seurat_obj_tmp <- subset(seurat_obj_tmp, cells = seurat_obj_tmp_subset)

res <- FindMarkers(seurat_obj_tmp, ident.1 = "snRNA", ident.2 = "scRNA", group.by = "platform", test.use = "MAST", latent.vars = "donor", logfc.threshold = 0, min.pct = 0, verbose = T)
saveRDS(res, paste0("res", "_monocyte_dc_macro", "_MAST_fixed_donor_hj.rds"))

## "T/NK cell" ----
celltype_tmp <- "T/NK cell"
seurat_obj_tmp <- subset(seurat_obj, subset = Annotation == celltype_tmp)
seurat_obj_tmp$platform %>% table()

set.seed(1234); seurat_obj_tmp_subset <- seurat_obj_tmp@meta.data %>% group_by(platform) %>% slice_sample(prop = 0.1, replace = F) %>% pull(cellname); seurat_obj_tmp <- subset(seurat_obj_tmp, cells = seurat_obj_tmp_subset)

res <- FindMarkers(seurat_obj_tmp, ident.1 = "snRNA", ident.2 = "scRNA", group.by = "platform", test.use = "MAST", latent.vars = "donor", logfc.threshold = 0, min.pct = 0, verbose = T)
saveRDS(res, paste0("res", "_T_NK", "_MAST_fixed_donor_hj.rds"))

## "SMC/Pericyte" ----
celltype_tmp <- "SMC/Pericyte"
seurat_obj_tmp <- subset(seurat_obj, subset = Annotation == celltype_tmp)
seurat_obj_tmp$platform %>% table()

set.seed(1234); seurat_obj_tmp_subset <- seurat_obj_tmp@meta.data %>% group_by(platform) %>% slice_sample(prop = 0.1, replace = F) %>% pull(cellname); seurat_obj_tmp <- subset(seurat_obj_tmp, cells = seurat_obj_tmp_subset)

res <- FindMarkers(seurat_obj_tmp, ident.1 = "snRNA", ident.2 = "scRNA", group.by = "platform", test.use = "MAST", latent.vars = "donor", logfc.threshold = 0, min.pct = 0, verbose = T)
saveRDS(res, paste0("res", "_SMC_Pericyte", "_MAST_fixed_donor_hj.rds"))

## "Endothelial cell" ----
celltype_tmp <- "Endothelial cell"
seurat_obj_tmp <- subset(seurat_obj, subset = Annotation == celltype_tmp)
seurat_obj_tmp$platform %>% table()

set.seed(1234); seurat_obj_tmp_subset <- seurat_obj_tmp@meta.data %>% group_by(platform) %>% slice_sample(prop = 0.1, replace = F) %>% pull(cellname); seurat_obj_tmp <- subset(seurat_obj_tmp, cells = seurat_obj_tmp_subset)

res <- FindMarkers(seurat_obj_tmp, ident.1 = "snRNA", ident.2 = "scRNA", group.by = "platform", test.use = "MAST", latent.vars = "donor", logfc.threshold = 0, min.pct = 0, verbose = T)
saveRDS(res, paste0("res", "_Endothelial_cell", "_MAST_fixed_donor_hj.rds"))

## "Fibroblast/fDC/FRC" ----
celltype_tmp <- "Fibroblast/fDC/FRC"
seurat_obj_tmp <- subset(seurat_obj, subset = Annotation == celltype_tmp)
seurat_obj_tmp$platform %>% table()

set.seed(1234); seurat_obj_tmp_subset <- seurat_obj_tmp@meta.data %>% group_by(platform) %>% slice_sample(prop = 0.1, replace = F) %>% pull(cellname); seurat_obj_tmp <- subset(seurat_obj_tmp, cells = seurat_obj_tmp_subset)

res <- FindMarkers(seurat_obj_tmp, ident.1 = "snRNA", ident.2 = "scRNA", group.by = "platform", test.use = "MAST", latent.vars = "donor", logfc.threshold = 0, min.pct = 0, verbose = T)
saveRDS(res, paste0("res", "_fibroblast_fDC_FRC", "_MAST_fixed_donor_hj.rds"))

## "Glandular cell" ----
celltype_tmp <- "Glandular cell"
seurat_obj_tmp <- subset(seurat_obj, subset = Annotation == celltype_tmp)
seurat_obj_tmp$platform %>% table()

set.seed(1234);sn_cell <-  seurat_obj_tmp@meta.data %>% dplyr::filter(platform == "snRNA") %>% dplyr::pull(cellname) %>% sample(., size = 244, replace = F)
set.seed(1234);sc_cell <-  seurat_obj_tmp@meta.data %>% dplyr::filter(platform == "scRNA") %>% dplyr::pull(cellname)
seurat_obj_tmp <- subset(seurat_obj_tmp, cells = c(sc_cell, sn_cell))

res <- FindMarkers(seurat_obj_tmp, ident.1 = "snRNA", ident.2 = "scRNA", group.by = "platform", test.use = "MAST", latent.vars = "donor", logfc.threshold = 0, min.pct = 0, verbose = T)
saveRDS(res, paste0("res", "_glandular_cell", "_MAST_fixed_donor_hj.rds"))

## squamous cell ----
celltype_tmp <- "Squamous cell"
seurat_obj_tmp <- subset(seurat_obj, subset = Annotation == celltype_tmp)
seurat_obj_tmp$platform %>% table()

set.seed(1234);sn_cell <-  seurat_obj_tmp@meta.data %>% dplyr::filter(platform == "snRNA") %>% dplyr::pull(cellname) %>% sample(., size = 1361, replace = F)
set.seed(1234);sc_cell <-  seurat_obj_tmp@meta.data %>% dplyr::filter(platform == "scRNA") %>% dplyr::pull(cellname)
seurat_obj_tmp <- subset(seurat_obj_tmp, cells = c(sc_cell, sn_cell))

res <- FindMarkers(seurat_obj_tmp, ident.1 = "snRNA", ident.2 = "scRNA", group.by = "platform", test.use = "MAST", latent.vars = "donor", logfc.threshold = 0, min.pct = 0, verbose = T)
saveRDS(res, paste0("res", "_squamous_cell", "_MAST_fixed_donor_hj.rds"))

res_final <- data.frame()

celltype_tmp <- "Plasma cell"
res <- readRDS(paste0("res", "_plasma_cell", "_MAST_fixed_donor_hj.rds"))
res %>% dplyr::filter(p_val_adj > 0.05) %>% dim()
res <- res[rownames(res)[!str_detect(rownames(res), pattern = "MT\\-")],]
res <- res %>% dplyr::mutate(DEG_category = case_when(p_val_adj < 0.05 & avg_log2FC < -1 ~ "scRNA high",
                                                      p_val_adj < 0.05 & avg_log2FC > 1 ~ "snRNA high", 
                                                      T ~ "No significant difference"))

res <- res %>% dplyr::group_by(DEG_category) %>% dplyr::summarise(n=n())
res$celltype <- celltype_tmp
res_final <- rbind(res_final, res)

celltype_tmp <- "Mast"
res <- readRDS(paste0("res", "_mastcell", "_MAST_fixed_donor_hj.rds"))
res %>% dplyr::filter(p_val_adj > 0.05) %>% dim()
res <- res[rownames(res)[!str_detect(rownames(res), pattern = "MT\\-")],]
res <- res %>% dplyr::mutate(DEG_category = case_when(p_val_adj < 0.05 & avg_log2FC < -1 ~ "scRNA high",
                                                      p_val_adj < 0.05 & avg_log2FC > 1 ~ "snRNA high", 
                                                      T ~ "No significant difference"))

res <- res %>% dplyr::group_by(DEG_category) %>% dplyr::summarise(n=n())
res$celltype <- celltype_tmp
res_final <- rbind(res_final, res)

celltype_tmp <- "B cell"
res <- readRDS(paste0("res", "_B_cell", "_MAST_fixed_donor_hj.rds"))
res %>% dplyr::filter(p_val_adj > 0.05) %>% dim()
res <- res[rownames(res)[!str_detect(rownames(res), pattern = "MT\\-")],]
res <- res %>% dplyr::mutate(DEG_category = case_when(p_val_adj < 0.05 & avg_log2FC < -1 ~ "scRNA high",
                                                      p_val_adj < 0.05 & avg_log2FC > 1 ~ "snRNA high", 
                                                      T ~ "No significant difference"))

res <- res %>% dplyr::group_by(DEG_category) %>% dplyr::summarise(n=n())
res$celltype <- celltype_tmp
res_final <- rbind(res_final, res)

celltype_tmp <- "Monocyte/DC/Macrophage"
res <- readRDS(paste0("res", "_monocyte_dc_macro", "_MAST_fixed_donor_hj.rds"))
res %>% dplyr::filter(p_val_adj > 0.05) %>% dim()
res <- res[rownames(res)[!str_detect(rownames(res), pattern = "MT\\-")],]
res <- res %>% dplyr::mutate(DEG_category = case_when(p_val_adj < 0.05 & avg_log2FC < -1 ~ "scRNA high",
                                                      p_val_adj < 0.05 & avg_log2FC > 1 ~ "snRNA high", 
                                                      T ~ "No significant difference"))

res <- res %>% dplyr::group_by(DEG_category) %>% dplyr::summarise(n=n())
res$celltype <- celltype_tmp
res_final <- rbind(res_final, res)

celltype_tmp <- "T/NK cell"
res <- readRDS(paste0("res", "_T_NK", "_MAST_fixed_donor_hj.rds"))
res %>% dplyr::filter(p_val_adj > 0.05) %>% dim()
res <- res[rownames(res)[!str_detect(rownames(res), pattern = "MT\\-")],]
res <- res %>% dplyr::mutate(DEG_category = case_when(p_val_adj < 0.05 & avg_log2FC < -1 ~ "scRNA high",
                                                      p_val_adj < 0.05 & avg_log2FC > 1 ~ "snRNA high", 
                                                      T ~ "No significant difference"))

res <- res %>% dplyr::group_by(DEG_category) %>% dplyr::summarise(n=n())
res$celltype <- celltype_tmp
res_final <- rbind(res_final, res)

celltype_tmp <- "SMC/Pericyte"
res <- readRDS(paste0("res", "_SMC_Pericyte", "_MAST_fixed_donor_hj.rds"))
res %>% dplyr::filter(p_val_adj > 0.05) %>% dim()
res <- res[rownames(res)[!str_detect(rownames(res), pattern = "MT\\-")],]
res <- res %>% dplyr::mutate(DEG_category = case_when(p_val_adj < 0.05 & avg_log2FC < -1 ~ "scRNA high",
                                                      p_val_adj < 0.05 & avg_log2FC > 1 ~ "snRNA high", 
                                                      T ~ "No significant difference"))

res <- res %>% dplyr::group_by(DEG_category) %>% dplyr::summarise(n=n())
res$celltype <- celltype_tmp
res_final <- rbind(res_final, res)

celltype_tmp <- "Endothelial cell"
res <- readRDS("res_Endothelial_cell_MAST_fixed_donor_hj.rds")
res %>% dplyr::filter(p_val_adj > 0.05) %>% dim()
res <- res[rownames(res)[!str_detect(rownames(res), pattern = "MT\\-")],]
res <- res %>% dplyr::mutate(DEG_category = case_when(p_val_adj < 0.05 & avg_log2FC < -1 ~ "scRNA high",
                                                      p_val_adj < 0.05 & avg_log2FC > 1 ~ "snRNA high", 
                                                      T ~ "No significant difference"))

res <- res %>% dplyr::group_by(DEG_category) %>% dplyr::summarise(n=n())
res$celltype <- celltype_tmp
res_final <- rbind(res_final, res)

celltype_tmp <- "Fibroblast/fDC/FRC"
res <- readRDS("res_fibroblast_fDC_FRC_MAST_fixed_donor_hj.rds")
res %>% dplyr::filter(p_val_adj > 0.05) %>% dim()
res <- res[rownames(res)[!str_detect(rownames(res), pattern = "MT\\-")],]
res <- res %>% dplyr::mutate(DEG_category = case_when(p_val_adj < 0.05 & avg_log2FC < -1 ~ "scRNA high",
                                                      p_val_adj < 0.05 & avg_log2FC > 1 ~ "snRNA high", 
                                                      T ~ "No significant difference"))

res <- res %>% dplyr::group_by(DEG_category) %>% dplyr::summarise(n=n())
res$celltype <- celltype_tmp
res_final <- rbind(res_final, res)

celltype_tmp <- "Glandular cell"
res <- readRDS("res_glandular_cell_MAST_fixed_donor_hj.rds")
res %>% dplyr::filter(p_val_adj > 0.05) %>% dim()
res <- res[rownames(res)[!str_detect(rownames(res), pattern = "MT\\-")],]
res <- res %>% dplyr::mutate(DEG_category = case_when(p_val_adj < 0.05 & avg_log2FC < -1 ~ "scRNA high",
                                                      p_val_adj < 0.05 & avg_log2FC > 1 ~ "snRNA high", 
                                                      T ~ "No significant difference"))

res <- res %>% dplyr::group_by(DEG_category) %>% dplyr::summarise(n=n())
res$celltype <- celltype_tmp
res_final <- rbind(res_final, res)

celltype_tmp <- "Squamous cell"
res <- readRDS("./res_squamous_cell_MAST_fixed_donor_hj.rds")
res %>% dplyr::filter(p_val_adj > 0.05) %>% dim()
res <- res[rownames(res)[!str_detect(rownames(res), pattern = "MT\\-")],]
res <- res %>% dplyr::mutate(DEG_category = case_when(p_val_adj < 0.05 & avg_log2FC < -1 ~ "scRNA high",
                                                      p_val_adj < 0.05 & avg_log2FC > 1 ~ "snRNA high", 
                                                      T ~ "No significant difference"))

res <- res %>% dplyr::group_by(DEG_category) %>% dplyr::summarise(n=n())
res$celltype <- celltype_tmp
res_final <- rbind(res_final, res)

res_final <- res_final %>% dplyr::group_by(celltype) %>% dplyr::mutate(total = sum(n)) %>% dplyr::mutate(prop = n/total)

res_final$DEG_category <- factor(res_final$DEG_category, levels = c("scRNA high", "snRNA high", "No significant difference"))


# final figure
res_final %>% ggplot(aes(x = celltype, y = prop, fill = DEG_category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("deepskyblue", "coral2", "navy")) +
  theme_classic(base_family = "Arial") +
  theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(family = "Arial"))

# Fig 3F ----

setwd("~/Dropbox/project/ESCC/submit/analysis/platform_deg")

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250731_escc_global_final_annotation_hj.rds")

seurat_obj_tmp <- subset(seurat_obj, subset = Annotation == "B cell")
seurat_obj_tmp1 <- subset(seurat_obj_tmp, subset = platform == "snRNA")
seurat_obj_tmp2 <- subset(seurat_obj_tmp, subset = platform == "scRNA")
mygenes <- rownames(seurat_obj_tmp1@assays$RNA$data)
df <- data.frame()
for (gene_tmp in mygenes) {
  gene_mean1 <- mean(seurat_obj_tmp1@assays$RNA$data[gene_tmp, ])
  gene_mean2 <- mean(seurat_obj_tmp2@assays$RNA$data[gene_tmp, ])
  df_tmp <- data.frame(gene = gene_tmp, snRNA_mean = gene_mean1, scRNA_mean = gene_mean2)
  df <- rbind(df, df_tmp)
}
saveRDS(df, "cor_tmp1_hj.rds")

mygenes <- mygenes[11185:length(mygenes)]
for (gene_tmp in mygenes) {
  gene_mean1 <- mean(seurat_obj_tmp1@assays$RNA$data[gene_tmp, ])
  gene_mean2 <- mean(seurat_obj_tmp2@assays$RNA$data[gene_tmp, ])
  df_tmp <- data.frame(gene = gene_tmp, snRNA_mean = gene_mean1, scRNA_mean = gene_mean2)
  df <- rbind(df, df_tmp)
}

saveRDS(df, "cor_hj.rds")

head(df)

df <- readRDS("cor_hj.rds")

df <- df[!str_detect(df$gene, pattern = "RPS|RPL|MT\\-"),]

cor.test(df$snRNA_mean, df$scRNA_mean) # p-value < 2.2e-16, cor = 0.446726 

df %>% ggplot(aes(x = snRNA_mean, y = scRNA_mean)) +
  geom_point(size = 1) +
  xlab("scRNA") +
  ylab("snRNA") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_classic(base_family = "Arial") +
  theme(aspect.ratio = 1)

# Fig S1H ----

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
size_pt <- 7
cm_to_inch <- function(cm) cm / 2.54

date = "20250810"
project_name = "escc"

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig1")

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


library(colorspace)
#hcl_palettes(plot = TRUE)

mylevel <- c("Squamous cell", "Glandular cell", "Endothelial cell", "Fibroblast/fDC/FRC", "SMC/Pericyte", "B cell", "Plasma cell", "T/NK cell", "Monocyte/DC/Macrophage", "Mast")
escc$Annotation <- factor(escc$Annotation, levels = mylevel)

df <- escc@meta.data %>% dplyr::select(Annotation, platform, donor)
df$donor <- as.factor(df$donor)
df$Annotation <- factor(df$Annotation, levels = rev(mylevel))
df <- df %>% dplyr::group_by(Annotation, donor) %>% dplyr::summarise(n=n())
df <- df %>% dplyr::group_by(Annotation) %>% dplyr::mutate(total = sum(n)) %>% dplyr::mutate(prop = n/total)

p1 <- df %>% ggplot(aes(x = Annotation, y = prop,  fill = donor)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", size = 0.1) +
  xlab("Cell type") +
  ylab("Proportion of cells") +
  coord_flip() +
  scale_fill_manual(values = c("darkgreen", "green2", "#26A63A", "#B4B61A", "#704D9E", "#CF63A6", "#F7A086", "#F3E79A", "#273871", "blue", "#5087C1", "#ACCCE4")) +
  theme_set(theme_classic(base_family = "Arial")) +
  theme(line = element_line(linewidth = 0.3),
        axis.text = element_text(size = size_pt-1, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = size_pt-2),
        axis.title = element_text(size = size_pt-1),
        legend.text = element_text(size = size_pt-3),
        legend.position = "none",
        legend.title = element_text(size = size_pt-3),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width  = unit(0.1, "cm"),
        legend.spacing.x  = unit(0.02, "cm"),
        plot.margin = unit(c(0, 0.1, 0, 0), "cm")) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

df <- escc@meta.data %>% dplyr::select(Annotation, platform, donor, pathology)
df$Annotation <- factor(df$Annotation, levels = rev(mylevel))
df$pathology <- factor(df$pathology, levels = rev(c("Normal", "Dysplasia", "Microinvasive", "Macroinvasive")))
df <- df %>% dplyr::group_by(Annotation, pathology) %>% dplyr::summarise(n=n())
df <- df %>% dplyr::group_by(Annotation) %>% dplyr::mutate(total = sum(n)) %>% dplyr::mutate(prop = n/total)

p2 <- df %>% ggplot(aes(x = Annotation, y = prop,  fill = pathology)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", size = 0.1) +
  xlab("Cell type") +
  ylab("Proportion of cells") +
  coord_flip() +
  scale_fill_manual(values = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive" = "#ff6600", "Macroinvasive" = "#ff0099"),
                    breaks = c("Normal", "Dysplasia", "Microinvasive", "Macroinvasive")) +
  theme_set(theme_classic(base_family = "Arial")) +
  theme(line = element_line(linewidth = 0.3),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = size_pt-2),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = size_pt-1, colour = "black"),
        axis.title = element_text(size = size_pt-1),
        legend.text = element_text(size = size_pt-3),
        legend.position = "none",
        legend.title = element_text(size = size_pt-3),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width  = unit(0.1, "cm"),
        legend.spacing.x  = unit(0.02, "cm"),
        plot.margin = unit(c(0, 0.1, 0, 0), "cm")) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))


df <- escc@meta.data %>% dplyr::select(Annotation, platform, donor, pathology)
df$Annotation <- factor(df$Annotation, levels = rev(mylevel))
df$pathology <- factor(df$pathology, levels = rev(c("Normal", "Dysplasia", "Microinvasive", "Macroinvasive")))
df <- df %>% dplyr::group_by(Annotation, pathology) %>% dplyr::summarise(n=n())
#df <- df %>% dplyr::group_by(Annotation) %>% dplyr::mutate(total = sum(n)) %>% dplyr::mutate(prop = n/total)

p3 <- df %>% ggplot(aes(x = Annotation, y = n,  fill = pathology)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black",  size = 0.1) +
  xlab("Cell type") +
  ylab("Number of cells") +
  coord_flip() +
  scale_fill_manual(values = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive" = "#ff6600", "Macroinvasive" = "#ff0099"),
                    breaks = c("Normal", "Dysplasia", "Microinvasive", "Macroinvasive")) +
  theme_set(theme_classic(base_family = "Arial")) +
  theme(line = element_line(linewidth = 0.3),
        #axis.text.y = element_blank(),
        #axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = size_pt-2),
        #axis.ticks.y = element_blank(),
        axis.text = element_text(size = size_pt-1, colour = "black"),
        axis.title = element_text(size = size_pt-1),
        legend.text = element_text(size = size_pt-3),
        legend.position = "none",
        legend.title = element_text(size = size_pt-3),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width  = unit(0.1, "cm"),
        legend.spacing.x  = unit(0.02, "cm"),
        plot.margin = unit(c(0, 0.1, 0, 0), "cm")) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

df <- escc@meta.data %>% dplyr::select(Annotation, platform, donor, pathology)
df$donor <- as.factor(df$donor)
df$donor <- factor(df$donor, levels = rev(levels(df$donor)))
df$Annotation <- factor(df$Annotation, levels = mylevel)
df$pathology <- factor(df$pathology, levels = c("Normal", "Dysplasia", "Microinvasive", "Macroinvasive"))
df <- df %>% dplyr::group_by(donor, Annotation) %>% dplyr::summarise(n=n())
df <- df %>% dplyr::group_by(donor) %>% dplyr::mutate(total = sum(n)) %>% dplyr::mutate(prop = n/total)
df$Annotation_bar <- factor(df$Annotation, levels = rev(mylevel))

mycol <- c("#E41A1C", "#F781BF", "darkgreen", "#4DAF4A", "green", "#222F75","blue", "#377EB8", "#54B0E4", "cyan")
names(mycol) <- mylevel

p4 <- df %>% ggplot(aes(x = donor, y = prop,  fill = Annotation_bar)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", size = 0.1) +
  xlab("Patient") +
  ylab("Proportion of cells") +
  coord_flip() +
  scale_fill_manual(values = mycol, breaks = levels(df$Annotation)) +
  theme_set(theme_classic(base_family = "Arial")) +
  theme(line = element_line(linewidth = 0.3),
        axis.text.x = element_text(angle = 45, hjust = 1, size = size_pt-2),
        axis.text = element_text(size = size_pt-1, colour = "black"),
        axis.title = element_text(size = size_pt-1),
        legend.text = element_text(size = size_pt-3),
        legend.position = "none",
        legend.title = element_text(size = size_pt-3),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width  = unit(0.1, "cm"),
        legend.spacing.x  = unit(0.02, "cm"),
        plot.margin = unit(c(0, 0.1, 0, 0), "cm")) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

df <- escc@meta.data %>% dplyr::select(Annotation, platform, donor, pathology)
df$donor <- as.factor(df$donor)
df$donor <- factor(df$donor, levels = rev(levels(df$donor)))
df$Annotation <- factor(df$Annotation, levels = mylevel)
#df$pathology <- factor(df$pathology, levels = c("Normal", "Dysplasia", "Microinvasive", "Macroinvasive"))
df <- df %>% dplyr::group_by(donor, Annotation) %>% dplyr::summarise(n=n())
df$Annotation_bar <- factor(df$Annotation, levels = rev(mylevel))

p5 <- df %>% ggplot(aes(x = donor, y = n,  fill = Annotation_bar)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", size = 0.1) +
  xlab("Patient") +
  ylab("Number of cells") +
  coord_flip() +
  scale_fill_manual(values = mycol, breaks = levels(df$Annotation)) +
  theme_set(theme_classic(base_family = "Arial")) +
  theme(line = element_line(linewidth = 0.3),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = size_pt-2),
        #axis.title.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.text = element_text(size = size_pt-1, colour = "black"),
        axis.title = element_text(size = size_pt-1),
        legend.text = element_text(size = size_pt-3),
        legend.position = "none",
        legend.title = element_text(size = size_pt-3),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width  = unit(0.1, "cm"),
        legend.spacing.x  = unit(0.02, "cm"),
        plot.margin = unit(c(0, 0.1, 0, 0), "cm")) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

combined_plot <- p1 + p2 + p4 + plot_layout(ncol = 3)

print(combined_plot)
dev.off()

combined_plot <- p3 + p5 + plot_layout(ncol = 2)

print(combined_plot)

# Fig S1I ----
library(Seurat)
library(future)
library(ggplot2)
library(CellChat)
visium <- readRDS('~/Dropbox/escc_figure_share/dl/20250812_escc_visium_C2F_latest_dl.rds')

visium@active.assay <- 'Spatial'

visium$niche %>% table

visium$niche <- factor(visium$niche, levels = paste0('niche',1:12))

library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54
mycol <- scPalette(12)
names(mycol) <- paste0('niche',1:12)

#width = 38, height = 90

cairo_pdf("~/Dropbox/escc_figure_share/FigS1/20250828/20250828_escc_FigS1f_Bayesspace_km.pdf", width = 7, height = 90, family="Arial")
SpatialDimPlot(visium, group.by = 'niche', ncol = 1, cols = mycol) & theme(text = element_text(family = "Arial")) & NoLegend()
dev.off()

# Fig S1J ----

library(CellChat)
library(dplyr)
library(Seurat)
library(tidyverse)

merged <- readRDS("~/Dropbox/project/ESCC/submit/figures/escc_figure_share/dl/20250816_escc_visium_add_TLS_feature_dl.rds")

c2l_data <- GetAssayData(merged, assay = "C2L", slot='data')
c2l_data.t <- t(as.matrix(c2l_data))

cluster_info <- merged$niche

cluster_means <- c2l_data.t %>%
  as.data.frame() %>%
  mutate(cluster = cluster_info) %>%
  group_by(cluster) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE))

cluster_means_long <- cluster_means %>%
  pivot_longer(cols = -cluster, names_to = "CellType", values_to = "Expression")
cluster_means_long$cluster <- factor(x = cluster_means_long$cluster, levels = c("niche8", "niche1", "niche5", "niche4", "niche12", "niche7", "niche10", "niche11", "niche3", "niche9", "niche2", "niche6"))
c2l_data <- cluster_means_long

global <- c2l_data %>%
  as.data.frame() %>%
  mutate(
    Cell_Type = case_when(
      CellType %in% c("BEC", "LEC") ~ "Endothelial",
      CellType %in% c("Basal", "Suprabasal", 'Tumor')~ "Epithelial",
      CellType %in% c('APOD high CAF', 'NAF', 'NMF', 'Pericyte', 'SMC', 'CFD high CAF', 'fDC', 'iCAF', 'myCAF', 'Progenitor iCAF') ~ "Stromal",
      TRUE ~ 'Immune'
    )
  )
df <- aggregate(Expression ~ cluster + Cell_Type, data = global, sum, na.rm = TRUE)


keep_ct <- c('Epithelial', 'Immune', 'Stromal', 'Endothelial')

df <- df %>%
  dplyr::filter(Cell_Type %in% keep_ct) %>%
  dplyr::mutate(Cell_Type = factor(Cell_Type, levels = keep_ct))

cell_colors <- setNames(scPalette(length(keep_ct)), keep_ct)

cairo_pdf('/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/FigS1/celltype_barplot/20250817_escc_global_niche_barplot_hj.pdf',width = cm_to_inch(9), height = cm_to_inch(10), family = "Arial")

ggplot(df, aes(x = cluster, y = Expression, fill = Cell_Type)) +
  geom_bar(stat = "identity", position='fill') +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(x = "Niche", y = "Celltype Proportion", fill = "Cell_Type") +
  theme_classic(base_size = 5, base_family = "Arial") +
  theme(
    aspect.ratio = 0.7,
    line = element_line(linewidth = 0.3),
    plot.title = element_text(size = 7, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 7, face = "bold"),
    axis.title.y = element_text(size = 7, face = "bold"),
    legend.position = "bottom",
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size=6),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin  = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  ) +
  guides(fill = guide_legend(nrow = 4)) +
  scale_fill_manual(values = cell_colors,
                    limits = keep_ct,
                    breaks = keep_ct,
                    drop = FALSE)
dev.off()


c2l_data <- GetAssayData(merged, assay = "C2L", slot='data')
c2l_data.t <- t(as.matrix(c2l_data))

cluster_info <- merged$niche
cluster_info <- merged$niche

cluster_means <- c2l_data.t %>%
  as.data.frame() %>%
  mutate(cluster = cluster_info) %>%
  group_by(cluster) %>%
  dplyr::summarise(across(everything(), mean, na.rm = TRUE))

cluster_means_long <- cluster_means %>%
  pivot_longer(cols = -cluster, names_to = "CellType", values_to = "Expression")
cluster_means_long$cluster <- factor(x = cluster_means_long$cluster, levels = c("niche8", "niche1", "niche5", "niche4", "niche12", "niche7", "niche10", "niche11", "niche3", "niche9", "niche2", "niche6"))
top3_with_others <- cluster_means_long %>% 
  group_by(cluster) %>% 
  arrange(desc(Expression), .by_group = TRUE) %>% 
  mutate(rank = row_number()) %>% 
  mutate(CellType = ifelse(rank > 3, "Others", as.character(CellType))) %>% 
  group_by(cluster, CellType) %>% 
  summarise(Expression = sum(Expression), .groups = "drop") %>% 
  ungroup() %>%  
  group_by(cluster) %>% 
  mutate(Percent = Expression / sum(Expression) * 100)

keep_ct <- c(
  "Tumor","Basal","CD4 Tn","CD8 CXCR6+ Trm","CD8 Inflammatory Trm",
  "GC DZ","GC LZ","Glandular cell","Memory B","NMF","Pericyte",
  "BEC","Plasma cell","SMC","Suprabasal","Others"
)

top3_with_others <- top3_with_others %>%
  dplyr::filter(CellType %in% keep_ct) %>%
  dplyr::mutate(CellType = factor(CellType, levels = keep_ct))

cell_colors <- setNames(c(scPalette(length(keep_ct))[-length(keep_ct)], "darkgrey"), keep_ct)

cairo_pdf('/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/FigS1/celltype_barplot/20250817_escc_intricate_celltype_niche_barplot_hj.pdf',width = cm_to_inch(9), height = cm_to_inch(10), family = "Arial")

ggplot(top3_with_others, aes(x = cluster, y = Expression, fill = CellType)) +
  geom_bar(stat = "identity", position='fill') +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(x = "Niche", y = "Celltype Proportion", fill = "Cell Type") +
  theme_classic(base_size = 5, base_family = "Arial") +
  theme(
    aspect.ratio = 0.7,
    line = element_line(linewidth = 0.3),
    plot.title = element_text(size=7, hjust=0.5, face="bold"),
    axis.text.x = element_text(size=6, angle=45, hjust=1),
    axis.text.y = element_text(size=6),
    axis.title.x = element_text(size = 7, face = "bold"),
    axis.title.y = element_text(size = 7, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size=6),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin  = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  ) +
  scale_fill_manual(values = cell_colors,
                    limits = keep_ct,
                    breaks = keep_ct,
                    drop = FALSE)
dev.off()

# Fig S1K ----

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/celltype_heterogenity/")
library(ggalluvial)
library(Seurat)
library(dplyr)
library(cowplot)

##myeloid cells
myeloid <- readRDS('../../object_share/20250722_escc_myeloid_flt_km.rds')
myeloid$annotation_myeloid <- gsub('Monocytes','monocytes',myeloid$annotation_myeloid)

Idents(myeloid) <- 'pathology'
table_myeloid <- myeloid@meta.data %>% select(pathology, annotation_myeloid) %>% table()
table_myeloid_v2 <- ( table_myeloid / table_myeloid %>% apply(1,sum) ) %>% data.frame()
table_myeloid_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',table_myeloid_v2$pathology)
table_myeloid_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',table_myeloid_v2$pathology)
table_myeloid_v2$pathology <- factor(table_myeloid_v2$pathology,  levels = c('Normal', 'Dysplasia', 'Microinvasive carcinoma', 'Macroinvasive carcinoma'))
table_myeloid_v2$annotation_myeloid %>% table
rm(myeloid)

##fibroblasts
fib <- readRDS("~/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250728_escc_Fib_fDC_FRC_hj_km.rds")

Idents(fib) <- 'pathology'
table_fib <- fib@meta.data %>% select(pathology, Annotation_v2) %>% table()
table_fib_v2 <- ( table_fib / table_fib %>% apply(1,sum) ) %>% data.frame()
table_fib_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',table_fib_v2$pathology)
table_fib_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',table_fib_v2$pathology)
table_fib_v2$pathology <- factor(table_fib_v2$pathology,  levels = c('Normal', 'Dysplasia', 'Microinvasive carcinoma', 'Macroinvasive carcinoma'))
table_fib_v2$Annotation_v2 %>% table
rm(fib)


##CD8T
cd8t <- readRDS("../../object_share/20250728_escc_cd8_annotation_dl.rds")

Idents(cd8t) <- 'pathology'
table_cd8t <- cd8t@meta.data %>% select(pathology, Annotation_v2) %>% table()
table_cd8t_v2 <- ( table_cd8t / table_cd8t %>% apply(1,sum) ) %>% data.frame()
table_cd8t_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',table_cd8t_v2$pathology)
table_cd8t_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',table_cd8t_v2$pathology)
table_cd8t_v2$pathology <- factor(table_cd8t_v2$pathology,  levels = c('Normal', 'Dysplasia', 'Microinvasive carcinoma', 'Macroinvasive carcinoma'))
table_cd8t_v2$Annotation_v2 %>% table
rm(cd8t)


##CD4T
cd4t <- readRDS("../../object_share/20250729_escc_cd4_celltype_annotated_dl.rds")

Idents(cd4t) <- 'pathology'
table_cd4t <- cd4t@meta.data %>% select(pathology, Annotation_v2) %>% table()
table_cd4t_v2 <- ( table_cd4t / table_cd4t %>% apply(1,sum) ) %>% data.frame()
table_cd4t_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',table_cd4t_v2$pathology)
table_cd4t_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',table_cd4t_v2$pathology)
table_cd4t_v2$pathology <- factor(table_cd4t_v2$pathology,  levels = c('Normal', 'Dysplasia', 'Microinvasive carcinoma', 'Macroinvasive carcinoma'))
table_cd4t_v2$Annotation_v2 %>% table
rm(cd4t)


##B cells
b_cells <- readRDS("../../object_share/20250723_escc_B_GCB_Plasma_flt_km_v2.rds")
b_cells$annotation_B <- gsub('naive','Naive',b_cells$annotation_B)

Idents(b_cells) <- 'pathology'
table_B <- b_cells@meta.data %>% select(pathology, annotation_B) %>% table()
table_B_v2 <- ( table_B / table_B %>% apply(1,sum) ) %>% data.frame()
table_B_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',table_B_v2$pathology)
table_B_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',table_B_v2$pathology)
table_B_v2$pathology <- factor(table_B_v2$pathology,  levels = c('Normal', 'Dysplasia', 'Microinvasive carcinoma', 'Macroinvasive carcinoma'))
table_B_v2$annotation_B %>% table
rm(b_cells)

colnames(table_B_v2) <- c("pathology", "annotation", "Freq")
colnames(table_cd4t_v2) <- c("pathology", "annotation", "Freq")
colnames(table_cd8t_v2) <- c("pathology", "annotation", "Freq")
colnames(table_fib_v2) <- c("pathology", "annotation", "Freq")
colnames(table_myeloid_v2) <- c("pathology", "annotation", "Freq")

mytable <- rbind(table_B_v2, table_cd4t_v2, table_cd8t_v2, table_fib_v2, table_myeloid_v2)

mytable %>% dplyr::group_by(pathology) %>% dplyr::summarise(var = var(Freq)) %>% ggplot(aes(x = pathology, y = var, fill = pathology)) +
  geom_bar(stat = "identity") +
  theme_classic(base_family = "Arial") +
  labs(y = "Variance", x = "Pathology") +
  theme(aspect.ratio = 0.7, axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_manual(values = c("Normal"= "#83f52c", "Dysplasia" = "#f3f315", "Microinvasive carcinoma" = "#ff6600", "Macroinvasive carcinoma" = "#ff0099"))
