# Fig 4A ----
library(dplyr)
library(Seurat)
library(CellChat)
library(ggplot2)
library(showtext)
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
# fileidentity <- "normal_dysplasia_cellchat_change"
# cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 9, height = 9, family="Arial")
# change_p1
# dev.off()

# Fig 4A ----
mat <- change_p1@matrix

groupSize <- mat[, colnames(mat) == "Tumor"]
groupSize[groupSize < 0 | is.na(groupSize) | is.infinite(groupSize)] <- 0
mat[, colnames(mat) != "Tumor"] <- 0
#mat[, colnames(mat) == "Tumor"]
mat[is.na(mat)] <- 0

fileidentity <- "cellchat_tumor_receiver_normal_to_dys"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(20), height = cm_to_inch(20), family = "Arial")
par(family = "Arial", mar = c(0, 0, 0, 0))
netVisual_circle(mat, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, vertex.label.cex = 0.2) 
dev.off()

# Fig 4B ----

library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100))


df.net_dys <- subsetCommunication(object.list$Dysplasia, slot.name = "net", 
                                  sources.use = "myCAF", targets.use = "Tumor")
df.net_norm <- subsetCommunication(object.list$Normal, slot.name = "net", 
                                   sources.use = "myCAF", targets.use = "Tumor")

df.net_dys$stage <- "Dysplasia"
df.net_norm$stage <- "Normal"

df.net_dys = df.net_dys %>% dplyr::select(ligand, prob, receptor,stage,interaction_name_2)
min_prob = min(df.net_dys$prob)
max_prob = max(df.net_dys$prob)
df.net_dys = df.net_dys %>% dplyr::mutate(prob_norm = (prob - min_prob)/(max_prob - min_prob))

df.net_norm = df.net_norm %>% dplyr::select(ligand, prob, receptor,stage,interaction_name_2)
min_prob = min(df.net_norm$prob)
max_prob = max(df.net_norm$prob)
df.net_norm = df.net_norm %>% dplyr::mutate(prob_norm = (prob - min_prob)/(max_prob - min_prob))

mydf = rbind(df.net_dys, df.net_norm)

common_cci = mydf %>% dplyr::group_by(interaction_name_2) %>% dplyr::summarise(n=n()) %>% dplyr::filter(n == 2) %>% dplyr::pull(interaction_name_2)
private_cci = mydf %>% dplyr::group_by(interaction_name_2) %>% dplyr::summarise(n=n()) %>% dplyr::filter(n == 1) %>% dplyr::pull(interaction_name_2)

mydf = mydf %>% dplyr::mutate(Pair_info = case_when(
  interaction_name_2 %in% common_cci ~ "Common pair",
  interaction_name_2 %in% private_cci ~ "Private pair"))
mydf$stage <- factor(mydf$stage, levels = c("Normal", "Dysplasia"))

sc <- scale_fill_gradientn(colours = myPalette(100))

ligand_level <- unique(sort(mydf$ligand))
mydf$ligand <- factor(mydf$ligand, levels = rev(ligand_level))

fileidentity <- "normal_dysplasia_cellchat_dotplot"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(9), height = cm_to_inch(9), family="Arial")
mydf %>%
  ggplot(aes(y = ligand, x = receptor, fill = prob_norm), color = "black") +
  geom_point(size = 1, aes(shape = Pair_info), stroke = 0.2) +
  facet_wrap(~ stage, ncol = 2) + 
  labs(y = "Ligand", x = "Receptor") +
  scale_shape_manual(values = c(21,24)) + 
  theme_bw(base_family = "Arial")+
  scale_y_discrete(position = "right") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(colour = "black", size = 5),
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


# Fig 4C ----
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(showtext)
library(ggsignif)

col_myCAF <- "#00CC00"
col_iCAF  <- "#0000FF"

# load data
visium <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/escc_figure_share/object/20250806_escc_misty_dl.rds")

EGS_19_N2        <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Analysis/Visium/Cottrazm/Object/EGS_19_N2_BoundaryDefine.rds")
EGS_19_d2_hgd    <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Analysis/Visium/Cottrazm/output/EGS_19_D2_HGD_v2/EGS-19-D2_HGD_BoundaryDefine.rds")
EGS_19_t2_hgd    <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Analysis/Visium/Cottrazm/output/EGS-19-T2_HGD_v2/EGS-19-T2_HGD_BoundaryDefine.rds")
EGS_19_t2_mic    <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Analysis/Visium/Cottrazm/output/EGS-19-T2-MIC_3p_v2/EGS-19-T2-MIC_3p_BoundaryDefine.rds")
EGS_31_st_1      <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Analysis/Visium/Cottrazm/output/EGS-31-ST-1_v2/EGS-31-ST-1_BoundaryDefine.rds")
EGS_19_d2_mic    <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Analysis/Visium/Cottrazm/output/EGS-19-D2-5p_v2/EGS-19-D2-5p_BoundaryDefine.rds")
EGS_26_d3        <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Analysis/Visium/Cottrazm/output/EGS-26-D3_v2/EGS-26-D3_BoundaryDefine.rds")
EGS_33_st_n      <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Analysis/Visium/Cottrazm/output/EGS-33-ST-N_3p_v2/EGS-33-ST-N_3p_BoundaryDefine.rds")
EGS_33_st_m      <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Analysis/Visium/Cottrazm/output/EGS-33-ST-M_3p_v2/EGS-33-ST-M_3p_BoundaryDefine.rds")
EGS_19_t2_mac    <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Analysis/Visium/Cottrazm/output/EGS_19_T2_MAC_v2/EGS_19_T2_MAC_BoundaryDefine.rds")

# load metadata
read_meta <- function(path) {
  read.table(path, sep = ",", header = TRUE) %>%
    setNames(c("Barcode", "Location2")) %>%
    tibble::column_to_rownames("Barcode")
}

EGS_19_D2_HGD_meta <- read_meta(".../EGS_19_D2_HGD_Location.csv")
EGS_19_D2_5p_meta  <- read_meta(".../EGS_19_D2_5p_Location.csv")
EGS_19_t2_hgd_meta <- read_meta(".../EGS_19_T_2_HGD_Location.csv")
EGS_19_t2_mic_meta <- read_meta(".../EGS_19_T_2_MIC_Location.csv")
EGS_26_d_3_meta    <- read_meta(".../EGS_26_D_3_Location.csv")
EGS_31_st_1_meta   <- read_meta(".../EGS_31_ST_1_Location.csv")
EGS_33_st_m_meta   <- read_meta(".../EGS_33_ST_MIC_Location.csv")
EGS_33_st_n_meta   <- read_meta(".../EGS_33_ST_N_Location.csv")

# add metadata and set idents
EGS_19_d2_hgd <- AddMetaData(EGS_19_d2_hgd, EGS_19_D2_HGD_meta) %>% SetIdent(value = "Location2")
EGS_19_d2_mic <- AddMetaData(EGS_19_d2_mic, EGS_19_D2_5p_meta)  %>% SetIdent(value = "Location2")
EGS_19_t2_hgd <- AddMetaData(EGS_19_t2_hgd, EGS_19_t2_hgd_meta) %>% SetIdent(value = "Location2")
EGS_19_t2_mic <- AddMetaData(EGS_19_t2_mic, EGS_19_t2_mic_meta) %>% SetIdent(value = "Location2")
EGS_26_d3     <- AddMetaData(EGS_26_d3,     EGS_26_d_3_meta)    %>% SetIdent(value = "Location2")
EGS_31_st_1   <- AddMetaData(EGS_31_st_1,   EGS_31_st_1_meta)   %>% SetIdent(value = "Location2")
EGS_33_st_n   <- AddMetaData(EGS_33_st_n,   EGS_33_st_n_meta)   %>% SetIdent(value = "Location2")
EGS_33_st_m   <- AddMetaData(EGS_33_st_m,   EGS_33_st_m_meta)   %>% SetIdent(value = "Location2")

# integrate Cottrazm metadata
make_df <- function(obj, prefix, col) {
  data.frame(
    Barcode    = paste0(prefix, colnames(obj)),
    Annotation = obj[[col]][, 1]
  )
}

combined_df <- rbind(
  make_df(EGS_19_N2,     "s1_",  "Location"),
  make_df(EGS_19_d2_hgd, "s2_",  "Location2"),
  make_df(EGS_19_t2_hgd, "s3_",  "Location2"),
  make_df(EGS_19_d2_mic, "s4_",  "Location2"),
  make_df(EGS_19_t2_mic, "s5_",  "Location2"),
  make_df(EGS_26_d3,     "s6_",  "Location2"),
  make_df(EGS_31_st_1,   "s7_",  "Location2"),
  make_df(EGS_33_st_m,   "s8_",  "Location2"),
  make_df(EGS_33_st_n,   "s9_",  "Location2"),
  make_df(EGS_19_t2_mac, "s10_", "Location")
) %>%
  tibble::column_to_rownames("Barcode") %>%
  setNames("Cottrazm")

visium <- AddMetaData(visium, combined_df)
visium$Cottrazm <- factor(visium$Cottrazm)
visium <- SetIdent(visium, value = "Cottrazm")

# subset features
features <- c("c2l_Tumor", "c2l_myCAF", "c2l_iCAF", "c2l_progenitor.iCAF")

df <- FetchData(visium, vars = features) %>%
  tibble::rownames_to_column("barcode") %>%
  left_join(
    visium@meta.data %>%
      dplyr::select(Cottrazm) %>%
      tibble::rownames_to_column("barcode"),
    by = "barcode"
  ) %>%
  dplyr::filter(grepl("^(HGD|MIC)_", Cottrazm)) %>%
  tidyr::separate(Cottrazm, into = c("Condition", "Region"), sep = "_", remove = TRUE) %>%
  tidyr::pivot_longer(all_of(features), names_to = "CellType", values_to = "score") %>%
  dplyr::mutate(
    Region    = factor(Region,    levels = c("Mal", "Bdy", "nMal")),
    CellType  = factor(CellType,  levels = features),
    Condition = factor(Condition, levels = c("HGD", "MIC"))
  )

summ <- df %>%
  group_by(Condition, Region, CellType) %>%
  summarize(
    mean_score = mean(score, na.rm = TRUE),
    se         = sd(score, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups    = "drop"
  )

# filtering myCAF, iCAF
summ_2 <- summ %>%
  dplyr::filter(CellType %in% c("c2l_myCAF", "c2l_iCAF")) %>%
  dplyr::mutate(
    CellType  = factor(CellType,
                       levels = c("c2l_myCAF", "c2l_iCAF"),
                       labels = c("myCAF", "iCAF")),
    Region    = factor(Region,
                       levels = c("Mal", "Bdy", "nMal"),
                       labels = c("intratumor", "tumor border", "peritumor")),
    Condition = factor(Condition,
                       levels = c("HGD", "MIC"),
                       labels = c("Dysplasia", "Microinvasive carcinoma"))
  )

# annotation data
annot_df <- data.frame(
  Condition   = factor(rep(c("Dysplasia", "Microinvasive carcinoma"), each = 4),
                       levels = c("Dysplasia", "Microinvasive carcinoma")),
  CellType    = rep(c("myCAF", "myCAF", "iCAF", "iCAF"), 2),
  xmin        = rep(c("intratumor", "tumor border", "intratumor", "tumor border"), 2),
  xmax        = rep(c("tumor border", "peritumor",  "tumor border", "peritumor"),  2),
  annotation  = c("***", "ns",  "***", "***",    # Dysplasia
                  "***", "***", "***", "**"),      # Microinvasive carcinoma
  y_position  = c(0.0108, 0.0112, 0.0028, 0.0030, # Dysplasia
                  0.0135, 0.0140, 0.0025, 0.0027)  # MIC     
)

# plot
p_final <- ggplot(
  summ_2,
  aes(x = Region, y = mean_score, color = CellType, group = CellType)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8, stroke = 0.2) +
  geom_errorbar(
    aes(ymin = mean_score - se, ymax = mean_score + se),
    width = 0.15, linewidth = 0.35
  ) +
  # statistics bracket
  ggsignif::geom_signif(
    data        = annot_df,
    aes(xmin        = xmin,
        xmax        = xmax,
        annotations = annotation,
        y_position  = y_position),
    manual      = TRUE,
    tip_length  = 0.01,
    color       = "black",
    size        = 0.3,
    textsize    = 2.2,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Condition, nrow = 2, scales = "free_y") +     
  scale_color_manual(
    values = c("myCAF" = col_myCAF, "iCAF" = col_iCAF)
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) + 
  labs(x = "Region", y = "Mean score") +
  guides(color = guide_legend(override.aes = list(size = 2.2, linewidth = 0.9))) +
  theme_classic(base_size = size_pt, base_family = "Arial") +
  theme(
    text             = element_text(family = "Arial"),
    line             = element_line(linewidth = 0.3),
    axis.text.x      = element_text(hjust = 0.5, size = size_pt - 2, color = "black"),
    axis.text.y      = element_text(size = size_pt - 2, color = "black"),
    axis.title.x     = element_text(size = size_pt),
    axis.title.y     = element_text(size = size_pt),
    legend.title     = element_blank(),
    legend.text      = element_text(size = size_pt - 2),
    legend.position        = "inside",
    legend.position.inside = c(0.88, 0.78), 
    legend.key.height      = unit(0.25, "cm"),
    legend.key.width       = unit(0.35, "cm"),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text       = element_text(face = "bold", size = size_pt)
  )


# Fig 4D ----
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

# Fig 4E ----

escc <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250806_escc_squamouscells_umap_tuning_sk_hj.rds")
escc_subset_Tumor = subset(escc, subset = Annotation_v2 == "Tumor")
escc_subset_Tumor = subset(escc_subset_Tumor, subset = pathology == "Dysplasia")

pw <- gmtPathways("/Users/hojin/Dropbox/db/h.all.v2025.1.Hs.symbols.gmt")

# to remove effect of ITGB1, ITGB1 is removed from each gene set.

for (pw_tmp in names(pw)) {
  gene_list <- list(pw[[pw_tmp]])
  gene_list <- list(setdiff(gene_list[[1]], "ITGB1"))
  escc_subset_Tumor <- AddModuleScore(escc_subset_Tumor, features = gene_list, name = pw_tmp)
  colnames(escc_subset_Tumor@meta.data) <- gsub(paste0(pw_tmp, "1"), pw_tmp, colnames(escc_subset_Tumor@meta.data))
}

gex <- escc_subset_Tumor@assays$RNA@data["ITGB1",]
escc_subset_Tumor$ITGB1 <- gex
metadat <- escc_subset_Tumor@meta.data

metadat <- metadat[,str_detect(colnames(metadat), pattern = "HALLMARK|ITGB1")]
gex <- metadat$ITGB1
metadat <- metadat %>% dplyr::select(-c("ITGB1"))
df <- data.frame()
for (geneset_tmp in colnames(metadat)) {
  cor_res <- cor.test(x = gex, y = metadat[[geneset_tmp]], method = "pearson")
  pval_tmp <- cor_res$p.value
  cor_tmp <- cor_res$estimate
  df_tmp <- data.frame(geneset = geneset_tmp, cor = cor_tmp, pval = pval_tmp)
  df <- rbind(df, df_tmp)
}

df$Padj <- p.adjust(df$pval, method = "BH")

hallmark_res <- df

colnames(hallmark_res) <- c("pathway", "cor", "cor_pval", "cor_Padj")

write.table(hallmark_res, "20250807_escc_tumor_ITGB1_hallmark_correlation_hj.txt", col.names = T, sep = "\t", row.names = F, quote = F)

hallmark_res_filtered <- hallmark_res %>% dplyr::filter(cor_Padj < 0.05) %>% dplyr::filter(abs(cor) > 0.1) %>% dplyr::arrange(cor)
hallmark_res_filtered$pathway <- factor(hallmark_res_filtered$pathway, levels = hallmark_res_filtered$pathway)
# plotting 

max_size_val <- max(-log10(hallmark_res_filtered$cor_Padj), na.rm = TRUE)
min_size_val <- min(-log10(hallmark_res_filtered$cor_Padj), na.rm = TRUE)

fileidentity <- "dysplasia_ITGB1"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(10), height = cm_to_inch(6), family="Arial")
hallmark_res_filtered %>% ggplot(aes(y = pathway, x = cor, size = -log10(cor_Padj))) +
  geom_point() +
  theme_bw(base_family = "Arial") +
  xlab("Correlation") +
  ylab("Pathway") +
  xlim(c(0.1, 0.25)) +
  theme(panel.border = element_rect(linewidth = 0.5), aspect.ratio = 2,
        axis.title = element_text(size = 5),
        axis.ticks = element_line(linewidth = 0.3),
        line = element_line(linewidth = 0.3),
        axis.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 3),
        legend.ticks = element_blank(),
        legend.key.height = unit(0.2, "cm"),
        legend.key.width  = unit(0.1, "cm"),
        legend.spacing.x  = unit(0.02, "cm"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -0.05, unit = "cm"))+
  guides(size = guide_legend(title = expression(-log[10] ~ "(padj)"))) + scale_size(range = c(0.5, 3), breaks = c(20, 30, 40), limits = c(min_size_val, max_size_val))
dev.off()

# Fig 4F ----

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

## 1. make viper objects -------
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

## Regulon analysis was conducted in mccleary ----

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

## 2. viper wilcoxon -------
escc_subset <- readRDS(file = "escc_tumor_viper_hj.rds")

### 2.1 normal -------
pathology_target = "Normal"
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
  y <- df_tmp %>% dplyr::filter(!pathology %in% pathology_target) %>% dplyr::pull(protein_activity)
  wilres <- wilcox.test(x = x, y = y)
  
  wilres_df_tmp <- data.frame(TF = reg_tmp,
                              mean_group1 = mean(x),
                              mean_group2 = mean(y),
                              diff = mean(x) - mean(y),
                              pvalue = wilres$p.value)
  wilres_df <- rbind(wilres_df, wilres_df_tmp)
}

wilres_df$adj_pvalue <- p.adjust(wilres_df$pvalue, method = "bonferroni")

write.table(wilres_df, "20250807_escc_viper_normal_others_wilcoxon_hj.txt", col.names = T, sep = "\t", row.names = F, quote = F)

### 2.2 dysplasia ------
pathology_target = "Dysplasia"
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
  y <- df_tmp %>% dplyr::filter(!pathology %in% pathology_target) %>% dplyr::pull(protein_activity)
  wilres <- wilcox.test(x = x, y = y)
  
  wilres_df_tmp <- data.frame(TF = reg_tmp,
                              mean_group1 = mean(x),
                              mean_group2 = mean(y),
                              diff = mean(x) - mean(y),
                              pvalue = wilres$p.value)
  wilres_df <- rbind(wilres_df, wilres_df_tmp)
}

wilres_df$adj_pvalue <- p.adjust(wilres_df$pvalue, method = "bonferroni")

write.table(wilres_df, "20250807_escc_viper_dysplasia_others_wilcoxon_hj.txt", col.names = T, sep = "\t", row.names = F, quote = F)

## micro ------
pathology_target = c("Microinvasive carcinoma")
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
  y <- df_tmp %>% dplyr::filter(!pathology %in% pathology_target) %>% dplyr::pull(protein_activity)
  wilres <- wilcox.test(x = x, y = y)
  
  wilres_df_tmp <- data.frame(TF = reg_tmp,
                              mean_group1 = mean(x),
                              mean_group2 = mean(y),
                              diff = mean(x) - mean(y),
                              pvalue = wilres$p.value)
  wilres_df <- rbind(wilres_df, wilres_df_tmp)
}

wilres_df$adj_pvalue <- p.adjust(wilres_df$pvalue, method = "bonferroni")

write.table(wilres_df, "20250807_escc_viper_micro_others_wilcoxon_hj.txt", col.names = T, sep = "\t", row.names = F, quote = F)

## macro ------
pathology_target = c("Macroinvasive carcinoma")
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
  y <- df_tmp %>% dplyr::filter(!pathology %in% pathology_target) %>% dplyr::pull(protein_activity)
  wilres <- wilcox.test(x = x, y = y)
  
  wilres_df_tmp <- data.frame(TF = reg_tmp,
                              mean_group1 = mean(x),
                              mean_group2 = mean(y),
                              diff = mean(x) - mean(y),
                              pvalue = wilres$p.value)
  wilres_df <- rbind(wilres_df, wilres_df_tmp)
}

wilres_df$adj_pvalue <- p.adjust(wilres_df$pvalue, method = "bonferroni")

write.table(wilres_df, "20250807_escc_viper_macro_others_wilcoxon_hj.txt", col.names = T, sep = "\t", row.names = F, quote = F)

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

## 3. dot plots ------
tb1 = read.table("20250807_escc_viper_normal_others_wilcoxon_hj.txt", header = T, sep = "\t")
tb2 = read.table("20250807_escc_viper_dysplasia_others_wilcoxon_hj.txt", header = T, sep = "\t")
tb3 = read.table("20250807_escc_viper_micro_others_wilcoxon_hj.txt", header = T, sep = "\t")
tb4 = read.table("20250807_escc_viper_macro_others_wilcoxon_hj.txt", header = T, sep = "\t")
tb5 = read.table("20250807_escc_viper_dys_normal_wilcoxon_hj.txt", header = T, sep = "\t")


tb1$stage <- "normal"
tb2$stage <- "dysplasia"
tb3$stage <- "micro"
tb4$stage <- "macro"


total_pos_tb = tb5 %>% dplyr::filter(diff > 0) %>% dplyr::filter(adj_pvalue < 0.05) %>% dplyr::arrange(adj_pvalue)
total_pos_tb %>% dplyr::filter(total_pos_tb$TF == "HIF1A")

# ORA

library(org.Hs.eg.db)
library(msigdbr)
library(clusterProfiler)
library(stringr)

##enricher
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
hs <- org.Hs.eg.db
ENTREZID <- AnnotationDbi::select(hs, keys = total_pos_tb$TF, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>% dplyr::select(ENTREZID) %>% na.omit() %>% pull()

res <- clusterProfiler::enricher(ENTREZID,
                                 TERM2GENE = m_t2g,
                                 pvalueCutoff = 1,
                                 qvalueCutoff = 1,
                                 minGSSize = 1)
res <- res@result
res <- res %>% dplyr::filter(p.adjust < 0.05)

saveRDS(res, "20250810_tumor_viper_normal_dys_hallmark_hj.rds")

res <- readRDS("20250810_tumor_viper_normal_dys_hallmark_hj.rds")

df <- data.frame()
for (i in 1:3) {
  mygenes <- str_split(res$geneID[i], pattern = "/")[[1]]
  SYMBOL <- AnnotationDbi::select(hs, keys = mygenes, columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")
  df_tmp <- data.frame(gene = SYMBOL$SYMBOL, hallmark = res$Description[i])
  df <- rbind(df, df_tmp)
}

df$hallmark <- gsub(df$hallmark, pattern = "HALLMARK_", replacement = "") %>% gsub(df$hallmark, pattern = "_", replacement = " ")

saveRDS(df, "20250810_tumor_viper_hallmark_genes_hj.rds")

df <- readRDS("20250810_tumor_viper_hallmark_genes_hj.rds")
mymat <- matrix(rep(NA, length(unique(df$gene))*3), nrow = length(unique(df$gene)), ncol = 3)
rownames(mymat) <- unique(df$gene)
colnames(mymat) <- gsub(res$ID, pattern = "HALLMARK_", replacement = "") %>% gsub(res$ID, pattern = "_", replacement = " ")

mytb <- tb5 %>% dplyr::filter(TF %in% df$gene)

for (i in 1:nrow(df)) {
  
  gene_tmp <- df[i,]$gene
  hallmark_tmp <- df[i,]$hallmark
  
  mymat[gene_tmp, hallmark_tmp] <- mytb %>% dplyr::filter(TF == gene_tmp) %>% dplyr::pull(diff)
  
}

df_long <- as.data.frame(mymat) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "hallmark", values_to = "score")
df_long$gene <- factor(df_long$gene, levels = rownames(mymat))


fileidentity <- "viper_tumor_normal_dys_dotplot"
cairo_pdf(paste0("20250818", "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(11), height = cm_to_inch(7), family="Arial")
ggplot(df_long, aes(x = gene, y = hallmark)) +
  geom_point(aes(size = abs(score), color = score)) +
  scale_size_continuous(range = c(0.5, 2)) +
  scale_color_gradient2(low = "ivory", high = "red") +
  theme_classic(base_family = "Arial") +
  theme(line = element_line(linewidth = 0.3), aspect.ratio = 0.25,
        axis.title = element_text(size = 5),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 4.5, colour = "black"),
        axis.text.y = element_text(size = 4.5, colour = "black"),
        legend.key.size = unit(0.15, "cm"),  
        legend.position = "bottom",
        legend.title = element_text(size = 4),
        legend.text  = element_text(size = 3),
  ) +
  labs(size = "Strength", color = "Difference of protein activity")
dev.off()

# Fig 4G ----

library(Seurat)
library(ggplot2)

fibroblast_final_sk <- readRDS("~/Dropbox/20250724_ESCC_scRNAseq_annotation_T_cells/RData/20250812_fibroblast_final_sk.rds")
colors <- CellChat::scPalette(30)
DimPlot(fibroblast_final_sk, cols = colors, label = TRUE) + theme(aspect.ratio = 1) + NoLegend()

# Fig 4H ----

library(monocle3)
library(SeuratWrappers)

fibroblast_final_sk <- readRDS("~/Dropbox/20250724_ESCC_scRNAseq_annotation_T_cells/RData/20250812_fibroblast_final_sk.rds")
DimPlot(fibroblast_final_sk, cells.highlight = 'EGS-22-N-SC-5p_CTGAAACAGGACGAAA-1') + theme(aspect.ratio = 1)
fibroblast_final_sk_cds <- as.cell_data_set(fibroblast_final_sk)
fibroblast_final_sk_cds <- cluster_cells(fibroblast_final_sk_cds,reduction_method = "UMAP")
fibroblast_final_sk_cds <- learn_graph(fibroblast_final_sk_cds)
fibroblast_final_sk_cds <- order_cells(fibroblast_final_sk_cds, root_cells = 'EGS-22-N-SC-5p_CTGAAACAGGACGAAA-1')
plot_cells(fibroblast_final_sk_cds, color_cells_by = "pseudotime", label_cell_groups = F,label_groups_by_cluster=F,show_trajectory_graph = TRUE, group_label_size=0,graph_label_size = 0,labels_per_group=0) + theme(aspect.ratio = 1)

# Fig 4I ----
library(fgsea)
library(ggplot2)
library(dplyr)
library(Seurat)
library(dplyr)
library(AUCell)
library(SCENIC)
library(ComplexHeatmap)
library(stringr)
library(showtext)
library(org.Hs.eg.db)
library(msigdbr)
library(clusterProfiler)
library(openxlsx)
library(circlize)

setwd("/Users/jihyunkim/Library/CloudStorage/Dropbox/escc_figure_share/Fig4")

Stromal <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Latest_Figures/Figures_20250806/object/20250728_escc_Fib_fDC_FRC_hj.rds")
Stromal <- subset(Stromal, Annotation_v2 %in% c("NMF","progenitor iCAF", "iCAF", "myCAF"))
log2fc_df <- FoldChange(Stromal, ident.1 = "myCAF", group.by = "Annotation_v2")
res <- readRDS("/Users/jihyunkim/Library/CloudStorage/Dropbox/ESCC/Latest_Figures/Figures_20250806/object/20250809_myCAF_DEG_mast_re.rds")
res$gene <- rownames(res); log2fc_df$gene <- rownames(log2fc_df)
res <- merge(res, log2fc_df, by = "gene", all.x = T)
colnames(res) <- c("gene", "p_val", "p_val_adj", "avg_log2FC", "pct.1", "pct.2") 


source("/Users/jihyunkim/Library/CloudStorage/Dropbox/treg_2025/20230816_fgsea_function.R")
myCAF.hallmark <- fgsea_hallmark_function("/Users/jihyunkim/Library/CloudStorage/Dropbox/escc_figure_share/Fig4/20250809_myCAF_MAST_RE.csv") %>% filter(padj < 0.05)
myCAF.kegg <- fgsea_kegg_function("/Users/jihyunkim/Library/CloudStorage/Dropbox/escc_figure_share/Fig4/20250809_myCAF_MAST_RE.csv") %>% filter(padj < 0.05)
# ====== I/O & Parameters ======
gmt_file   <- "/Users/jihyunkim/gsea/c5.go.bp.v2025.1.Hs.symbols.gmt"
topN_each  <- 25 
y_breaks   <- seq(-2, 2, by = 1) 
y_limits   <- NULL  

# ====== Pathways (GMT) ======
pathways <- fgsea::gmtPathways(gmt_file)   

# ====== Rank vector (from 'res') ======
# res: data.frame with columns gene, p_val, avg_log2FC
myCAF_MAST_RE <- res %>%
  mutate(
    statistics = qnorm(p_val/2, lower.tail = FALSE) * sign(avg_log2FC),
    max2 = max(statistics[statistics !=  Inf]),
    min2 = min(statistics[statistics != -Inf]),
    stat_corr = case_when(
      statistics ==  Inf ~ max2 + 1,
      statistics == -Inf ~ min2 - 1,
      TRUE               ~ statistics
    )
  )

ranks_vec <- myCAF_MAST_RE$stat_corr
names(ranks_vec) <- res$gene
ranks_vec <- ranks_vec[is.finite(ranks_vec)]

# ====== FGSEA ======
mycaf_gobp <- fgsea(pathways = pathways, stats = ranks_vec)

# ====== Group mapping (no vessel) ======
tgfbeta <- c(
  'GOBP_NEGATIVE_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY',
  'GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA_STIMULUS',
  'GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY',
  'GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA'
)
ecm <- c('GOBP_EXTRACELLULAR_MATRIX_DISASSEMBLY')
adhesion <- c(
  'GOBP_FOCAL_ADHESION_ASSEMBLY',
  'GOBP_CELL_MATRIX_ADHESION',
  'GOBP_REGULATION_OF_CELL_ADHESION'
)
collagen <- c(
  'GOBP_REGULATION_OF_COLLAGEN_FIBRIL_ORGANIZATION',
  'GOBP_POSITIVE_REGULATION_OF_COLLAGEN_METABOLIC_PROCESS',
  'GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS',
  'GOBP_COLLAGEN_CATABOLIC_PROCESS',
  'GOBP_COLLAGEN_FIBRIL_ORGANIZATION',
  'GOBP_COLLAGEN_METABOLIC_PROCESS'
)

mycaf_gobp <- mycaf_gobp %>%
  mutate(
    Groups = case_when(
      pathway %in% ecm      ~ "ECM remodelling",
      pathway %in% tgfbeta  ~ "TGF beta",
      pathway %in% adhesion ~ "Focal adhesion",
      pathway %in% collagen ~ "Collagen formation",
      TRUE                  ~ "Others"
    ),
    Groups = factor(Groups,
                    levels = c("ECM remodelling","TGF beta","Focal adhesion","Collagen formation","Others")
    )
  )

# ====== Plot prep ======
# fgsea 결과 컬럼: pathway, pval, padj, NES, size, etc.
plot_df <- mycaf_gobp %>%
  transmute(
    pathway, NES, padj, Groups,
    ranks = rank(-NES, ties.method = "first"),
    sig   = ifelse(padj < 0.05, "< 0.05", ">= 0.05")
  )

bg_curve <- arrange(plot_df, ranks)

grp_cols <- c(
  "Collagen formation" = "#F39C12",
  "TGF beta"           = "#E74C3C",
  "Focal adhesion"     = "#7D3C98",
  "ECM remodelling"    = "#6E2C00",
  "Others"             = "grey60"
)


fg <- plot_df %>%
  filter(Groups != "Others") %>%
  group_by(Groups) %>%
  slice_min(ranks, n = topN_each, with_ties = FALSE) %>%
  ungroup()

# ====== Plot ======
p <- ggplot() +

  geom_point(data = bg_curve, aes(ranks, NES),
             color = "grey82", size = 0.6, alpha = 0.35) +
  geom_path (data = bg_curve, aes(ranks, NES),
             color = "grey85", linewidth = 0.5, alpha = 0.25) +

  geom_point(data = fg, aes(ranks, NES, shape = sig),
             size = 4.3, color = "white", stroke = 0.8, alpha = 0.95) +

  geom_point(data = fg,
             aes(ranks, NES, color = Groups, shape = sig),
             size = 3.1, stroke = 0.5) +
  scale_shape_manual(values = c("< 0.05" = 16, ">= 0.05" = 8), name = "p.adj") +
  scale_color_manual(values = grp_cols, drop = FALSE, name = NULL) +
  labs(x = "ranks", y = "Normalized Enrichment Score") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = c(0.83, 0.75),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  ) + scale_y_continuous(breaks=seq(-2, 3, 1),trans = scales::exp_trans())

setwd("/Users/jihyunkim/Library/CloudStorage/Dropbox/escc_figure_share/Fig4/")
ggsave("20250810_mycaf_enriched_gobp_jk.svg", p, width = 5, height = 4)


