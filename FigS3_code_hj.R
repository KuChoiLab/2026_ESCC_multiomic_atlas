# Fig S3A ----

library(miloR)
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(showtext)
library(ggbeeswarm)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

date <- "20250807"
project_name <- "escc"
fileidentity <- "miloR_plotDAbeeswarm"
myname <- "jk"
size_pt <- 7

output_dir <- "/Users/jihyunkim/Library/CloudStorage/Dropbox/escc_figure_share/Fig3/miloR_abundance_plot"
setwd(output_dir)

ESCC_results <- readRDS(file.path(output_dir, "20250807_miloR_obj.rds"))
miloR_plot <- readRDS(file.path(output_dir, "20250807_miloR_plot.rds"))

miloR_plot_final <- miloR_plot +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(colour = "black", size = size_pt, angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = size_pt),
    axis.title.x = element_text(colour = "black", size = size_pt),
    axis.title.y = element_blank(),
    legend.text = element_text(size = size_pt - 1),
    legend.title = element_blank(),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.1, "cm"),
    legend.spacing.x = unit(0.02, "cm"),
    legend.position = "right",
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    plot.title = element_blank()
  )

pdf_file <- paste0(date, "_", project_name, "_", fileidentity, "_", myname, ".pdf")
svg_file <- paste0(date, "_", project_name, "_", fileidentity, "_", myname, ".svg")
rds_file <- paste0(date, "_", project_name, "_", fileidentity, "_", myname, ".rds")

width_cm <- 6 * 1.5 
height_cm <- 6 * 1.2

cairo_pdf(pdf_file, width = cm_to_inch(width_cm), height = cm_to_inch(height_cm), family = "Arial")
print(miloR_plot_final)
dev.off()

# Fig S3B ----
library(ggplot2)
library(Seurat)

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250806_escc_squamouscells_umap_tuning_sk_hj.rds")

old <- GeomViolin$default_aes$linewidth
update_geom_defaults("violin", list(linewidth = 0.3))

VlnPlot(seurat_obj, features = c("COL17A1", "KRT13", "GNGT1", "POLQ", "LGR5"), group.by = "Annotation_v3", cols = c("#F8766D", "#0CB702", "blue4", "magenta2", "gold2", "#619CFF") , pt.size = 0, ncol = 5) & theme_classic(base_family = "Arial") & theme(aspect.ratio = 0.7, axis.text = element_text(size = 5), axis.title = element_text(size = 5), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, size = 6), line = element_line(linewidth = 0.3)) & NoLegend()


# Fig S3C ----

library(CytoTRACE2)

cytotrace2_squamous_final <- cytotrace2(squamous_final, species = 'human', is_seurat = TRUE, slot_type = 'counts', batch_size = 10000, parallelize_models = TRUE, ncores = 4, seed = 1234)
cytotrace2_squamous_final_plots <- plotData(cytotrace2_squamous_final, is_seurat = TRUE, seed = 1234, pc_dims = 30)
cytotrace2_squamous_final_plots$CytoTRACE2_UMAP


# Fig S3D ----

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(CellChat)
library(showtext)
library(MAST)
library(SCENIC)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

date = "20250804"
project_name = "escc"

size_pt <- 7

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/scenic/TIC/")

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250806_escc_squamouscells_umap_tuning_sk_hj.rds")

seurat_obj@meta.data <- seurat_obj@meta.data %>%
  dplyr::mutate(pathology = case_when(pathology == "Microinvasive" ~ "Microinvasive carcinoma",
                                      pathology == "Macroinvasive" ~ "Macroinvasive carcinoma",
                                      T ~ pathology
  ))

seurat_obj$pathology <- factor(seurat_obj$pathology, levels = c("Normal", "Dysplasia", "Microinvasive carcinoma", "Macroinvasive carcinoma"))

set.seed(1234); seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)

# SCENIC 
library(Seurat)
library(dplyr)
library(AUCell)
library(SCENIC)
library(ComplexHeatmap)
library(stringr)
library(showtext)
library(ComplexHeatmap)
library(colorRamp2)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

## 1. wilcoxon ----

wilcoxon_res <- list()
regulonAUC <- readRDS('/Users/hojin/Dropbox/project/ESCC/submit/analysis/scenic/sq/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")
metadat <- seurat_obj@meta.data
regulonmat <- regulonAUC@assays@data$AUC

res_df <- data.frame()
TFs <- rownames(regulonmat) %>% unique()
#TFs <- TFs[!str_detect(TFs, "_extended")]
metadat <- seurat_obj@meta.data
metadat$cellname <- rownames(metadat)

res_df <- data.frame()
for (i in 1:length(TFs)) {
  tf_tmp <- TFs[i]
  react <- regulonmat[tf_tmp, ]
  df_tmp <- data.frame(cellname = names(react), regulon_activity = react)
  df_tmp <- merge(df_tmp, metadat, by.x = "cellname", by.y = "cellname", all.x = T)
  
  for (celltype_tmp in c("Tumor initiating cell", "Transitional suprabasal", "LGR5+ cancer stem cell", "Tumor")) {
    df1 <- df_tmp %>% dplyr::filter(Annotation_v3 == celltype_tmp)
    print(celltype_tmp)
    #background <- setdiff(c("Transitional suprabasal", "Tumor initiating cell", "Tumor", "LGR5+ cancer stem cell"), celltype_tmp)
    #print(background)
    df2 <- df_tmp %>% dplyr::filter(Annotation_v3 != celltype_tmp)
    res <- wilcox.test(df1$regulon_activity, df2$regulon_activity)
    pval_tmp <- res$p.value
    res_tmp <- data.frame(TF = tf_tmp,
                          group1 = celltype_tmp,
                          group1_mean = mean(df1$regulon_activity),
                          group2_mean = mean(df2$regulon_activity),
                          pval = pval_tmp)
    res_df <- rbind(res_df, res_tmp)
  }
}

res_df$padj <- p.adjust(res_df$pval, method = "bonferroni")

saveRDS(res_df, "tic_test_extended_hj.rds")

res_df <- readRDS("tic_test_extended_hj.rds")

res_df1 <- res_df %>% dplyr::filter(group1 == "Tumor initiating cell")
res_df2 <- res_df %>% dplyr::filter(group1 == "Transitional suprabasal")
res_df3 <- res_df %>% dplyr::filter(group1 == "LGR5+ cancer stem cell")
#res_df4 <- res_df %>% dplyr::filter(group1 == "Tumor")

res_df1$padj <- p.adjust(res_df1$pval, method = "bonferroni")
res_df2$padj <- p.adjust(res_df2$pval, method = "bonferroni")
res_df3$padj <- p.adjust(res_df3$pval, method = "bonferroni")
#res_df4$padj <- p.adjust(res_df4$pval, method = "bonferroni")

res_df <- rbind(res_df1, res_df2)
#res_df <- rbind(res_df1)
res_df <- res_df %>% dplyr::filter(padj < 0.05)

quant <- res_df$padj %>% quantile()

res_df_filtered <- res_df %>% dplyr::mutate(log2FC = log2(group1_mean/group2_mean)) %>% dplyr::filter(group2_mean < group1_mean)

tmp <- res_df_filtered %>% dplyr::group_by(TF) %>% dplyr::summarise(n = n()) %>% dplyr::filter(n == 1) %>% dplyr::pull(TF)


tf1 <- res_df_filtered %>% dplyr::filter(group1 == "Tumor initiating cell") %>% dplyr::pull(TF)
tf2 <- res_df_filtered %>% dplyr::filter(group1 == "Transitional suprabasal") %>% dplyr::pull(TF)

mytfs <- setdiff(tf1, tf2)

res_df_filtered <- res_df_filtered %>% dplyr::filter((TF %in% mytfs) & (log2FC > 0.5))
res_df_filtered <- res_df_filtered %>% dplyr::filter(padj < 0.00001)
mytfs <- res_df_filtered$TF

regulonAUC_scaled <- getAUC(regulonAUC)[mytfs, colnames(seurat_obj)]
regulonAUC_scaled <- t(scale(t(regulonAUC_scaled), center = T, scale=T))

# plot

celltypes <- seurat_obj@meta.data %>% select(Annotation_v3)
celltypes$TIC <- factor(celltypes$Annotation_v3, levels = c("Basal", "Suprabasal", "Transitional suprabasal", "Tumor initiating cell", "LGR5+ cancer stem cell", "Tumor"))
selectedResolution <- "TIC"
cellsPerTypes <- split(rownames(celltypes), celltypes[,selectedResolution]) 
#cellorder <- unlist(cellsPerTypes)

set.seed(1234);cellorder <- unlist(lapply(names(cellsPerTypes),
                                          function(ct) {
                                            unlist(lapply(split(cellsPerTypes[[ct]], seurat_obj@meta.data[cellsPerTypes[[ct]], "platform"]),sample))}))

cell_type <- data.frame(
  "cell type" = unlist(mapply(rep, names(cellsPerTypes), lapply(cellsPerTypes, length))))
regulonAUC_scaled_order <- regulonAUC_scaled[mytfs, cellorder]
row.names(cell_type) <- colnames(regulonAUC_scaled_order)
palette_length = 100

date <- "20250810"
project_name <- "escc"
cm_to_inch <- function(cm) cm / 2.54
fileidentity <- "TIC_scenic_v2"

library(CellChat)
mycol <- c("Basal" = "#F8766D", "Suprabasal" = "#0CB702", "Transitional suprabasal" = "blue4", "Tumor initiating cell" = "magenta2", "LGR5+ cancer stem cell" = "gold2", "Tumor" = "#619CFF")

platform_col <- c("scRNA" = "grey", "snRNA" = "orange3")
platform_vec <- seurat_obj@meta.data[colnames(regulonAUC_scaled_order), "platform"] |> as.character()
platform_vec <- factor(platform_vec, levels = names(platform_col))

library(colorRamp2)
ha = HeatmapAnnotation(Platform = platform_vec,
                       Celltype = unlist(mapply(rep, names(cellsPerTypes), lapply(cellsPerTypes, length))),
                       col = list(
                         Platform = platform_col,
                         Celltype = mycol
                       ),
                       simple_anno_size = unit(0.3, "cm"),
                       show_legend = FALSE,
                       #annotation_name_side = NULL,
                       annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param = list(
                         title = NULL,
                         title_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                         labels_gp = gpar(fontsize = 5, fontfamily = "Arial")
                       ))


col_fun = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "red"))

png(filename = paste0(date, "_", project_name, "_", fileidentity, "_hj.png"),
    width = 5, height = 4, units = "in", res = 600, family = "Arial")
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
                             row_names_gp = gpar(fontsize = 7, fontfamily = "Arial"),
                             column_names_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                             heatmap_legend_param = list(
                               title = "Score",
                               title_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                               labels_gp = gpar(fontsize = 5, fontfamily = "Arial")
                             ),show_heatmap_legend = FALSE,
                             top_annotation = ha)
draw(ht, use_raster = TRUE,  show_annotation_legend = FALSE, annotation_legend_list = list(
  gpar(fontfamily = "Arial", fontsize = 5)
))
dev.off()

## violin plot ----
GeomViolin$default_aes <- modifyList(
  GeomViolin$default_aes,
  list(linewidth = 0.2)
)

df <- readRDS("20250820_TIC_regact_wilcoxon_tumor_ts_hj.rds")
df <- df %>% dplyr::filter(test_tumor_padj < 0.05) %>% dplyr::filter(test_ts_padj < 0.05)

df %>% dplyr::filter(test_tumor_padj < 0.0000000001)
df %>% dplyr::filter(test_ts_padj < 1e-50)

seurat_obj_subset_regulon_log2FC <- subset(seurat_obj_subset, subset = Annotation_v3 != "LGR5+ cancer stem cell")
seurat_obj_subset_regulon_log2FC$Annotation_v3 %>% table()

mat <- GetAssayData(seurat_obj_subset_regulon_log2FC, assay = "scenic", slot = "data")
groups <- seurat_obj_subset_regulon_log2FC$Annotation_v3
group1 <- rowMeans(mat[, groups == "Tumor initiating cell", drop = FALSE])
group2 <- rowMeans(mat[, groups != "Tumor initiating cell", drop = FALSE])
log2FC <- log2((group1 + 1e-6) / (group2 + 1e-6))

mylog2FC_sort_df <- as.data.frame(log2FC[c(df$TF)] %>% sort(decreasing = T))
colnames(mylog2FC_sort_df) <- "log2FC"
write.table(mylog2FC_sort_df, "20250820_escc_TIC_log2FC_table_hj_v2.txt", col.names = T, sep = "\t", row.names = T, quote = F)

mylog2FC_sort <- log2FC[c(df$TF)] %>% sort(decreasing = T) %>% names()
library(scales)
mycol <- c("Basal" = "#F8766D", "Suprabasal" = "#0CB702", "Transitional suprabasal" = "blue4", "Tumor initiating cell" = "magenta2", "LGR5+ cancer stem cell" = "gold2", "Tumor" = "#619CFF")

VlnPlot(seurat_obj_subset, features = c(mylog2FC_sort), pt.size = 0, group.by = "Annotation_v3", ncol = 4, cols = mycol) & ylab("Regulon activity") & theme_classic(base_family = "Arial") & theme(aspect.ratio = 0.5, line = element_line(linewidth = 0.2), axis.text = element_text(size = 5), axis.title.x = element_blank(), axis.title = element_text(size = 5), plot.title = element_text(size = 5, hjust = 0.5), axis.text.x = element_blank(),axis.title.y = element_text(size = 4), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 3)) & NoLegend() & scale_y_continuous(labels = number_format(accuracy = 0.01))

# Fig S3E ----
library(Seurat)
library(ggplot2)
library(cowplot)

# module gene lists
mod0 <- c(
  "MKI67","DIAPH3","NEURL1B","DHFR","FSCN1","KANK2","AURKB","ABCC5","CENPE","TPX2",
  "UBE2C","TOP2A","BIRC5","HMGB3","NUSAP1","ATAD2","LPAR2","CHAF1A","CENPF","PLAU",
  "ASPM","LRP8","GPI","POLQ","SAMD1","KLF7","TONSL","TPPP","MYBL1","BRIP1",
  "KIF18B","ADGRE5","SPSB1","PAG1","IFI27","SLC2A1","CEL","TAGLN2","ST6GALNAC4","PREX1",
  "ARFGAP1","RAB40B","RTEL1","ZC3H3","MFHAS1","NAV1","MOB3A","SULF2","MYH9","PPIF",
  "PSMB3","DOC2B","AKR1C2","RUNX1","EIF4A3","ID2","EPS8","ZNF593","INPP5D","EPAS1",
  "RASSF2","DAPK1","CCND1","CYBA","ARHGAP23","SMAD3","IKZF3","PTTG1","RAN","RAB40C",
  "CTSD","GSTM3","ALDH1A1","DOT1L","FAAP100","WDR34","MRPL21","TMEM44","RPL22L1","TUBA1B",
  "NEDD4L","CUX1","RDH10","GPX2","MED1","ITGA4","SOX13","APOL2","FRMD8","IGHMBP2",
  "CTTN","PDLIM7","HSP90AA1","GLUL","FAM83E","AGRN","HCLS1","LAPTM5","STARD3","PHLDB2"
)

mod1 <- c(
  "KRT14","SERPINE2","DLL1","LAMB4","PCDH7","COL8A1","LGR6","SLC1A3",
  "NOVA1","IL1R2","DCN","GEM","PREX2","VIT","SLC7A5","COL14A1","BOC",
  "MAATS1","PAPPA","LINC02232","ZBTB16","SLC6A6","AKR1C1","PTPRD","CPA6",
  "CXCL14","SOX5","IL11RA","AC055758.2","FGFR1","CLIP2","AC124947.1",
  "COL7A1","CDH13","FHOD3","KRT15","TENM1","SNCA","GRIK4","PTHLH",
  "ACSL4","CACNB4","TNS3","LINC01127","PCDH11X","HIP1","ADAMTS9","BNC2",
  "TMTC1","WDR49","FLRT2","CDH22","AHR","NRG1","TRO","CYYR1","SFRP1",
  "TNNI2","SLC25A48","ANOS1","BACH2","PDE1A","NKD1","ST6GAL2","CCBE1",
  "TNC","ARMC9","CCDC178","KREMEN2","NAV3","SH3RF3","KCNMA1","ABCC3",
  "FGF14","SLC35F3","DSEL","SYNPO2","STXBP6","VWC2","NLGN4Y","GUCY1A2",
  "KLHL14","C1R","WIF1","AJAP1","MAGI2","CCSER1","CDON","MAN1C1","BICC1",
  "NWD1","DOCK8","POSTN","SUCLG2-AS1","APCDD1","SOBP","DISC1"
)

mod2 <- c(
  "AC024230.1","AC106799.2","TEX41","MME","AC005355.1","SLC52A1","CASC9",
  "AL078590.2","LINC01980","IGFBP3","GNGT1","CECR2","IQCM","VCAM1",
  "BANK1","IGHG1","RGS6","AC093895.1","AL139327.2","AC233296.1","ANO1",
  "SLC44A5","AC009264.1","IGF2BP2","MS4A1","AC011287.1","DYNC1I1",
  "MIR2052HG","WNK2","IDO1","AC116049.2","EXOC3L4","AL583808.1","FIRRE",
  "AC022031.2","ODC1","PTPRQ","LINC02428","AL117329.1","EPPK1",
  "AC011632.1","GALNT13","FDCSP","IGSF1","POU6F2","CACNA1B","RUNX3",
  "AC005993.1","AC109454.3","PIF1","KIF14","CHST1","HRK","SHOX2","KYNU",
  "AC009262.1","POTEF","LTF","LSAMP","PIK3R3","PAX5","AC019197.1",
  "RANBP17","ABCA13","PPP1R9A","COL19A1","SCN3A","SCN2A","CHRM3",
  "ADAM28","DTNA","SPIB","PPM1H","ATP13A5","AL033530.1","MIR9-3HG",
  "DNAH5","VNN1","FAM20C","TICRR","BBOX1-AS1","IGHG4","AC120193.1","CIT",
  "AC018697.1","BIRC3","NKAIN2","LINC01524","EPHB2","IGHG2","IGHG3",
  "UHRF1","AC010343.3","AC011997.1","SLC43A2","SYT1","CCL20","MUC16",
  "CENPW","RIMS2","ADAMDEC1","FGF12"
)

mod3 <- c(
  "KRT1","GBP6","BARX2","GSTA1","MAPT","LY6D","TSPAN5","LINC01605","ADH7","FABP5",
  "DNASE1L3","NMU","FA2H","FGD3","AKR1C3","CES1","ESPN","SERPINB4","SERPINB3","SAPCD2",
  "CTTNBP2","FRMD6-AS2","EPB41L1","AQP3","RHOV","RHCG","ITPRIP","CPNE4","BICDL1","KRT6C",
  "FOSB","RGMA","ARL4D","ALDH3B2","MGLL","AKR1B10","HS6ST2","ADCY7","OCA2","BTBD11",
  "SLC4A4","UNC5B","NOTCH3","ANKRD33B","KRT13","FYB1","DSP","NEBL","ZMIZ1","SAMD5",
  "PLD1","AL033504.1","PTGR1","AC019349.1","NECTIN1","WNK4","KLF5","HS3ST3A1","KRT6B","AHNAK2",
  "SERPINB5","KRT6A","NTN1","EGFL6","BAIAP2","PLA2G4E","RNF152","DUSP2","CAMK1D","SYT7",
  "SDC4","MTSS1","AC138305.1","CRABP2","C3orf67","XYLT1","AHNAK","GRAP2","SEMA4B","UGDH",
  "DSC2","MX2","HES1","SEMA3F","PERP","CDKN2A","PHLDA2","RECQL4","PHACTR1","TXN",
  "EEPD1","SEMA5A","ADAM19","DDIT4","PLA2G4E-AS1","EDAR","SRCIN1","ARHGAP23","HSPB1","OSBPL10"
)

mod4 <- c(
  "TTTY14","CDA","CYP4F22","LINC01170","TGM1","GJB2","SERPINB11","IVL","SLC5A1","CIDEA",
  "RRAD","GJB6","SERPINB2","CLCA4","LINC00536","SCD","CLDN4","MBOAT2","TNFAIP2","LINC00278",
  "CLVS1","KRT4","SULT2B1","AC078923.1","PCSK6","NDRG1","VNN3","CAPN14","PGBD5","HS3ST1",
  "KLK7","IFFO2","KLK13","DDIT4-AS1","SERPINB1","CSTA","ID1","LYPD3","CD24","GRHL1",
  "PDZRN3","IL1RN","TMPRSS11F","GPR78","ADM","GRHL3","TMEM184A","SPINK5","MN1","GAS7",
  "MAB21L3","PAX1","FMO2","KRT6B","SLC24A3","ATP10B","B4GALNT3","IGFL2-AS1","MALL","IL36G",
  "MCF2L2","UNC5B-AS1","KRT16","PLEKHN1","PRIMA1","HK2","ATP13A4","PLEKHA7","HBEGF","AL160408.1",
  "DUOX2","IGFL1","SPRR3","PLA2G4E","FGFBP1","SBSN","DMKN","PLA2G4E-AS1","APOBEC3A","CRNN",
  "B3GNT5","PHLDA2","FUT6","RAB11FIP1","SHANK2","CDH19","CD55","A2ML1","KRT6C","TMPRSS11A",
  "MUC21","DSC2","IL34","KRT6A","CRYBG2","EPHA4","EVPL","NBEAL2","TMEM45B","WNK4"
)

mod5 <- c(
  "CEACAM7","FETUB","ATP12A","CNFN","AC093515.1","GNG4","KLK6","LINC02487","LINC01705","HOPX",
  "KRT24","FAM3D","CLDN10-AS1","PLA2G4D","EREG","CEACAM5","CLDN10","CRISP3","GPRC5A","NRXN1",
  "SLC8A1-AS1","PRSS3","SHROOM3","SLC6A14","MUC4","PLEKHA6","PRSS27","SLCO4A1","F3","CDH26",
  "IL12A-AS1","SERPINB9","MSLN","TMPRSS11B","ARNTL2","ITGAX","CXCL8","SLC15A1","AC024597.1","LIPH",
  "EMP1","ABLIM3","AIF1L","NUAK2","PLAC8","MAL","ALOX12","AC068631.1","FAM83A","CD177",
  "ECM1","PRDM1","PNLIPRP3","EPS8L1","TGM3","NCCRP1","CALB2","KIFC3","LYPD2","IGF2BP3",
  "SULT1B1","PAX1","TMPRSS11D","FAM25A","MROH6","TRPV4","TMPRSS2","GPR78","KRT23","LRMP",
  "HCG22","USP6NL","ADIRF","KLHL2","AP006621.4","MYEOV","CAPN14","CRCT1","BCAT1","ADGRL2",
  "C15orf48","AQP5","CEACAM6","TMPRSS11E","BCAS1","HEPHL1","ST8SIA6","GSTT2B","CRNN","DOCK3",
  "S100P","AP005233.2","ATP6V1C2","ARHGAP27","KRT78","KIF21A","UPP1","WFDC12","AREG","CEACAM1",
  "AMTN","TPRG1","AL157938.2","EGLN3","MXD1","PADI1","MUCL3","RANBP9","DIO2-AS1","HILPDA"
)

mod6 <- c(
  "LINC02303","KRT79","IL1A","MT1G","DHRS9","FLG","MT1H","SFTA2","IL36A","SPRR2D",
  "MUC22","LMO7DN","ATG9B","LCN2","EDN2","LCE3D","SPNS2","PADI3","ABTB2","HMOX1",
  "RNASE7","KRT78","NABP1","MYOZ1","ADGRF1","ATP6V1C2","DUSP5","GPAT3","ERO1A","B3GNT6",
  "LPIN1","TIMP2","PSCA","KRT17","PADI2","S100P","SYNPO2L","ALDH1A3","LINC01559","FAM25A",
  "GRPEL2","WFDC12","TMC1","SPINK7","STEAP4","YOD1","C15orf48","LINC01214","C18orf25","LAYN",
  "MUCL3","GADD45B","LMO7","ECM1","AL133453.1","TMPRSS11B","AIF1L","ZNF426","SPAG1","CRCT1",
  "CGNL1","BTBD19","MUC1","ELF3","PADI1","ABLIM3","MAFF","SEMA7A","B3GALT5","MXD1",
  "QSOX1","LINC02099","SCARB1","EPS8L1","LNX1","NDRG2","KIFC3","IL18","RND3","KRT7",
  "TP53INP2","SULT1B1","SORT1","PPP1R3C","PNLIPRP3","TMPRSS11D","FCGBP","NCCRP1","TMPRSS11E","TJP3",
  "FAM83D","SNX9","HCG22","ZNF431","RSAD2","ADAM9","TP53BP2","BPGM","ANKRD37","TMPRSS2"
)

mod19 <- c(
  "IL24","S100A7A","GP2","KRT75","TLL2","SNTG2","AMTN","ORM1","SLC26A4","AC020741.1",
  "MT4","ICOS","ITK","CAMK2A","AC084816.1","DCLK1","NTNG1","AC093515.1","PDE10A","SULT1E1",
  "LY86-AS1","P2RY8","NTS","LOXL4","IL6","CALB1","DEFB4A","SLAMF1","AL110292.1","KCNIP4",
  "IFNG-AS1","CYP26A1","SLC38A11","PTPRH","NPSR1","CNBD1","CAPN8","AC004704.1","EPHA7","C17orf77",
  "LINC02484","GRM8","AC009088.3","DCC","GATA2-AS1","TRPM3","RGS3","CXCL1","AC073365.1","AC097460.1",
  "AADAC","LGR5","AP000561.1","LINC02432","CLEC2A","LINC00871","STXBP5L","COL26A1","NRG3","AC243829.4",
  "CTNND2","NCAM2","AC068580.3","XIST","CLEC2L","VTCN1","PTPRT","LINC01250","P2RY6","ENKUR",
  "CCL18","FREM2","ZNF114","B3GALT5","TRBV23-1","NKAIN3","INHBA","MYO3A","AL157871.4","XKR4",
  "IL1B","TRIM55","PAX7","AL138720.1","NECAB2","NEFL","MIR548XHG","LINC01194","SYN2","GDA",
  "SHISA9","ANKRD30A","ZNF831","FGF19","AL035401.1","LINC02099","FABP4","LINC01456","HSD11B1","GLIS1"
)

mod20 <- c(
  "MYBPC1","GRM8","PLEKHS1","KRT75","PTPRH","SLC38A11","XKR4","CTNND2","TIGIT","MYO3A",
  "ZNF114","DCC","IL24","SERPINE1","PKHD1","DEFB4A","LGR5","SLC26A9","LINC00871","ANKRD30A",
  "VIM-AS1","ENKUR","AC092164.1","SPIB","FANCA","INHBA","UHRF1","C17orf77","CNBD1","CLEC2L",
  "TLL2","LINC01250","SLC26A4","LINC02432","SNTG2","AC073365.1","ZNF831","EGR2","P2RY8","DTL",
  "CXCL9","ORM1","CXCL1","S100A7A","AC097460.1","COL26A1","LINC01194","AC009081.1","C11orf16","CSF3",
  "IFNG-AS1","AP000561.1","C10orf67","AC004704.1","NPSR1","FLT3","AL110292.1","KCNIP4","NRG3","RDH10-AS1",
  "PTPRT","LINC02484","VTCN1","CCL18","AP001636.3","NECAB2","AL035401.1","KIF14","XIST","SULT1E1",
  "CLEC2A","AC093515.1","AC084816.1","TRBV23-1","IL6","AC068580.3","AL157871.4","NKAIN3","ADAMDEC1","CD163L1",
  "MT4","AL138720.1","KCND2","SHISA9","TRIM55","FREM2","NEFL","MIR548XHG","AL162414.1","GATA2-AS1"
)

# module score on scRNA-seq
sq.sc <- AddModuleScore(sq.sc, features = list(mod20), name = "mod20_score")

p_umap <- FeaturePlot(
  sq.sc, reduction = "umap",
  features = "mod20_score1",
  cols = c("lightgray", "red3"),
  order = TRUE, pt.size = 1
) & theme(text = element_text(family = "Arial"))

p_vln <- VlnPlot(
  sq.sc,
  features = "mod20_score1",
  group.by = "Annotation_v3",
  pt.size  = 0
) + theme(axis.text.x = element_text(size = 13, angle = 45))

cairo_pdf("figure/Fig3/escc_mod20_module_score_umap_vln.pdf",
          width = 11, height = 4.8, family = "Arial")
plot_grid(p_umap, p_vln, rel_widths = c(1, 1.2))
dev.off()

# module score on Visium
cairo_pdf("figure/Fig3/escc_mod20_spatial.pdf",
          width = 18, height = 6, family = "Arial")
SpatialFeaturePlot(merged, features = "Suprabasal", ncol = 5) &
  theme(
    aspect.ratio  = 1,
    legend.text   = element_text(size = 5),
    legend.title  = element_text(size = 7),
    plot.title    = element_text(size = 9)
  )
dev.off()

# Fig S3F ----

library(Seurat)
library(Matrix)
library(RColorBrewer)

# aggregate epithelial subtypes into single "Epi" score
c2l_data <- GetAssayData(merged, assay = "C2L", slot = "data")

epi_agg <- Matrix::colSums(c2l_data[c("Tumor", "Basal", "Suprabasal"), , drop = FALSE])
epi_row <- Matrix(epi_agg, nrow = 1, sparse = TRUE)
rownames(epi_row) <- "Epi"
colnames(epi_row) <- colnames(c2l_data)

c2l_data <- rbind(c2l_data, epi_row)
c2l_data <- c2l_data[setdiff(rownames(c2l_data), c("Tumor", "Basal", "Suprabasal")), , drop = FALSE]
merged[["Epi"]] <- CreateAssayObject(counts = c2l_data)

# normalize & scale
obj <- merged
obj <- NormalizeData(obj)
obj <- ScaleData(obj)

# mask spots outside epithelial niches
bad_niches <- c("niche2", "niche3", "niche6", "niche7",
                "niche9", "niche10", "niche11", "niche12")
bad_idx <- which(obj$niche %in% bad_niches)

scale_mat <- as.matrix(GetAssayData(obj[["C2F_sq"]], slot = "scale.data"))
scale_mat[, bad_idx] <- NA
obj[["C2F_masked"]] <- SetAssayData(obj[["C2F_sq"]], slot = "scale.data", new.data = scale_mat)
DefaultAssay(obj) <- "C2F_masked"

# Module2 - Module0 difference score
mat <- GetAssayData(obj, slot = "scale.data")[c("Module.2", "Module.0"), ]
obj$diff_Mod2_Mod0 <- as.numeric(mat["Module.2", ] - mat["Module.0", ])

# spatial plot
SpatialColors <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

p <- SpatialFeaturePlot(
  obj,
  features = "diff_Mod2_Mod0",
  images = paste0("s", 1:10),
  ncol = 1, pt.size.factor = 1.8
) &
  scale_fill_gradientn(colours = SpatialColors(100), na.value = "black") &
  theme(plot.title = element_blank(),
        plot.margin = margin(t = 37, r = 0, b = 37, l = 0),
        text = element_text(family = "Arial"),
        legend.text = element_text(size = 6)) &
  NoLegend()

cairo_pdf("figure/Fig3/escc_Module2_minus_Module0_spatial.pdf",
          width = 5, height = 50, family = "Arial")
print(p)
dev.off()


# Fig S3G ----

library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(CellChat)
library(fgsea)
library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/cellchat/TIC/")

date = "20250813"
project_name = "escc"

size_pt <- 7
escc <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250731_escc_global_final_annotation_hj.rds")
seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/object_share/20250806_escc_squamouscells_umap_tuning_sk_hj.rds")

escc@meta.data <- escc@meta.data %>% dplyr::mutate(Annotation_v3 = case_when(cellname %in% b ~ "Basal",
                                                                             cellname %in% s ~ "Suprabasal",
                                                                             cellname %in% ts ~ "Transitional suprabasal",
                                                                             cellname %in% tic ~ "Tumor initiating cell",
                                                                             cellname %in% lgr5 ~ "LGR5+ cancer stem cell",
                                                                             cellname %in% tumor ~ "Tumor",
                                                                             T ~ Annotation_v2))

mycells <- unique(escc$Annotation_v3)[!unique(escc$Annotation_v3) %in% c("Mito-high CD8", "Mito-high CD4", ".")]
escc <- subset(escc, subset = Annotation_v3 %in% mycells)

escc$Annotation_v3 <- factor(escc$Annotation_v3)
escc$meta <- escc$Annotation_v3

library(Seurat)
old <- GeomViolin$default_aes$linewidth
update_geom_defaults("violin", list(linewidth = 0.3))

fileidentity <- "vlnplot_TNF"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(18), height = cm_to_inch(6), family="Arial")
VlnPlot(escc, group.by = "Annotation_v3", features = "TNF", pt.size = 0) & theme_classic(base_family = "Arial") & theme(axis.text = element_text(size = 5), axis.title = element_text(size = 5), axis.text.x = element_text(hjust = 1, angle = 45), axis.ticks.x = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 6), line = element_line(linewidth = 0.3)) & NoLegend()
dev.off()

escc_subset <- subset(escc, subset = Annotation_v3 %in% c("Basal", "Suprabasal", "Transitional suprabasal", "Tumor initiating cell", "Tumor"))
escc_subset$Annotation_v3 <- factor(escc_subset$Annotation_v3, levels = c("Basal", "Suprabasal", "Transitional suprabasal", "Tumor initiating cell", "Tumor"))

fileidentity <- "vlnplot_TNFRSF1A"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(4), height = cm_to_inch(6), family="Arial")
VlnPlot(escc_subset, group.by = "Annotation_v3", cols =  c("#F8766D", "#0CB702", "blue4", "magenta2", "gold2", "#619CFF"), features = "TNFRSF1A", pt.size = 0) & theme_classic(base_family = "Arial") & theme(aspect.ratio = 0.7, axis.text = element_text(size = 5), axis.title = element_text(size = 5), axis.text.x = element_text(hjust = 1, angle = 45), axis.ticks.x = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 6), line = element_line(linewidth = 0.3)) & NoLegend()
dev.off()

# Fig S3I ----

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

slide_tmp <- "s1"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p1 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

slide_tmp <- "s2"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p2 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

slide_tmp <- "s3"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p3 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

slide_tmp <- "s4"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p4 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

slide_tmp <- "s5"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p5 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

slide_tmp <- "s6"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p6 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

slide_tmp <- "s7"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p7 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

slide_tmp <- "s8"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p8 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

slide_tmp <- "s9"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p9 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

slide_tmp <- "s10"
q <- seurat_obj@meta.data %>% dplyr::filter(orig.ident == slide_tmp) %>% dplyr::filter(diff > 0) %>% dplyr::pull(diff) %>% quantile(na.rm = T)
p10 <- SpatialFeaturePlot(seurat_obj, features = c("diff"), images = slide_tmp) + scale_fill_gradient2(low = "green3", mid = "ivory", high = "red3", midpoint = 0, limits = c(-1*q[4], q[4]*(1)), oob = scales::squish, na.value = "black")

library(patchwork)
fileidentity <- "monocyte_fraction_diff_v2"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4*5, height = 8, family = "Arial")
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+plot_layout(ncol = 5) & theme(legend.ticks = element_blank(), legend.text = element_blank())
dev.off()

# Fig S3J ----

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

#wilres_nonclassical <- readRDS("20250823_escc_nonclasscial_test_hj.rds")
totalres_nonclassical <- readRDS("20250823_escc_nonclasscial_total_hj.rds")

#wilres_nonclassical %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FC > 0)

totalres$slide <- factor(totalres$slide, levels = paste0("s", 1:10))
totalres %>% head()

library(scales)
#totalres_nonclassical$variable <- factor(totalres$variable, levels = c("non_classical_monocyte", "classical_monocyte"))

totalres_nonclassical$slide <- factor(totalres_nonclassical$slide, levels = paste0("s", 1:10))
totalres_nonclassical %>% ggplot(aes(x = slide, y = value, fill = niche)) +
  geom_boxplot(outliers = F) +
  theme_classic(base_family = "Arial") +
  ylab("log10(Fraction for each spot") +
  xlab("Slide") +
  theme(aspect.ratio = 0.5, axis.text = element_text(size = 5, colour = "black"), line = element_line(linewidth = 0.2), axis.title = element_text(size = 5),
        legend.key.height = unit(0.2,  "cm"), legend.key.width = unit(0.2,  "cm"), legend.title = element_text(size = 3), legend.text = element_text(size = 3)    ,
        axis.ticks.length = unit(0.05, "cm")) +
  scale_y_log10(labels = label_number()) +
  scale_fill_manual(values = c("red3", "blue3"))


