library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(CellChat)
library(showtext)
library(dittoSeq)
library(MAST)

library(showtext)
library(MAST)

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

size_pt <- 7

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/v2/Fig2/")

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")

# Fig2A ----

p1 <- DimPlot(seurat_obj, group.by = "Annotation_v2") +
  theme_set(theme_classic(base_family = "Arial")) +
  theme_void() +
  theme(aspect.ratio = 1, plot.title = element_blank(), legend.text = element_text(size = size_pt-3))

print(p1)

# Fig2B ----

library(Seurat)
library(dplyr)
library(showtext)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54


date = "20250816"
project_name = "escc"

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/step02_subclustering/squamous_cells/MAST_FIXED_subcluster/")

# volcano ----
# library(EnhancedVolcano)
# res <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/analysis/step02_subclustering/squamous_cells/MAST_FIXED_subcluster/test_res_basal_hj.rds")
# res$gene <- rownames(res)
# 
# posgenes <- res %>% dplyr::filter(avg_log2FC > 1) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(p_val_adj) %>% head(n=10) %>% pull(gene)
# neggenes <- res %>% dplyr::filter(avg_log2FC < -1) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(p_val_adj) %>% head(n=10) %>% pull(gene)
# 
# posgenes_basal <- posgenes
# 
# fileidentity <- "MAST_FIXED_basal"
# cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 7, height = 7, family = "Arial")
# EnhancedVolcano(res, x= 'avg_log2FC', y='p_val_adj', lab = res$gene, pCutoff = 0.05, FCcutoff = 1,
#                 selectLab = c(posgenes, neggenes), drawConnectors = TRUE,
#                 max.overlaps = Inf,
#                 title = "Basal vs Suprabasal+Tumor", subtitle = "") +
#   theme(text = element_text(family = "Arial"), aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
# dev.off()
# 
# library(EnhancedVolcano)
# res <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/analysis/step02_subclustering/squamous_cells/MAST_FIXED_subcluster/test_res_suprabasal_hj.rds")
# res$gene <- rownames(res)
# 
# posgenes <- res %>% dplyr::filter(avg_log2FC > 1) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(p_val_adj) %>% head(n=5) %>% pull(gene)
# neggenes <- res %>% dplyr::filter(avg_log2FC < -1) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(p_val_adj) %>% head(n=5) %>% pull(gene)
# 
# posgenes_suprabasal <- posgenes
# 
# fileidentity <- "MAST_FIXED_suprabasal"
# cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 7, height = 7, family = "Arial")
# EnhancedVolcano(res, x= 'avg_log2FC', y='p_val_adj', lab = res$gene, pCutoff = 0.05, FCcutoff = 1,
#                 selectLab = c(posgenes, neggenes), drawConnectors = TRUE,
#                 max.overlaps = Inf,
#                 title = "Suprabasal vs Basal+Tumor", subtitle = "") +
#   theme(text = element_text(family = "Arial"), aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
# dev.off()
# 
# library(EnhancedVolcano)
# res <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/analysis/step02_subclustering/squamous_cells/MAST_FIXED_subcluster/test_res_tumor_hj.rds")
# res$gene <- rownames(res)
# 
# posgenes <- res %>% dplyr::filter(avg_log2FC > 1) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(p_val_adj) %>% head(n=10) %>% pull(gene)
# neggenes <- res %>% dplyr::filter(avg_log2FC < -1) %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::arrange(p_val_adj) %>% head(n=10) %>% pull(gene)
# 
# posgenes_tumor <- posgenes
# 
# fileidentity <- "MAST_FIXED_tumor"
# cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 7, height = 7, family = "Arial")
# EnhancedVolcano(res, x= 'avg_log2FC', y='p_val_adj', lab = res$gene, pCutoff = 0.05, FCcutoff = 1,
#                 selectLab = c(posgenes, neggenes), drawConnectors = TRUE,
#                 max.overlaps = Inf,
#                 labSize = 4 ,
#                 title = "Tumor vs Basal+Suprabasal", subtitle = "") +
#   theme(text = element_text(family = "Arial"), aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
#   ylim(0, 400)
# dev.off()
# 
# seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")
# seurat_obj$Annotation_v2 <- factor(seurat_obj$Annotation_v2, levels = rev(c('Basal','Suprabasal','Tumor')))

library(dittoSeq)
size_pt <- 7
p <- dittoDotPlot(seurat_obj, vars = c(posgenes_basal, posgenes_suprabasal, posgenes_tumor), group.by = "Annotation_v2", size = size_pt) + 
  theme_bw(base_family = "Arial") +
  theme(line = element_line(linewidth = 0.3),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 45, size = size_pt-2, vjust = 0.9, hjust = 0.95),
        legend.text = element_text(size = size_pt-3), legend.position = "right", 
        legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.1, "cm"), legend.spacing.x = unit(0.02, 'cm'),
        legend.spacing.y = unit(0.1, "cm"),
        legend.ticks = element_blank(),
        legend.title = element_text(size= size_pt-3))+
  scale_size_area(max_size = 2) +
  ylab("") +scale_color_gradientn(colours=c("blue","white", "red"), name = 'Average\nExpression',limits= c(-1.2, 1.2))


fileidentity <- "DEtest_dotplot"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(10), height = cm_to_inch(4), family = "Arial")
print(p)
dev.off()

# Fig2C ----

library(Seurat)
library(ggplot2)
library(SeuratObject)
library(patchwork)
library(fgsea)
library(dplyr)
library(fgsea)
library(showtext)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig2/visium_b_s_t/")

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

date = "20250818"
project_name <- "escc"

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")
fileidentity <- "visium_sq"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4*5, height = 8*3, family = "Arial")
SpatialFeaturePlot(seurat_obj, features = c("Basal", "Suprabasal", "Tumor"), ncol = 5)
dev.off()

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")
fileidentity <- "visium_sq_main"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4, height = 4*3, family = "Arial")
SpatialFeaturePlot(seurat_obj, features = c("Basal", "Suprabasal", "Tumor"), ncol = 1, images = "s4") & NoLegend()
dev.off()

p <- SpatialFeaturePlot(seurat_obj, features = c("Basal"), ncol = 1, images = "s4") & NoLegend()
p2 <- SpatialFeaturePlot(seurat_obj, features = c("Suprabasal"), ncol = 1, images = "s4") & NoLegend()
p3 <- SpatialFeaturePlot(seurat_obj, features = c("Tumor"), ncol = 1, images = "s4") & NoLegend()


dev.off()

# Fig2D ----

#Fig.2D####
setwd('~/Dropbox/ESCC_snrna/202507/infercnv/local/v1/')
library(Seurat)
library(infercnv)
library(dplyr)

seurat_obj <- readRDS('epi0.33_mono_infercnv_input.rds') #kakaocloud: /home/rocky/escc/infercnv/km/v5/20250806_escc_infercnv_km.R

mycounts <- as.matrix(seurat_obj@assays$RNA@counts)

#annotations_file
Mono_Epi_Annotations <- seurat_obj@meta.data %>% select(Annotation_v2)
Mono_Epi_Annotations[,2] <- row.names(Mono_Epi_Annotations)
Mono_Epi_Annotations <- Mono_Epi_Annotations[,c("V2", "Annotation_v2")]
colnames(Mono_Epi_Annotations) <- c("Cellname", "Annotation_v2")

write.table(Mono_Epi_Annotations, file = "20250815_escc_infercnv_monocyte_epi_km.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

library(Seurat)

Mono_highlight_Epi_infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = mycounts,
                                                        annotations_file = "20250815_escc_infercnv_monocyte_epi_km.txt",
                                                        delim = "\t",
                                                        gene_order_file = "~/Dropbox/database/infercnv/hg38_gencode_v27.txt",
                                                        ref_group_names = c("Classical Monocyte"))


Mono_highlight_Epi_infercnv_obj <- infercnv::run(Mono_highlight_Epi_infercnv_obj,
                                                 cutoff = 0.1,
                                                 out_dir = "infercnvR_Epi_Mono_filtered",
                                                 cluster_by_groups = T,
                                                 plot_steps = F,
                                                 denoise = T,
                                                 sd_amplifier = 1.5,
                                                 HMM = F,
                                                 output_format = 'pdf')

##pathology annotation bar####
setwd('~/Dropbox/ESCC_snrna/202507/')
obs <- read.csv('infercnv/local/v1/infercnvR_Epi_Mono_filtered/infercnv.observations.txt', sep = ' ')
epi <- readRDS('R_objects/20250730_escc_squamouscells_only_hj.rds')
tumor_cellid <- epi@meta.data %>% filter(Annotation_v2 == 'Tumor') %>% row.names() %>% gsub('-','.',.,fixed = T)
obs_tumor_cellid <- intersect(colnames(obs),tumor_cellid)
obs_tumor <- t(obs['MYC',obs_tumor_cellid])

epi@meta.data <- epi@meta.data %>% mutate(pathology_v2 = case_when(
  pathology == 'Normal' ~ 'Normal',
  pathology == 'Dysplasia' ~ 'Dysplasia',
  pathology %in% c('Microinvasive','Macroinvasive') ~ 'Carcinoma',
  TRUE ~ 'Others'
))

anno_pathology <- epi@meta.data %>% filter(rownames(.) %in% gsub('.','-',obs_tumor_cellid,fixed = T)) %>% select(pathology_v2)
row.names(anno_pathology) <- gsub('-','.',row.names(anno_pathology),fixed = T)
pdf('infercnv/local/v1/downstream/20250815_escc_scrna_pathology_annotation_bar_km.pdf', width = 7, height = 7)
pheatmap(
  obs_tumor,
  annotation_row = anno_pathology,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_colors = list(
    pathology_v2 = c("Normal" = "#83f52c", 
                     "Dysplasia" = "#f3f315", 
                     "Carcinoma" = "#B22222")
  )
)
dev.off()

##Tumor burden test between pathologies within tumor####
references <- read.table("infercnv/local/v1/infercnvR_Epi_Mono_filtered/infercnv.references.txt",sep=" ",header=TRUE)
Mono_highlight_Epi_infercnv_obj <- readRDS("infercnv/local/v1/infercnvR_Epi_Mono_filtered/Mono_highlight_Epi_infercnv_obj.rds")

#pathology별 cell ID를 뽑는다.
Mono_Epi_Annotations2 <- read.table(file = "infercnv/local/v1/20250815_escc_infercnv_monocyte_epi_km.txt", sep = "\t")

Tumor_id <- Mono_Epi_Annotations2 %>% filter(V2 == 'Tumor') %>% pull(V1)
normal <- epi@meta.data %>% filter(row.names(.) %in% Tumor_id, pathology == 'Normal') %>% row.names()
dysplasia <- epi@meta.data %>% filter(row.names(.) %in% Tumor_id, pathology == 'Dysplasia') %>% row.names()
carcinoma <- epi@meta.data %>% filter(row.names(.) %in% Tumor_id, pathology %in% c('Microinvasive','Macroinvasive')) %>% row.names()

#infercnv matrix를 pathology별로 나눈다.
normal_mat <- Mono_highlight_Epi_infercnv_obj@expr.data[,normal]
dysplasia_mat <- Mono_highlight_Epi_infercnv_obj@expr.data[,dysplasia]
carcinoma_mat <- Mono_highlight_Epi_infercnv_obj@expr.data[,carcinoma]

#pathology별 cell들에 대해 cnv score가 0.99미만, 1.01초과인 cell들을 뽑는다.
normal_cnv = c()
for (i in 1:ncol(normal_mat)) {
  count = length(c(which(normal_mat[,i] < 0.99), which(normal_mat[,i] > 1.01)))
  normal_cnv = c(normal_cnv, count)
}

dysplasia_cnv = c()
for (i in 1:ncol(dysplasia_mat)) {
  count = length(c(which(dysplasia_mat[,i] < 0.99), which(dysplasia_mat[,i] > 1.01)))
  dysplasia_cnv = c(dysplasia_cnv, count)
}

carcinoma_cnv = c()
for (i in 1:ncol(carcinoma_mat)) {
  count = length(c(which(carcinoma_mat[,i] < 0.99), which(carcinoma_mat[,i] > 1.01)))
  carcinoma_cnv = c(carcinoma_cnv, count)
}

wilcox.test(normal_cnv, dysplasia_cnv)$p.value #0.1706954
wilcox.test(normal_cnv, carcinoma_cnv)$p.value #0.003322012
wilcox.test(dysplasia_cnv, carcinoma_cnv)$p.value #6.548891e-57
wilcox.test(c(normal_cnv,dysplasia_cnv), carcinoma_cnv)$p.value #6.180953e-57


# Fig2E ----

#Fig.2E####
## tumor clone diversity####
setwd('~/Dropbox/ESCC_snrna/202507/infercnv/local/v1/')
library(Seurat)
library(dplyr)
library(ape)


tree <- read.tree("downstream/infercnv.tumor_dendrogram.txt")
hc <- as.hclust.phylo(tree)
k <- 8
clusters <- cutree(hc, k = k)
cell_clusters <- split(names(clusters), clusters)
str(cell_clusters)

cell_clusters

library(ggplot2)
library(tidyverse)
library(dplyr)
library(Seurat)
tumor_k1 <- read.table("downstream/infercnv.observation_groupings_tumor_c1.txt")
tumor_k2 <- read.table("downstream/infercnv.observation_groupings_tumor_c2.txt")
tumor_k3 <- read.table("downstream/infercnv.observation_groupings_tumor_c3.txt")
tumor_k4 <- read.table("downstream/infercnv.observation_groupings_tumor_c4.txt")
tumor_k5 <- read.table("downstream/infercnv.observation_groupings_tumor_c5.txt")
tumor_k6 <- read.table("downstream/infercnv.observation_groupings_tumor_c6.txt")
tumor_k7 <- read.table("downstream/infercnv.observation_groupings_tumor_c7.txt")
tumor_k8 <- read.table("downstream/infercnv.observation_groupings_tumor_c8.txt")

tumor_k1[,2] <- "Clone 1"
tumor_k2[,2] <- "Clone 2"
tumor_k3[,2] <- "Clone 3"
tumor_k4[,2] <- "Clone 4"
tumor_k5[,2] <- "Clone 5"
tumor_k6[,2] <- "Clone 6"
tumor_k7[,2] <- "Clone 7"
tumor_k8[,2] <- "Clone 8"

total <- rbind(tumor_k1, tumor_k2, tumor_k3, tumor_k4, tumor_k5, tumor_k6, tumor_k7, tumor_k8)
colnames(total) <- c("cell", "dendrogram")
row.names(total) <- total$cell
total[,1] <- NULL

epi <- readRDS('epi0.33_mono_infercnv_input.rds')
Idents(epi) <- 'Annotation_v2'
tumor <- subset(epi, subset = Annotation_v2 == 'Tumor')
tumor <- AddMetaData(tumor, total)
tumor <- SetIdent(tumor, value="dendrogram")
tumor_clones_per_donor <- tumor@meta.data %>% select(donor, dendrogram) %>% group_by(donor, dendrogram) %>% summarise(n = n()) %>% as.data.frame()
tumor_clones_per_donor <- tumor_clones_per_donor %>% group_by(donor) %>% mutate(prop = n / sum(n)) %>% as.data.frame()
tumor_clones_per_donor$prop <- round(tumor_clones_per_donor$prop, 3)
tumor_clones_per_donor <- plyr::ddply(tumor_clones_per_donor, "donor", transform, label_y=cumsum(prop))
tumor_clones_per_donor$dendrogram <- factor(tumor_clones_per_donor$dendrogram, 
                                            levels = c("Clone 1","Clone 2","Clone 3","Clone 4","Clone 5","Clone 6","Clone 7", "Clone 8"))
library(openxlsx)
write.xlsx(tumor_clones_per_donor, file = 'downstream/shannon/20250905_escc_tumor_clones_per_donor.xlsx')
tumor_clones_per_donor <- read.xlsx('downstream/shannon/20250905_escc_tumor_clones_per_donor.xlsx')

library(vegan)
shannon <- c()
patient <- c()
for (i in unique(tumor_clones_per_donor$donor)) {
  shannon <- c(shannon, round(diversity(tumor_clones_per_donor[tumor_clones_per_donor$donor == i,]$n), 3))
  patient <- c(patient, i)
}
shannon_patient <- data.frame(patient, shannon)
shannon_order <- order(shannon, decreasing = T)
shannon_patient$patient[shannon_order]
shannon_patient$patient <- factor(shannon_patient$patient, levels = c(shannon_patient$patient[shannon_order]))
shannon_patient

write.xlsx(shannon_patient, file = 'downstream/shannon/20250905_escc_tumor_clones_per_donor_shannon_diversity.xlsx')
shannon_patient <- read.xlsx('downstream/shannon/20250905_escc_tumor_clones_per_donor_shannon_diversity.xlsx')

shannon_patient <- shannon_patient %>% filter(patient %in% c('19','20','21','24','26','27','31','33'))

library(ggplot2)

#19: darkgreen
#20: green2
#21: #26A63A
#22: #B4B61A
#24: #704D9E
#26: #CF63A6
#27: #F7A086
#28: #F3E79A
#31: #273871
#33: blue
#36: #5087C1
#37: #ACCCE4

library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

cairo_pdf("downstream/shannon/20250914_escc_clones_per_donor_shannon_plot_km.pdf", width = 4, height = 7, family = "Arial")
ggplot(shannon_patient) + 
  geom_bar(stat = "identity", aes(x=shannon, y = reorder(patient, shannon), fill=patient), color = "black", size = 0.3) +
  theme_classic(base_family = 'Arial') + scale_fill_manual(values = c(c("darkgreen", "green2", "#26A63A", "#704D9E", "#CF63A6", "#F7A086", "#273871", "blue"))) +
  theme(axis.title.x = element_text(size=18, vjust = 0, family="Arial", margin = margin(0.28, 0, 0, 0, "cm")), 
        axis.text.x = element_text(family="Arial", size = 18, angle = 45, hjust = 1, vjust = 1, color = 'black'), 
        axis.title.y = element_text(size=18, vjust = 0, family="Arial", margin = margin(0, 0.28, 0, 0, "cm")), 
        axis.text.y = element_text(family="Arial", size = 18, hjust = 1, vjust = 1, color = 'black'), panel.grid = element_blank(),
        legend.position=c("bottom"),
        legend.title = element_blank(),
        legend.text = element_text(family="Arial", size=18),
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.key.size = unit(0.4, 'cm')) + 
  labs(x = "Shannon index", y = "Patient") + 
  geom_text(aes(x=shannon, y=patient, label=shannon), colour="black", family="Arial", size = 6, position = position_stack(vjust=0.5))+
  guides(fill = guide_legend(override.aes = list(colour = "black", size = 5), ncol = 4))
#scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4), labels= c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "1.4"))
dev.off()

# Extract Shannon diversity vector and name it with patient IDs
shannon_vector <- setNames(shannon_patient$shannon, shannon_patient$patient)
shannon_vector <- shannon_vector[c(1,2,4,3,7,5,8,6)]
# Compute distance matrix
dist_matrix <- dist(shannon_vector, method = "euclidean")

# Hierarchical clustering
hc <- hclust(dist_matrix, method = "average")  # you can also use "complete", "single", etc.

# Plot dendrogram
cairo_pdf("downstream/shannon/20250914_escc_clones_per_donor_shannon_plot_dendrogram_km.pdf", width = 6, height = 3, family = "Arial")
plot(hc, 
     main = "Dendrogram of Shannon Diversity Indices",
     xlab = "Patients",
     sub = "",
     ylab = "Distance",
     lwd = 2,
     ylim = c(0, 0.3))
dev.off()

# Fig2F ----

library(devtools)
install_github("lhe17/nebula")
library(nebula)
library(showtext)
library(MAST)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

setwd("/Users/hojin/Dropbox/project/ESCC/submit/analysis/NEBULA/")

# tumor vs suprabasal vs basal

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")

mycells1 <- readRDS("~/Dropbox/project/ESCC/submit/analysis/step02_subclustering/squamous_cells/MAST_FIXED_subcluster/test2_tumor_random_cells_hj.rds")
mycells2 <- readRDS("~/Dropbox/project/ESCC/submit/analysis/step02_subclustering/squamous_cells/MAST_FIXED_subcluster/test2_suprabasal_random_cells_hj.rds")
mycells3 <- readRDS("~/Dropbox/project/ESCC/submit/analysis/step02_subclustering/squamous_cells/MAST_FIXED_subcluster/test2_basal_random_cells_hj.rds")

mycells <- c(mycells1, mycells2, mycells3)

seurat_obj_subset <- subset(seurat_obj, cells = mycells)

length(mycells)

seurat_obj_subset@meta.data <- seurat_obj_subset@meta.data %>% dplyr::mutate(nebula_anno = case_when(Annotation_v2 == "Tumor" ~ "Tumor",
                                                                                                     T ~ "Others"))

seuratdata <- scToNeb(obj = seurat_obj_subset, assay = "RNA", id = "donor", pred = c("nebula_anno", "donor"), offset = "nCount_RNA")
df = model.matrix(~nebula_anno+donor, data = seuratdata$pred)
## include only the first two cell types in the model to avoid separation due to too many binary variables
data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoTumor")],offset=seuratdata$offset)
re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset)

saveRDS(re, "nebula_tumor_hj.rds")

seurat_obj_subset@meta.data <- seurat_obj_subset@meta.data %>% dplyr::mutate(nebula_anno = case_when(Annotation_v2 == "Suprabasal" ~ "Suprabasal",
                                                                                                     T ~ "Others"))

seuratdata <- scToNeb(obj = seurat_obj_subset, assay = "RNA", id = "donor", pred = c("nebula_anno", "donor"), offset = "nCount_RNA")
df = model.matrix(~nebula_anno+donor, data = seuratdata$pred)
colnames(df)
## include only the first two cell types in the model to avoid separation due to too many binary variables
data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoSuprabasal")],offset=seuratdata$offset)
re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset)

saveRDS(re, "nebula_suprabasal_hj.rds")

seurat_obj_subset@meta.data <- seurat_obj_subset@meta.data %>% dplyr::mutate(nebula_anno = case_when(Annotation_v2 == "Basal" ~ "Basal",
                                                                                                     T ~ "Others"))

seurat_obj_subset$nebula_anno <- factor(seurat_obj_subset$nebula_anno, levels = rev(c("Basal", "Others")))

seuratdata <- scToNeb(obj = seurat_obj_subset, assay = "RNA", id = "donor", pred = c("nebula_anno", "donor"), offset = "nCount_RNA")
df = model.matrix(~nebula_anno+donor, data = seuratdata$pred)
colnames(df)
## include only the first two cell types in the model to avoid separation due to too many binary variables
data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoBasal")],offset=seuratdata$offset)
re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset)

saveRDS(re, "nebula_basal_hj.rds")


## GSEA ##

library(fgsea)
library(tibble)

mygsea_res <- list()

# Basal
celltype <- "Basal"
re <- readRDS("nebula_basal_hj.rds")
res <- re$summary %>% dplyr::select(gene, logFC_nebula_annoBasal, se_nebula_annoBasal, p_nebula_annoBasal)
res$Padj <- p.adjust(res$p_nebula_annoBasal, method = "bonferroni")
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
fgseaRes_final_df$celltype <- celltype
mygsea_res[[celltype]] <- fgseaRes_final_df

# Suprabasal
celltype <- "Suprabasal"
re <- readRDS("nebula_suprabasal_hj.rds")
res <- re$summary %>% dplyr::select(gene, logFC_nebula_annoSuprabasal, se_nebula_annoSuprabasal, p_nebula_annoSuprabasal)
res$Padj <- p.adjust(res$p_nebula_annoSuprabasal, method = "bonferroni")
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
fgseaRes_final_df$celltype <- celltype
mygsea_res[[celltype]] <- fgseaRes_final_df

# Tumor
celltype <- "Tumor"
re <- readRDS("nebula_tumor_hj.rds")
res <- re$summary %>% dplyr::select(gene, logFC_nebula_annoTumor, se_nebula_annoTumor, p_nebula_annoTumor)
res$Padj <- p.adjust(res$p_nebula_annoTumor, method = "bonferroni")
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
fgseaRes_final_df$celltype <- celltype
mygsea_res[[celltype]] <- fgseaRes_final_df

mygsea_res <- do.call(rbind, mygsea_res)

mygsea_res$pathway <- gsub(mygsea_res$pathway, pattern = "HALLMARK_", replacement = "")
mygsea_res$pathway <- gsub(mygsea_res$pathway, pattern = "_", replacement = " ")
mygsea_res$pathway <- tolower(mygsea_res$pathway)
mygsea_res$pathway <- paste0(toupper(substr(mygsea_res$pathway, 1, 1)), substr(mygsea_res$pathway, 2, nchar(mygsea_res$pathway)))

test <- mygsea_res %>% dplyr::filter(padj < 0.05)
test$celltype <- factor(test$celltype, levels = c("Basal","Suprabasal","Tumor"))
test <- test %>% dplyr::arrange(celltype, padj)

myindex <- unique(test$pathway)

test$celltype <- factor(test$celltype, levels = rev(c("Tumor", "Suprabasal", "Basal")))

count_sort <- test %>% dplyr::group_by(pathway) %>% dplyr::summarise(n=n()) %>% dplyr::arrange(n)
test <- merge(count_sort, test, by = "pathway", all.y = T)
test$n <- factor(test$n, levels = c(1,2,3))
myindex <- test %>% dplyr::arrange(celltype, n, desc(padj)) %>% dplyr::pull(pathway) %>% unique()

test$pathway <- factor(test$pathway, levels = myindex)

test <- test %>% dplyr::arrange(pathway)

test <- test %>% dplyr::mutate(significance = case_when(padj < 0.05 ~ "significant",
                                                        T ~ "not significant"))

max_size_val <- max(-log10(test$padj), na.rm = TRUE)
min_size_val <- min(-log10(test$padj), na.rm = TRUE)

fileidentity <- "tumor_basal_suprabasal_gsea_dotplot"
date <- "20250930"
project_name <- "escc"

cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(6), height = cm_to_inch(7), family = "Arial")

test %>% ggplot(aes(x = celltype, y = pathway)) + 
  geom_point(aes(
    color = NES,
  ), na.rm = TRUE) + 
  geom_point(aes(color = NES, size = -log10(padj))) + ylab(NULL)+ theme_classic() + 
  theme_classic(base_family = "Arial") +
  theme(
    line = element_line(linewidth = 0.3),
    axis.text.x = element_text(hjust = 1,size = 5, angle = 30, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"), 
    axis.title.x = element_blank(), 
    legend.text = element_text(size = 4),
    legend.title = element_text(size = 4),
    legend.ticks = element_blank(),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width  = unit(0.1, "cm"),
    legend.spacing.x  = unit(0.02, "cm")) + 
  guides(size = guide_legend(title = expression(-log[10] ~ "(padj)"))) + 
  scale_size_area(max_size = 3) +
  scale_color_gradientn(colours=c("blue","white","red"))
dev.off()

# Fig2H ----

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

setwd(paste0("/Users/hojin/Dropbox/project/ESCC/submit/analysis/scenic/sq/"))

# 1. wilcoxon

wilcoxon_res <- list()
regulonAUC <- readRDS('/Users/hojin/Dropbox/project/ESCC/submit/analysis/scenic/sq/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")
metadat <- seurat_obj@meta.data
#scenicOptions <- readRDS("int/scenicOptions.Rds")
#regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonmat <- regulonAUC@assays@data$AUC
TFs <- rownames(regulonmat) %>% unique()

## 1.1 split cells

#regulonmat_sc <- regulonmat[, intersect(colnames(regulonmat), colnames(seurat_obj_subset1))]
res_df <- data.frame()
TFs <- rownames(regulonmat) %>% unique()
#TFs <- TFs[!str_detect(TFs, "_extended")]
metadat <- seurat_obj@meta.data
metadat$cellname <- rownames(metadat)

donors <- seurat_obj@meta.data$donor %>% unique()

res_df <- data.frame()

for (donor_tmp in donors) {
  
  cellnames <- seurat_obj@meta.data %>% dplyr::filter(donor == donor_tmp) %>% rownames()
  
  for (i in 1:length(TFs)) {
    tf_tmp <- TFs[i]
    react <- regulonmat[tf_tmp, cellnames]
    df_tmp <- data.frame(cellname = names(react), regulon_activity = react)
    df_tmp <- merge(df_tmp, metadat, by.x = "cellname", by.y = "cellname", all.x = T)
    
    for (celltype_tmp in c(c("Tumor", "Basal", "Suprabasal"))) {
      df1 <- df_tmp %>% dplyr::filter(Annotation_v2 == celltype_tmp)
      df2 <- df_tmp %>% dplyr::filter(Annotation_v2 != celltype_tmp)
      res <- wilcox.test(df1$regulon_activity, df2$regulon_activity)
      pval_tmp <- res$p.value
      res_tmp <- data.frame(TF = tf_tmp,
                            sample = donor_tmp,
                            group1 = celltype_tmp,
                            group1_mean = mean(df1$regulon_activity),
                            group2_mean = mean(df2$regulon_activity),
                            pval = pval_tmp)
      res_df <- rbind(res_df, res_tmp)
    }
  }
  
}

#res_df$padj <- p.adjust(res_df$pval, method = "bonferroni")

saveRDS(res_df, "b_sq_t_test_extended2_hj.rds")

res_df <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/analysis/scenic/sq/b_sq_t_test_extended2_hj.rds")

sig_res <- data.frame()
for (donor_tmp in donors) {
  
  res_df_tmp <- res_df %>% dplyr::filter(sample == donor_tmp) %>% dplyr::filter(group1 == "Tumor")
  res_df_tmp$padj <- p.adjust(res_df_tmp$pval, method = "bonferroni")
  res_df_tmp <- res_df_tmp %>% dplyr::filter(padj < 0.0000000001) %>% dplyr::mutate(log2FC = log2(group1_mean/group2_mean)) %>% dplyr::filter(group2_mean < group1_mean) %>% dplyr::filter(log2FC > 0.5)
  sig_res <- rbind(sig_res, res_df_tmp)
  
}

mytfs <- sig_res %>% dplyr::group_by(TF) %>% dplyr::summarise(n=n()) %>% dplyr::filter(n > 4) %>% dplyr::pull(TF)

library(scales)
regulonAUC_scaled <- list()

for (donor_tmp in donors) {
  
  seurat_obj_subset <- subset(seurat_obj, subset = donor == donor_tmp)
  regulonAUC_scaled_tmp <- getAUC(regulonAUC)[mytfs, colnames(seurat_obj_subset)]
  regulonAUC_scaled_tmp <- t(scale(t(regulonAUC_scaled_tmp), center = T, scale=T))
  regulonAUC_scaled[[donor_tmp]] <- regulonAUC_scaled_tmp
  
}

regulonAUC_scaled <- do.call(cbind, regulonAUC_scaled)

# plot
celltypes <- seurat_obj@meta.data %>% select(Annotation_v2)
celltypes$New_Identity <- factor(celltypes$Annotation_v2, levels = c("Basal", "Suprabasal", "Tumor"))
selectedResolution <- "New_Identity"
cellsPerTypes <- split(rownames(celltypes), celltypes[,selectedResolution]) 

set.seed(1234); cellorder <- c(
  unlist(lapply(split(cellsPerTypes$Basal, 
                      seurat_obj@meta.data[cellsPerTypes$Basal, "platform"]), sample)),
  unlist(lapply(split(cellsPerTypes$Suprabasal, 
                      seurat_obj@meta.data[cellsPerTypes$Suprabasal, "platform"]), sample)),
  unlist(lapply(split(cellsPerTypes$Tumor, 
                      seurat_obj@meta.data[cellsPerTypes$Tumor, "platform"]), sample))
)

celltypes <- c(rep(x = "Basal", times = length(cellsPerTypes$Basal)), rep("Suprabasal", times = length(cellsPerTypes$Suprabasal)), rep("Tumor", times = length(cellsPerTypes$Tumor)))
cell_type <- data.frame("cell type" = rep(c('Basal', 'Suprabasal','Tumor'), c(length(cellsPerTypes$Basal), length(cellsPerTypes$Suprabasal), length(cellsPerTypes$Tumor))))
regulonAUC_scaled_order <- regulonAUC_scaled[mytfs, cellorder]
row.names(cell_type) <- colnames(regulonAUC_scaled_order)
palette_length = 100

setwd(paste0("/Users/hojin/Dropbox/project/ESCC/submit/figures/v2/Fig2/"))

date <- "20250818"
project_name <- "escc"
cm_to_inch <- function(cm) cm / 2.54
fileidentity <- "Fig2F"

platform_col <- c("scRNA" = "grey", "snRNA" = "orange3")
platform_vec <- seurat_obj@meta.data[colnames(regulonAUC_scaled_order), "platform"] |> as.character()
platform_vec <- factor(platform_vec, levels = names(platform_col))

library(colorRamp2)
ha = HeatmapAnnotation(Platform = platform_vec,
                       Celltype = c(rep("Basal",length(cellsPerTypes$Basal)),
                                    rep("Suprabasal", length(cellsPerTypes$Suprabasal)),
                                    rep("Tumor", length(cellsPerTypes$Tumor))),
                       col = list(
                         Platform = platform_col,
                         Celltype = c("Basal" = "#F8766D", "Suprabasal" = "#0CB702", "Tumor" = "#619CFF")
                       ),
                       simple_anno_size = unit(0.3, "cm"),
                       show_legend = FALSE,
                       annotation_name_gp = gpar(fontsize = 7),
                       annotation_legend_param = list(
                         title_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                         labels_gp = gpar(fontsize = 5, fontfamily = "Arial")
                       ))

col_fun = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "red"))

png(filename = paste0(date, "_", project_name, "_", fileidentity, "_hj.png"),
    width = 5, height = 4, units = "in", res = 600, family = "Arial")
set.seed(1234);ht <- Heatmap(matrix = regulonAUC_scaled_order, name = "mat", col = col_fun,
                             cluster_columns = F,
                             cluster_column_slices = F,
                             clustering_method_rows = "ward.D2",
                             show_column_names = F,
                             use_raster = TRUE,
                             raster_device = "png",
                             raster_quality = 5,
                             raster_device_param = list(res = 600),
                             show_heatmap_legend = FALSE,
                             #show_annotation_legend = FALSE,
                             #width = unit(ncol(regulonAUC_scaled_order) * 0.000001, "in"),
                             row_names_gp = gpar(fontsize = 7, fontfamily = "Arial"),
                             column_names_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                             # heatmap_legend_param = list(
                             #   title = "Score",
                             #   
                             #   title_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                             #   labels_gp = gpar(fontsize = 5, fontfamily = "Arial")
                             # ),
                             top_annotation = ha)
draw(ht, show_heatmap_legend = FALSE, annotation_legend_list = NULL)
dev.off()


library(colorRamp2)
ha = HeatmapAnnotation(Platform = platform_vec,
                       Celltype = c(rep("Basal",length(cellsPerTypes$Basal)),
                                    rep("Suprabasal", length(cellsPerTypes$Suprabasal)),
                                    rep("Tumor", length(cellsPerTypes$Tumor))),
                       col = list(
                         Platform = platform_col,
                         Celltype = c("Basal" = "#F8766D", "Suprabasal" = "#0CB702", "Tumor" = "#619CFF")
                       ),
                       simple_anno_size = unit(0.3, "cm"),
                       #show_legend = FALSE,
                       annotation_name_gp = gpar(fontsize = 7),
                       annotation_legend_param = list(
                         title_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                         labels_gp = gpar(fontsize = 5, fontfamily = "Arial")
                       ))

col_fun = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "red"))

png(filename = paste0(date, "_", project_name, '_', "legend", "_", fileidentity, "_hj.png"),
    width = 6, height = 4, units = "in", res = 600, family = "Arial")
set.seed(1234);ht <- Heatmap(matrix = regulonAUC_scaled_order, name = "mat", col = col_fun,
                             cluster_columns = F,
                             cluster_column_slices = F,
                             clustering_method_rows = "ward.D2",
                             show_column_names = F,
                             use_raster = TRUE,
                             raster_device = "png",
                             raster_quality = 5,
                             raster_device_param = list(res = 600),
                             #show_heatmap_legend = FALSE,
                             #show_annotation_legend = FALSE,
                             #width = unit(ncol(regulonAUC_scaled_order) * 0.000001, "in"),
                             row_names_gp = gpar(fontsize = 7, fontfamily = "Arial"),
                             column_names_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                             # heatmap_legend_param = list(
                             #   title = "Score",
                             #   
                             #   title_gp = gpar(fontsize = 5, fontfamily = "Arial"),
                             #   labels_gp = gpar(fontsize = 5, fontfamily = "Arial")
                             # ),
                             top_annotation = ha)
draw(ht)
dev.off()


## ORA ##
setwd(paste0("/Users/hojin/Dropbox/project/ESCC/submit/analysis/scenic/sq/"))
library(org.Hs.eg.db)
library(msigdbr)
library(clusterProfiler)

##enricher
regulons <- readRDS('./2.6_regulons_asGeneSet.Rds')
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
hs <- org.Hs.eg.db
row_order_vec <- row_order(ht)
ordered_row_names <- rownames(regulonAUC_scaled_order)[row_order_vec]

tfs <- gsub(ordered_row_names, pattern = " \\([0-9]*g\\)", replacement = "")

saveRDS(tfs, "20250818_tfs_hj.rds")

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

WriteXLS(enricher_res, "20250818_escc_sq_targetgenes_hallmark_hj.xlsx")

enricher_res <- readxl::read_excel("20250818_escc_sq_targetgenes_hallmark_hj.xlsx")
enricher_res <- enricher_res %>% dplyr::filter(p.adjust < 0.05)

enricher_res <- enricher_res %>% dplyr::group_by(TFs) %>% dplyr::arrange(p.adjust, .by_group = T) %>% dplyr::slice(1)
WriteXLS(enricher_res, "20250818_escc_sq_targetgenes_hallmark_top1_hj.xlsx")

saveRDS(mytfs, "20250818_escc_sq_tumor_tfs_final_list_hj.rds")

# Fig2G ----
library(Seurat)
library(ggplot2)
library(dplyr)

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")

date <- "20250818"
project_name <- "escc"

hallmark <- gmtPathways("/Users/hojin/Dropbox/db/h.all.v2025.1.Hs.symbols.gmt")

seurat_obj@active.assay <- c("Spatial")

for (tmp in names(hallmark)) {
  mygenes <- hallmark[[tmp]]
  seurat_obj <- AddModuleScore(seurat_obj, features = list(mygenes), name = tmp)
}

fileidentity <- "singlecell_gsea_on_visium"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4*3, height = 8, family = "Arial")
SpatialFeaturePlot(seurat_obj, features = c("HALLMARK_MYC_TARGETS_V11", "HALLMARK_MTORC1_SIGNALING1"), ncol = 3, images = c("s3", "s7", "s10")) & theme(plot.title = element_blank()) & NoLegend()
dev.off()

fileidentity <- "singlecell_gsea_on_visium_myc1_total"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), height = 5*10, width = 5, family = "Arial")
p <- SpatialFeaturePlot(
  seurat_obj,
  features = "HALLMARK_MYC_TARGETS_V11",
  images = paste0("s", 1:10),
  ncol = 1
) & theme(plot.title = element_blank(), plot.margin = margin(t = 37, r = 0, b = 37, l = 0)) & NoLegend()
print(p)
dev.off()

fileidentity <- "singlecell_gsea_on_visium_mtorc1_total"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), height = 5*10, width = 5, family = "Arial")
p <- SpatialFeaturePlot(
  seurat_obj,
  features = "HALLMARK_MTORC1_SIGNALING1",
  images = paste0("s", 1:10),
  ncol = 1
) & theme(plot.title = element_blank(), plot.margin = margin(t = 37, r = 0, b = 37, l = 0)) & NoLegend()
print(p)
dev.off()

# Fig2I ----
#Fig.2I####
epi <- readRDS('R_objects/20250730_escc_squamouscells_only_hj.rds')
epi$pathology <- factor(epi$pathology, levels = c('Normal', 'Dysplasia', 'Microinvasive', 'Macroinvasive'))
DimPlot(epi, group.by = 'Annotation_v2', split.by = 'pathology')

epi <- DietSeurat(epi, dimreducs = c('pca','harmony','umap'))
epi <- subset(epi, Annotation_v2 == 'Tumor')
#all remaining genes should have a nonzero sum when summing across cells
counts <- GetAssayData(epi, slot = "counts")
gene_cell_counts <- Matrix::rowSums(counts > 0)
genes_expressed_in_10 <- names(gene_cell_counts[gene_cell_counts >= 10])
epi <- subset(epi, features = genes_expressed_in_10)

library(SeuratDisk)
SaveH5Seurat(epi, filename = "R_objects/20250730_escc_squamouscells_only_tumor_km.h5Seurat")
Convert("R_objects/20250730_escc_squamouscells_only_tumor_km.h5Seurat", dest = "h5ad")

#module scores
epi <- readRDS('R_objects/20250730_escc_squamouscells_only_hj.rds')
tumor_modules <- read.csv('hotspot/outs/20250813_escc_sc.snrna_tumor_hotspot_min_gene25_results_modules_km.csv')
tumor_modules <- na.omit(tumor_modules) %>% filter(Module != '-1')

for (i in unique(tumor_modules$Module)) {
  geneset <- tumor_modules %>% filter(Module == i) %>% pull(Gene)
  set.seed(20250825)
  epi <- AddModuleScore(epi, features = list(geneset), name=paste0('Module',i))
}
epi@meta.data %>% 
  write.xlsx(results_wide, "hotspot/outs/20250825_escc_sc.snrna_epi_hotspot_min_gene25_module_scores_wilcoxon_pathology_km.xlsx")

#Fig.2J####
library(enrichR)
library(openxlsx)
dbs <- c("MSigDB_Hallmark_2020")
df_enriched_flt_all <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Term","Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes", "Module"))
for(i in seq(1:9)) {
  genes <- tumor_modules %>% filter(Module == i) %>% pull(Gene)
  df_enriched <- as.data.frame(enrichr(genes, dbs))
  colnames(df_enriched) <- gsub("MSigDB_Hallmark_2020.", "", colnames(df_enriched))
  df_enriched$Genes <- gsub(";", ", ", df_enriched$Genes)
  write.xlsx(df_enriched, paste0("hotspot/outs/20250813_escc_sc.snrna_tumor_module",i,'_enrichR_hallmark_km.xlsx'))
  df_enriched_flt <- df_enriched %>% select(Term, Adjusted.P.value, Odds.Ratio, Combined.Score, Genes) %>% filter(Adjusted.P.value < 0.05) %>% arrange(desc(Combined.Score))
  if (nrow(df_enriched_flt) > 0) {
    df_enriched_flt$Module <- i
    df_enriched_flt_all <- rbind(df_enriched_flt_all,df_enriched_flt)
  }
}

write.xlsx(df_enriched_flt_all, "hotspot/outs/20250813_escc_sc.snrna_tumor_module_all_enrichR_hallmark_km.xlsx")

#hypoxia LEGs
VlnPlot(tumor, c('CA12', 'LDHA', 'EDN2', 'SERPINE1', 'BHLHE40', 'SLC2A1', 'MYH9', 'FOS', 'NDRG1', 'HK2', 'PFKP'), group.by = 'pathology', ncol = 3, pt.size = 0)

#dotplot
df_enriched_flt_all <- read.xlsx("hotspot/outs/20250813_escc_sc.snrna_tumor_module_all_enrichR_hallmark_km.xlsx")
df_enriched_flt_all$Module <- paste("Module", df_enriched_flt_all$Module)
df_enriched_flt_all$Module <- factor(df_enriched_flt_all$Module, levels = c('Module 2', 'Module 7', 'Module 1', 'Module 9'))
df_enriched_flt_all$Term <- factor(df_enriched_flt_all$Term, levels = unique(c(df_enriched_flt_all$Term[df_enriched_flt_all$Module == 'Module 2'],
                                                                               df_enriched_flt_all$Term[df_enriched_flt_all$Module == 'Module 7'],
                                                                               df_enriched_flt_all$Term[df_enriched_flt_all$Module == 'Module 1'],
                                                                               df_enriched_flt_all$Term[df_enriched_flt_all$Module == 'Module 9'])))

library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cairo_pdf("hotspot/outs/20250818_escc_sc.snrna_tumor_hotspot_min_gene25_module_scores_ORA_dotplot_figure_km.pdf", width = cm_to_inch(14), height = cm_to_inch(15), family = "Arial")
ggplot(df_enriched_flt_all, aes(Module,Term)) + geom_point(aes(color = Odds.Ratio, size = -log10(Adjusted.P.value))) + 
  ylab(NULL) + theme_minimal() +
  theme_bw(base_family = 'Arial') +
  theme(text = element_text(family = "Arial"),
        axis.text.x = element_text(colour = "black", angle = 45, size = 13, vjust=1, hjust=1, family = 'Arial'), 
        axis.text.y = element_text(colour = "black", size = 13, family = 'Arial'), plot.title = element_blank(), 
        legend.position = "right", legend.text = element_text(size = 13),
        legend.title = element_text(size = 13), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black", size=0.5),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.spacing = unit(1, 'cm'), legend.spacing.x = unit(0.5,'cm')) + 
  guides(size = guide_legend(title = expression(-log[10] ~ "(padj)"))) +
  #scale_size(range = c(2,6), breaks = c(10,20,30)) +
  scale_color_gradientn(colours=c("orange","red")) +
  scale_shape_manual(values = c(13,16))
dev.off()

