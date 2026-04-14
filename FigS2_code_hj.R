# Fig S2A ----

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/v2/Fig2/")

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")

p2 <- FeaturePlot(seurat_obj, features = c("COL17A1"), cols = c("grey", "red"), ncol = 1, order = F) + theme_set(theme_classic(base_family = "Arial")) + theme_void() + theme(aspect.ratio = 1,  plot.title = element_text(hjust = 0.5)) & NoLegend()
p3 <- FeaturePlot(seurat_obj, features = c("KRT13"), cols = c("grey", "red"), ncol = 1) + theme_set(theme_classic(base_family = "Arial")) + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) & NoLegend()
p4 <- FeaturePlot(seurat_obj, features = c("GNGT1"), cols = c("grey", "red"), ncol = 1) + theme_set(theme_classic(base_family = "Arial")) + theme_void() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) & NoLegend()

p2+p3+p4+plot_layout(ncol = 3)

# Fig S2B ----

# Fig S2C (KM) ----

# Fig S2D ----
#Fig.S2D####
setwd('~/Dropbox/ESCC_snrna/202507/infercnv/local/patient_pathology/')
seurat_obj <- readRDS('../../km/patient/epi_mono_highlight_infercnv_input_cells_km.rds')
seurat_obj$pathology <- factor(seurat_obj$pathology, levels = c('Normal','Dysplasia','Microinvasive','Macroinvasive'))
seurat_obj@meta.data %>% filter(Annotation_v2 == 'Tumor') %>% select(donor, pathology) %>% table
pt19_mic <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(19) & pathology == c("Microinvasive") | Annotation_v2 == c("Classical Monocyte") & donor == c(19))
pt19_mac <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(19) & pathology == c("Macroinvasive") | Annotation_v2 == c("Classical Monocyte") & donor == c(19))
pt20_dys <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(20) & pathology == c("Dysplasia") | Annotation_v2 == c("Classical Monocyte") & donor == c(20))
pt20_mic <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(20) & pathology == c("Microinvasive") | Annotation_v2 == c("Classical Monocyte") & donor == c(20))
pt20_mac <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(20) & pathology == c("Macroinvasive") | Annotation_v2 == c("Classical Monocyte") & donor == c(20))
pt21_dys <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(21) & pathology == c("Dysplasia") | Annotation_v2 == c("Classical Monocyte") & donor == c(21))
pt21_mac <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(21) & pathology == c("Macroinvasive") | Annotation_v2 == c("Classical Monocyte") & donor == c(21))
pt24_dys <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(24) & pathology == c("Dysplasia") | Annotation_v2 == c("Classical Monocyte") & donor == c(24))
pt24_mic <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(24) & pathology == c("Microinvasive") | Annotation_v2 == c("Classical Monocyte") & donor == c(24))
pt26_dys <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(26) & pathology == c("Dysplasia") | Annotation_v2 == c("Classical Monocyte") & donor == c(26))
pt27_nor <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(27) & pathology == c("Normal") | Annotation_v2 == c("Classical Monocyte") & donor == c(27))
pt27_dys <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(27) & pathology == c("Dysplasia") | Annotation_v2 == c("Classical Monocyte") & donor == c(27))
pt27_mac <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(27) & pathology == c("Macroinvasive") | Annotation_v2 == c("Classical Monocyte") & donor == c(27))
pt31_mic <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(31) & pathology == c("Microinvasive") | Annotation_v2 == c("Classical Monocyte") & donor == c(31))
pt33_dys <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(33) & pathology == c("Dysplasia") | Annotation_v2 == c("Classical Monocyte") & donor == c(33))
pt33_mic <- subset(seurat_obj, Annotation_v2 == c("Tumor") & donor == c(33) & pathology == c("Microinvasive") | Annotation_v2 == c("Classical Monocyte") & donor == c(33))


cases <- c("pt19_mic", "pt19_mac", "pt20_dys", "pt20_mic", "pt20_mac", "pt21_dys", "pt21_mac", "pt24_dys", "pt24_mic", "pt26_dys", "pt27_nor", "pt27_dys", "pt27_mac", "pt31_mic", "pt33_dys", "pt33_mic")
for (i in cases) {
  dir.create(paste0('~/Dropbox/ESCC_snrna/202507/infercnv/local/patient_pathology/',i))
  setwd(paste0('~/Dropbox/ESCC_snrna/202507/infercnv/local/patient_pathology/',i))
  seurat_obj <- get(i)
  mycounts <- as.matrix(seurat_obj@assays$RNA@counts)
  #annotations_file
  Mono_Epi_Annotations <- seurat_obj@meta.data %>% select(Annotation_v2)
  Mono_Epi_Annotations[,2] <- row.names(Mono_Epi_Annotations)
  Mono_Epi_Annotations <- Mono_Epi_Annotations[,c("V2", "Annotation_v2")]
  colnames(Mono_Epi_Annotations) <- c("Cellname", "Annotation_v2")
  write.table(Mono_Epi_Annotations, file = paste0("20250909_escc_infercnv_monocyte_epi_",i,"_km.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  pt_infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = mycounts,
                                          annotations_file = paste0("20250909_escc_infercnv_monocyte_epi_",i,"_km.txt"),
                                          delim = "\t",
                                          gene_order_file = "~/Dropbox/database/infercnv/hg38_gencode_v27.txt",
                                          ref_group_names = c("Classical Monocyte"))
  rm(Mono_Epi_Annotations,mycounts,seurat_obj)
  
  pt_infercnv_obj = infercnv::run(pt_infercnv_obj,
                                  cutoff=1,
                                  out_dir=paste0("output_",i),
                                  cluster_by_groups=FALSE,
                                  plot_steps=T,
                                  scale_data=T,
                                  denoise=T,
                                  noise_filter=0.12,
                                  analysis_mode='subclusters',
                                  HMM = T,
                                  HMM_type='i6',
                                  tumor_subcluster_partition_method = 'random_trees')
  rm(pt_infercnv_obj)
}


##Uphyloplot2####
#iMAC terminal
#pwd: /Users/yoo/Dropbox/ESCC_snrna/202507/infercnv/local/patient_pathology
for pt in pt19_mic pt19_mac pt20_dys pt20_mic pt20_mac pt21_dys pt21_mac pt24_dys pt24_mic pt26_dys pt27_nor pt27_dys pt27_mac pt31_mic pt33_dys pt33_mic; do
mkdir -p $pt/Uphyloplot2/Inputs
cp ~/Dropbox/database/uphyloplot2-2.3/uphyloplot2.py $pt/Uphyloplot2
sed '/^Classical Monocyte.Classical Monocyte/d' \
< $pt/output_$pt/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings \
> $pt/Uphyloplot2/Inputs/$pt.rand_trees.hmm_mode-subclusters.cell_groupings
done

for pt in pt19_mic pt19_mac pt20_dys pt20_mic pt20_mac pt21_dys pt21_mac pt24_dys pt24_mic pt26_dys pt27_nor pt27_dys pt27_mac pt31_mic pt33_dys pt33_mic; do
cd /Users/yoo/Dropbox/ESCC_snrna/202507/infercnv/local/patient_pathology/$pt/Uphyloplot2
python uphyloplot2.py
done

#changing percentage cutoff
for pt in pt19_mic pt19_mac pt20_dys pt20_mic pt20_mac pt21_dys pt21_mac pt24_dys pt24_mic pt26_dys pt27_nor pt27_dys pt27_mac pt31_mic pt33_dys pt33_mic; do
cd /Users/yoo/Dropbox/ESCC_snrna/202507/infercnv/local/patient_pathology/$pt/Uphyloplot2
python uphyloplot2.py -c 16
done


##correlation plot####
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpmisc)
library(ggpubr)

df_v1 <- cbind(c(rep('1',1), rep('2',6), rep('3',9)),
               c(1, 
                 2, 2, 2, 3, 2, 2, 
                 3, 2, 1, 2, 3,
                 2, 1, 2, 2)) %>% as.data.frame()
colnames(df_v1) <- c("stage", "subclones")
df_v1$stage <- factor(df_v1$stage, levels = c('1','2','3'))


library(showtext) 

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cairo_pdf("~/Dropbox/escc_figure_share/FigS2/20250914_escc_uphyloplot_subclone_by_pathology_km.pdf", width = 7, height = 4.5, family = "Arial")
ggplot(df_v1, aes(x=stage, y=subclones)) +
  geom_jitter(aes(color = stage, size = 3)) +
  xlab("") + ylab("# of subclones in tumor")+ scale_color_manual(values=c('#83f52c', '#f3f315', '#ff334d')) +  #'#83f52c', '#f3f315', '#ff6600', '#ff0099'
  geom_smooth(aes(x = as.numeric(stage), y = as.numeric(subclones)), inherit.aes = FALSE, color = 'red', method = 'lm') +
  #geom_smooth(aes(x = as.numeric(stage), y = as.numeric(subclones)),method = 'lm', color ='red', se = FALSE) +
  annotate("text", x = 1.5, y = 2, col = "black", size = 6.5,
           label = paste("r =", signif(cor.test(as.numeric(df_v1$stage), as.numeric(df_v1$subclones))$estimate,3), ", p =", signif(cor.test(as.numeric(df_v1$stage), as.numeric(df_v1$subclones))$p.value,3)),
           family = 'Arial') +
  theme_classic(base_family = 'Arial') +
  theme(legend.title = element_text(size = 18, family="Arial"),
        legend.text = element_text(size = 18, family="Arial"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        axis.title.x = element_text(size = 18, family = 'Arial', colour = 'black'),
        axis.text.x = element_text(size = 18, family = 'Arial', colour = 'black'),
        axis.title.y = element_text(size = 18, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 18, family = 'Arial', colour = 'black')) +
  scale_x_discrete(labels= c('Normal','Dysplasia', 'Carcinoma')) +
  guides(col=guide_legend(override.aes = list(size = 5)))
dev.off()

# Fig S2E ----

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

# Fig S2F ----

library(Seurat)
#library(mvtnorm)
library(ggplot2)
#library(QRM)
library(SeuratObject)
library(patchwork)
library(sctransform)
library(qusage)
#library(GSVA)
library(UCell)
library(fgsea)
library(dplyr)

library(fgsea)

library(showtext)

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig2/")

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

date = "20250811"
project_name <- "escc"

hallmark <- gmtPathways("/Users/hojin/Dropbox/db/h.all.v2025.1.Hs.symbols.gmt")

seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/20250806_escc_misty_dl.rds")

seurat_obj@active.assay <- c("Spatial")

for (tmp in names(hallmark)) {
  mygenes <- hallmark[[tmp]]
  seurat_obj <- AddModuleScore(seurat_obj, features = list(mygenes), name = tmp)
}

library(stringr)

mygeneset <- colnames(seurat_obj@meta.data)[str_detect(colnames(seurat_obj@meta.data), "HALLMARK_")]

SpatialDimPlot(seurat_obj, images = "s2", group.by = "niche")

df <- data.frame()
for (i in 1:length(mygeneset)) {
  group1 <- seurat_obj@meta.data %>% dplyr::filter(niche %in% c("niche8")) %>% dplyr::pull(mygeneset[i])
  group2 <- seurat_obj@meta.data %>% dplyr::filter(niche %in% c("niche1", "niche4")) %>% dplyr::pull(mygeneset[i])
  res <- wilcox.test(group1, group2)
  mypathway <- substring(mygeneset[i], 1, nchar(mygeneset[i]) - 1)
  df_tmp <- data.frame(pathway = mypathway, pvalue = res$p.value, group1_mean = mean(group1), group2_mean = mean(group2))
  df <- rbind(df, df_tmp)
}

df$Padj <- p.adjust(df$pvalue, method = "bonferroni")

#df <- df %>% dplyr::filter(Padj < 0.05)
myquant <- quantile(df$Padj)[3]
df <- df %>% dplyr::filter(Padj < myquant)
df <- df %>% dplyr::mutate(diff = group1_mean - group2_mean)
#df1<-df %>% dplyr::filter(diff > 0) %>% dplyr::arrange(Padj)
#df2<-df %>% dplyr::filter(diff < 0) %>% dplyr::arrange(Padj)

df <- df %>% dplyr::mutate(Category = case_when(diff > 0 ~ "Module score difference > 0",
                                                diff < 0 ~ "Module score difference < 0",
))


df <- df %>% dplyr::arrange(diff)
#df <- df %>% dplyr::group_by(color_group) %>% dplyr::arrange(desc(Padj), .by_group = T)

#df <- rbind(df2, df1)
#df$Padj[df$Padj == 0] = 1.244752e-304
first_word_cap <- function(x) {
  x <- tolower(x)  # 전체를 소문자로
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

df$pathway <- gsub(df$pathway, pattern = "HALLMARK_", replacement = "") %>% gsub(., pattern = "_", replacement = " ")
df$pathway <- first_word_cap(df$pathway)

df$pathway <- factor(df$pathway, levels = df$pathway)

fileidentity <- "Fig2F"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(10), height = cm_to_inch(6.8), family = "Arial")
ggplot(df, aes(diff, pathway)) + 
  geom_point(aes(size = -log10(Padj), color = Category)) +
  theme_classic(base_family = "Arial") +
  xlab("Module score difference") +
  theme(axis.text.x = element_text(hjust = 1,size = 5, angle = 45,color="black"),
        line = element_line(linewidth = 0.3),
        axis.text.y = element_text(size = 5, color="black"), 
        axis.title.x = element_text(size= 5, color="black"), 
        legend.key.height = unit(0.2, "cm"),
        axis.title.y = element_blank(),
        legend.key.width  = unit(0.2, "cm"),
        legend.spacing.x  = unit(0.02, "cm"),
        #legend.ticks = element_blank(),
        #axis.title = element_text(colour = "black", size = size_pt-1),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        #plot.title = element_text(hjust = 0.5, size=6, color="black"),
        legend.title = element_text(size = 4, color ="black"), 
        legend.text = element_text(size = 4, color = "black")) +
  scale_color_manual(values = rev(c("coral2", "deepskyblue2")) ) +
  guides(size = guide_legend(title = expression(-log[10] ~ "(padj)"))) + scale_size(range = c(1, 3)) 
dev.off()

# Fig S2G ----

epi <- readRDS('R_objects/20250730_escc_squamouscells_only_hj.rds')
epi$pathology <- factor(epi$pathology, levels = c('Normal', 'Dysplasia', 'Microinvasive', 'Macroinvasive'))
tumor <- subset(epi, Annotation_v2 == 'Tumor')

library(purrr)
set.seed(1234); tumor <- NormalizeData(tumor, normalization.method = "LogNormalize", scale.factor = 10000) # log norm
set.seed(1234); tumor@meta.data <- tumor@meta.data %>% dplyr::mutate(donor_platform = paste(donor, platform, sep = "_")) # batch

set.seed(1234); obj_list <- SplitObject(tumor, split.by = "donor_platform") # split
set.seed(1234); obj_list <- obj_list[sapply(obj_list, ncol) >= 10] # exclude small size samples
# find variable gene for each object
set.seed(1234); obj_list <- map(obj_list, function(obj) { 
  set.seed(1234); obj <- NormalizeData(obj)
  set.seed(1234); obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  return(obj)
})

hvg_all <- obj_list %>% map(VariableFeatures) %>% unlist() %>% unique() # pull hvg
VariableFeatures(tumor) <- hvg_all # set hvg
set.seed(1234); tumor <- ScaleData(tumor, features = hvg_all, split.by = "donor_platform") # scaling
set.seed(1234); tumor <- RunPCA(tumor, features = hvg_all) # pca

tumor <- RunPCA(tumor, npcs = 30, verbose = T)

tumor_modules <- read.csv('hotspot/outs/20250813_escc_sc.snrna_tumor_hotspot_min_gene25_results_modules_km.csv')
tumor_modules <- na.omit(tumor_modules) %>% filter(Module != '-1')

for (i in unique(tumor_modules$Module)) {
  geneset <- tumor_modules %>% filter(Module == i) %>% pull(Gene)
  tumor <- AddModuleScore(tumor, features = list(geneset), name=paste0('Module',i))
}

saveRDS(tumor, 'R_objects/20250813_escc_sc.snrna_tumor_hotspot_min_gene25_module_scores.rds')

modules_v2 <- paste0('Module',unique(tumor_modules$Module),'1')

scores <- tumor@meta.data %>% select(pathology,modules_v2) %>% filter(pathology != 'Normal')

library(dplyr)
library(pheatmap)
library(tidyverse)

module_scores_by_pathology <- scores %>%
  group_by(pathology) %>%
  summarise(across(starts_with("Module"), mean, na.rm = TRUE)) %>%
  column_to_rownames("pathology")

heatmap_matrix <- t(as.matrix(module_scores_by_pathology))
row.names(heatmap_matrix) <- sub('.$','',row.names(heatmap_matrix))
heatmap_matrix <- heatmap_matrix[c('Module2','Module5','Module4','Module8','Module7','Module1','Module6','Module3','Module9'),]

library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

cairo_pdf("hotspot/outs/20250813_escc_sc.snrna_tumor_hotspot_min_gene25_module_scores_heatmap_figure_km.pdf", width = cm_to_inch(6), height = cm_to_inch(8), family = "Arial")
pheatmap(heatmap_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         scale = "row",  # optional: z-score scaling per module
         fontsize_row = 7,
         fontsize_col = 7,
         fontsize = 7)
dev.off()

## T-test ----
library(dplyr)
library(ggplot2)
library(stringr)
library(Seurat)


seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/km/20250813_escc_sc.snrna_tumor_hotspot_min_gene25_module_scores.rds")

metadat <- seurat_obj@meta.data
seurat_obj <- subset(seurat_obj, subset = pathology != "Normal")

df <- data.frame()

for (mymodule in paste0("Module", 1:9, 1)) {
  
  for (stage_tmp in c("Dysplasia", "Microinvasive", "Macroinvasive")) {
    
    metadat_tmp <- metadat %>% dplyr::filter(pathology == stage_tmp)
    group1 <- metadat_tmp[, mymodule]
    
    metadat_tmp <- metadat %>% dplyr::filter(pathology != stage_tmp)
    group2 <- metadat_tmp[, mymodule]
    
    res <- t.test(x = group1, y = group2)
    pval_tmp <- res$p.value
    
    df_tmp <- data.frame(module = mymodule, group1 = stage_tmp, group2 = "others", pval = pval_tmp, group1_mean = mean(group1), group2_mean = mean(group2))
    df <- rbind(df, df_tmp)
    
  }
}

df$padj <- p.adjust(df$pval, method = "BH")

df <- df %>% dplyr::mutate(diff = group1_mean - group2_mean)
#df <- df %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(diff > 0)

library(WriteXLS)
WriteXLS(df, "/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig2/20250826_escc_hotspot_module_ttest_hj.xlsx")

library(readxl)
df <- read_xlsx("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig2/20250826_escc_hotspot_module_ttest_hj.xlsx")


mycomb <- combn(c( "Dysplasia", "Microinvasive", "Macroinvasive"), 2)
mycomb[,1]

df <- data.frame()

for (mymodule in paste0("Module", 1:9, 1)) {
  
  for (i in 1:3) {
    
    mycomb_tmp <- mycomb[,i]
    stage_tmp1 <- mycomb_tmp[1]
    stage_tmp2 <- mycomb_tmp[2]
    
    metadat_tmp <- metadat %>% dplyr::filter(pathology == stage_tmp1)
    group1 <- metadat_tmp[, mymodule]
    
    metadat_tmp <- metadat %>% dplyr::filter(pathology == stage_tmp2)
    group2 <- metadat_tmp[, mymodule]
    
    res <- t.test(x = group1, y = group2)
    pval_tmp <- res$p.value
    
    df_tmp <- data.frame(module = mymodule, group1 = stage_tmp1, group2 = stage_tmp2, pval = pval_tmp, group1_mean = mean(group1), group2_mean = mean(group2))
    df <- rbind(df, df_tmp)
    
  }
}

df$padj <- p.adjust(df$pval, method = "BH")

df <- df %>% dplyr::mutate(diff = group1_mean - group2_mean)
df <- df %>% dplyr::filter(padj < 0.05)

# Fig 2SH ----

visium <- readRDS("R_objects/20250805_escc_merged_rctd_c2l_dl.rds")
tumor_modules <- read.csv('hotspot/outs/20250803_escc_sc.snrna_tumor_hotspot_min_gene25_results_modules_km.csv')
tumor_modules <- na.omit(tumor_modules) %>% filter(Module != '-1')


for (i in unique(tumor_modules$Module)) {
  geneset <- tumor_modules %>% filter(Module == i) %>% pull(Gene)
  visium <- AddModuleScore(visium, features = list(geneset), name=paste0('Module',i), assay = 'Spatial')
}

modules_v2 <- paste0('Module',unique(tumor_modules$Module),'1')

library(ggplot2)
library(viridis)
SpatialFeaturePlot(visium, modules_v2[1], ncol = 3, pt.size = 0.1) & 
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))

# Fig 2SI ----

library(MAST)
library(Seurat)
library(dplyr)
library(showtext)
library(ggplot2)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

date = "20250816"
project_name = "escc"

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig2/SERPINE1/")
seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")

seurat_obj <- subset(seurat_obj, subset = Annotation_v2 == "Tumor")

Idents(seurat_obj) <- "pathology"
res <- FindAllMarkers(seurat_obj, only.pos = T, test.use = "MAST", latent.vars = "donor")

saveRDS(res, "dys_marker_MAST_hj.rds")

myres <- res %>% dplyr::group_by(cluster) %>% dplyr::arrange(p_val_adj) %>% dplyr::filter(avg_log2FC > 1.5) %>%  dplyr::group_by(cluster) %>% dplyr::slice(1:10)
myres <- res %>% dplyr::group_by(cluster) %>% dplyr::arrange(p_val_adj) %>% dplyr::filter(avg_log2FC > 1) %>%  dplyr::group_by(cluster) %>% dplyr::slice(1:10)

DefaultAssay(seurat_obj)

fileidentity <- "tumor_stage_DEtest_dotplot"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(12), height = cm_to_inch(4), family = "Arial")
DotPlot(seurat_obj, features = c(myres$gene), scale = T, dot.scale = 2) & scale_color_gradientn(na.value = "white", colours=c("blue","white", "red"), name = 'Average\nExpression') & theme_classic(base_family = "Arial", base_size = 6) & theme(axis.title = element_text(size = 6), line = element_line(linewidth = 0.3),aspect.ratio = 0.3, panel.grid = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.key.size = unit(0.2, units = "cm"), legend.title = element_text(size = 4),legend.key.height = unit(0.15, units = "cm"))
dev.off()

# Fig 2SJ ----
library(MAST)
library(Seurat)
library(dplyr)
library(showtext)
library(ggplot2)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

date = "20250816"
project_name = "escc"

setwd("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/Fig2/SERPINE1/")
seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")

seurat_obj <- subset(seurat_obj, subset = Annotation_v2 == "Tumor")

Idents(seurat_obj) <- "pathology"
res <- FindAllMarkers(seurat_obj, only.pos = T, test.use = "MAST", latent.vars = "donor")

saveRDS(res, "dys_marker_MAST_hj.rds")

myres <- res %>% dplyr::group_by(cluster) %>% dplyr::arrange(p_val_adj) %>% dplyr::filter(avg_log2FC > 1.5) %>%  dplyr::group_by(cluster) %>% dplyr::slice(1:10)
myres <- res %>% dplyr::group_by(cluster) %>% dplyr::arrange(p_val_adj) %>% dplyr::filter(avg_log2FC > 1) %>%  dplyr::group_by(cluster) %>% dplyr::slice(1:10)

DefaultAssay(seurat_obj)

fileidentity <- "tumor_stage_DEtest_dotplot"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = cm_to_inch(12), height = cm_to_inch(4), family = "Arial")
DotPlot(seurat_obj, features = c(myres$gene), scale = T, dot.scale = 2) & scale_color_gradientn(na.value = "white", colours=c("blue","white", "red"), name = 'Average\nExpression') & theme_classic(base_family = "Arial", base_size = 6) & theme(axis.title = element_text(size = 6), line = element_line(linewidth = 0.3),aspect.ratio = 0.3, panel.grid = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.key.size = unit(0.2, units = "cm"), legend.title = element_text(size = 4),legend.key.height = unit(0.15, units = "cm"))
dev.off()

# Fig 2SJ ----
seurat_obj <- readRDS("/Users/hojin/Dropbox/project/ESCC/submit/figures/escc_figure_share/object/sk/20250812_squamous_final_sk.rds")

FeaturePlot(seurat_obj, features = "SERPINE1", split.by = "pathology", cols = c("grey", "red"),order = F) & theme_void() & theme(aspect.ratio = 1, plot.title = element_blank()) & NoLegend()

FeaturePlot(seurat_obj, features = "ANKRD37", split.by = "pathology", cols = c("grey", "red"),order = F) & theme_void() & theme(aspect.ratio = 1, plot.title = element_blank()) & NoLegend()

VlnPlot(seurat_obj, features = "SERPINE1", group.by = "pathology", split.by = "Annotation_v2", cols = c("Basal" = "#F8766D", "Suprabasal" = "#0CB702", "Tumor" = "#619CFF"), pt.size = 0) & theme_classic(base_family = "Arial") & theme(line = element_line(linewidth = 0.2), axis.ticks.length = unit(0.5, "mm"), aspect.ratio = 0.7, plot.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), axis.text = element_text(size = 5), axis.title = element_text(size = 5), legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.2, "cm"), legend.text = element_text(size = 3)) 

VlnPlot(seurat_obj, features = "ANKRD37", group.by = "pathology", split.by = "Annotation_v2", cols = c("Basal" = "#F8766D", "Suprabasal" = "#0CB702", "Tumor" = "#619CFF"), pt.size = 0) & theme_classic(base_family = "Arial") & theme(line = element_line(linewidth = 0.2), axis.ticks.length = unit(0.5, "mm"), aspect.ratio = 0.7, plot.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), axis.text = element_text(size = 5), axis.title = element_text(size = 5), legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.2, "cm"), legend.text = element_text(size = 3)) 


