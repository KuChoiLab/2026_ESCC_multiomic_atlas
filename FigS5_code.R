# Fig. S5A -----
## 2. scaling ----
res_sc_df <- readRDS("scenic/cd4t/downstream/20250821_escc_scenic_cd4t_sc_test_km.rds")
res_sc_df$platform <- "sc"
res_sn_df <- readRDS("scenic/cd4t/downstream/20250821_escc_scenic_cd4t_sn_test_km.rds")
res_sn_df$platform <- "sn"

res_df <- rbind(res_sc_df, res_sn_df)
quant <- quantile(res_df$padj, 0.1)

filtered_df <- res_df %>%
  dplyr::group_by(group1, platform) %>%
  dplyr::filter(padj < quantile(padj, 0.15, na.rm = TRUE)) %>%
  ungroup()

res_df <- filtered_df %>% dplyr::mutate(log2FC = log2(group1_mean/group2_mean)) %>% 
  dplyr::filter(group2_mean < group1_mean, padj < 0.05) #%>% dplyr::filter(log2FC > 0.7)

res_counts <- res_df %>% dplyr::group_by(TF, group1) %>% dplyr::summarise(n=n()) %>% dplyr::arrange(desc(n))
mytfs <- res_counts %>% dplyr::filter(n == 2) %>% dplyr::pull(TF) %>% unique()


## 3. merge sc, sn ----
seurat_obj <- readRDS("R_objects/20250729_escc_cd4_celltype_annotated_dl.rds")

seurat_obj <- subset(seurat_obj, subset = Annotation_v2 %in% c("CD4 CTL","CD4 TNFRSF9- Treg","CD4 TNFRSF9+ Treg"))

seurat_obj <- subset(seurat_obj, subset = platform == "scRNA")
regulonAUC_scaled_sc <- getAUC(regulonAUC)[mytfs, colnames(seurat_obj)]

regulonAUC_scaled_sc <- t(scale(t(regulonAUC_scaled_sc), center = T, scale=T))
# regulonAUC_scaled_sc <- t(apply(regulonAUC_scaled_sc, 1, function(x) {
#   (x - min(x)) / (max(x) - min(x))
# }))
#regulonAUC_scaled_sc <- t(scale(t(getAUC(regulonAUC)[res_df_filtered$TF, colnames(seurat_obj)]), center = T, scale=T))

seurat_obj <- readRDS("R_objects/20250729_escc_cd4_celltype_annotated_dl.rds")


seurat_obj <- subset(seurat_obj, subset = Annotation_v2 %in% c("CD4 CTL","CD4 TNFRSF9- Treg","CD4 TNFRSF9+ Treg"))

seurat_obj <- subset(seurat_obj, subset = platform == "snRNA")
regulonAUC_scaled_sn <- getAUC(regulonAUC)[mytfs, colnames(seurat_obj)]
regulonAUC_scaled_sn <- t(scale(t(regulonAUC_scaled_sn), center = T, scale=T))
regulonAUC_scaled <- cbind(regulonAUC_scaled_sn, regulonAUC_scaled_sc)

seurat_obj <- readRDS("R_objects/20250729_escc_cd4_celltype_annotated_dl.rds")

seurat_obj <- subset(seurat_obj, subset = Annotation_v2 %in% c("CD4 CTL","CD4 TNFRSF9- Treg","CD4 TNFRSF9+ Treg"))

celltypes <- seurat_obj@meta.data %>% select(Annotation_v2)
celltypes$New_Identity <- factor(celltypes$Annotation_v2, levels = c("CD4 CTL","CD4 TNFRSF9- Treg","CD4 TNFRSF9+ Treg"))
selectedResolution <- "New_Identity"
cellsPerTypes <- split(rownames(celltypes), celltypes[,selectedResolution]) 

## 4. scenic heatmap -----

cellorder <- unlist(cellsPerTypes)

cell_type <- data.frame(
  "cell type" = unlist(mapply(rep, names(cellsPerTypes), lapply(cellsPerTypes, length))))

regulonAUC_scaled_order <- regulonAUC_scaled[mytfs, cellorder]
row.names(cell_type) <- colnames(regulonAUC_scaled_order)
palette_length = 100

date <- "20250823"
project_name <- "escc"
cm_to_inch <- function(cm) cm / 2.54
fileidentity <- "cd4t_scenic"

library(CellChat)
mycol <- scPalette(3)
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
                       annotation_name_gp = gpar(fontsize = 8),
                       annotation_legend_param = list(
                         title_gp = gpar(fontsize = 8, fontfamily = "Arial"),
                         labels_gp = gpar(fontsize = 8, fontfamily = "Arial")
                       ))

col_fun = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "red"))

png(filename = paste0("scenic/cd4t/downstream/",date, "_", project_name, "_", fileidentity, "_km.png"),
    width = 8, height = 2, units = "in", res = 600, family = "Arial")
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
                             column_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
                             heatmap_legend_param = list(
                               title = "Score",
                               title_gp = gpar(fontsize = 8, fontfamily = "Arial"),
                               labels_gp = gpar(fontsize = 8, fontfamily = "Arial")
                             ),
                             top_annotation = ha)
draw(ht, use_raster = TRUE)
dev.off()

cairo_pdf(paste0("scenic/cd4t/downstream/",date, "_", project_name, "_", fileidentity, "_km.pdf"),
          width = 8, height = 2, family = "Arial")

col_fun <- colorRamp2(c(-1, 0, 1), c("darkblue", "white", "red"))

ht <- Heatmap(
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
  row_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
  heatmap_legend_param = list(
    title = "Score",
    title_gp = gpar(fontsize = 8, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 8, fontfamily = "Arial")
  ),
  top_annotation = ha
)

draw(ht)
dev.off()


## 5. ORA ----
library(org.Hs.eg.db)
library(msigdbr)
library(clusterProfiler)

regulons <- readRDS('scenic/cd4t/int/2.6_regulons_asGeneSet.Rds')
#GOBP
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, entrez_gene)
#Hallmark
m_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

hs <- org.Hs.eg.db
row_order_vec <- row_order(ht)
ordered_row_names <- rownames(regulonAUC_scaled_order)[row_order_vec]

tfs <- gsub(ordered_row_names, pattern = " \\([0-9]*g\\)", replacement = "")

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

library(openxlsx)
write.xlsx(enricher_res, "scenic/cd4t/downstream/20250823_escc_cd4t_scenic_targetgenes_gobp_km.xlsx")
enricher_res <- read.xlsx("scenic/cd4t/downstream/20250823_escc_cd4t_scenic_targetgenes_gobp_km.xlsx")
enricher_res <- enricher_res %>% dplyr::filter(p.adjust < 0.05)

write.xlsx(enricher_res, "scenic/cd4t/downstream/20250823_escc_cd4t_scenic_targetgenes_gobp_sig_km.xlsx")
enricher_res <- enricher_res %>% dplyr::group_by(TFs) %>% dplyr::arrange(p.adjust, .by_group = T) %>% dplyr::slice(1)
rownames(enricher_res) <- enricher_res$TFs

tf_filtered <- intersect(tfs, enricher_res$TFs)
enricher_res <- enricher_res[tf_filtered,] # sorting top1 results

write.xlsx(enricher_res, "scenic/cd4t/downstream/20250823_escc_cd4t_scenic_targetgenes_gobp_top1_km.xlsx")
enricher_res <- read.xlsx("scenic/cd4t/downstream/20250823_escc_cd4t_scenic_targetgenes_gobp_top1_km.xlsx")


## 6. GOBP plot ----
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

rownames(mat) <- sub("^GOBP_", "", rownames(mat))        
rownames(mat) <- gsub("_", " ", rownames(mat))             
rownames(mat) <- str_to_lower(rownames(mat))      
rownames(mat) <- str_to_sentence(rownames(mat))
row_lbl <- rownames(mat)

rng <- range(val, na.rm = TRUE)
if (use_log) {
  col_fun <- colorRamp2(c(0, rng[2]), c("white", "red"))
} else {
  col_fun <- colorRamp2(c(rng[1], rng[2]), c("white", "red"))
}


ht <- Heatmap(
  mat,
  name = val_name,
  col  = col_fun,
  na_col = "darkgrey",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_labels = row_lbl,
  row_names_gp = gpar(fontsize = 9, fontfamily = "Arial"),
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(col = "black", fill = NA, lwd = 0.5))
  },
  width = unit(ncol(mat)/4, "cm"), 
  height = unit(nrow(mat)/4, "cm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 9, fontfamily = "Arial"),
    legend_height = unit(1, "cm"),                        
    grid_height = unit(0.2, "cm"),                        
    grid_width = unit(0.1, "cm")                          
  )
)

# save plots
fileidentity <- "scenic_top1_gobp"
cairo_pdf(paste0("scenic/cd4t/downstream/",date, "_", project_name, "_", fileidentity, "_km.pdf"),
          width = 8, height = 2, family = "Arial")
draw(
  ht,
  heatmap_legend_side = "left"
)
dev.off()


# Fig. S5B ----
##correlation btw CD4 CTL and CDH17+ cDC2####
#myeloid prop.
myeloid <- readRDS('R_objects/20250722_escc_myeloid_flt_km.rds')
# Extract current cell type info
celltypes <- myeloid$annotation_myeloid

# Logical index for the specific cell type
target_idx <- celltypes == "CD1C+ cDC2"

# Get gene expression for those cells
gene_expr <- FetchData(myeloid, vars = "CDH17")[,1]

# Define new labels for the target cell type
new_labels <- celltypes
new_labels[target_idx & gene_expr > 0] <- "CDH17+ cDC2"
new_labels[target_idx & gene_expr == 0] <- "CDH17- cDC2"

# Add back to metadata
myeloid$annotation_myeloid_v2 <- new_labels

Idents(myeloid) <- 'orig.ident'
table_myeloid <- myeloid@meta.data %>% select(orig.ident, annotation_myeloid_v2) %>% table()
table_myeloid_v2 <- ( table_myeloid / table_myeloid %>% apply(1,sum) ) %>% data.frame()
table_myeloid_v2 <- table_myeloid_v2 %>% filter(annotation_myeloid_v2 == 'CDH17+ cDC2')
table_myeloid_v2[is.na(table_myeloid_v2)] <- 0
colnames(table_myeloid_v2)[3] <- 'Freq_myeloid'

#CD4 CTL prop.
cd4t <-readRDS('R_objects/20250729_escc_cd4_celltype_annotated_dl.rds')
Idents(cd4t) <- 'orig.ident'
table_cd4t <- cd4t@meta.data %>% select(orig.ident, Annotation_v2) %>% table()
table_cd4t_v2 <- ( table_cd4t / table_cd4t %>% apply(1,sum) ) %>% data.frame()
table_cd4t_v2 <- table_cd4t_v2 %>% filter(Annotation_v2 == 'CD4 CTL')
table_cd4t_v2[is.na(table_cd4t_v2)] <- 0
colnames(table_cd4t_v2)[3] <- 'Freq_CD4T'

table_myeloid_cd4t_v2 <- merge(table_myeloid_v2, table_cd4t_v2, by = 'orig.ident', all = TRUE)
table_myeloid_cd4t_v2[is.na(table_myeloid_cd4t_v2)] <- 0
table_myeloid_cd4t_v2

cDC2_CTL <- data.frame(
  cDC2 = table_myeloid_cd4t_v2$Freq_myeloid,
  CTL = table_myeloid_cd4t_v2$Freq_CD4T,
  sample = table_myeloid_cd4t_v2$orig.ident)

sample_pathology <- cd4t@meta.data %>% select(orig.ident, pathology) %>% unique()
colnames(sample_pathology)[1] <- 'sample'

cDC2_CTL_v2 <- merge(cDC2_CTL, sample_pathology, by = 'sample')

cDC2_CTL_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',cDC2_CTL_v2$pathology)
cDC2_CTL_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',cDC2_CTL_v2$pathology)
cDC2_CTL_v2$pathology <- factor(cDC2_CTL_v2$pathology, levels = c('Normal','Dysplasia','Microinvasive carcinoma','Macroinvasive carcinoma'))

write.xlsx(cDC2_CTL_v2, '~/Dropbox/escc_figure_share/Fig5/20250816_escc_sc.snrna_CD4CTL_cDC2_corr_km.xlsx')
cDC2_CTL_v2 <- read.xlsx('~/Dropbox/escc_figure_share/Fig5/20250816_escc_sc.snrna_CD4CTL_cDC2_corr_km.xlsx')

library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

#V2
cairo_pdf("~/Dropbox/escc_figure_share/FigS5/20250822/20250823_escc_sc.snrna_CD4CTL_cDC2_corr_km.pdf", 
          width = 6, height = 5, family = "Arial")
ggplot(cDC2_CTL_v2, aes(x = cDC2, y = CTL)) +
  geom_point(aes(shape = pathology, color = pathology), size = 5) +
  scale_color_manual(values = c('#83f52c', '#f3f315', '#ff6600', '#ff0099')) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  geom_smooth(method = 'lm', color ='red', se = FALSE) +
  annotate("text", x = 0.1, y = 0.8, col = "black", size = 7,
           label = paste("r =", signif(cor.test(cDC2_CTL_v2$cDC2, cDC2_CTL_v2$CTL)$estimate,3), ", p =", signif(cor.test(cDC2_CTL_v2$cDC2, cDC2_CTL_v2$CTL)$p.value,3)),
           family = 'Arial')+
  theme_classic(base_family = 'Arial') +
  labs(x = 'CDH17+ cDC2 fraction', y = 'CD4 CTL fraction') +
  theme(axis.title.x = element_text(size = 20, family = 'Arial', colour = 'black'),
        axis.text.x = element_text(size = 20, family = 'Arial', colour = 'black'),
        axis.title.y = element_text(size = 20, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 20, family = 'Arial', colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(size = 20, family="Arial"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.title = element_blank()) &
  guides(shape=guide_legend(override.aes = list(size = 5)),
         col=guide_legend(ncol = 2))
dev.off()



##correlation btw CD4 CTL and IL15+ cDC3####
#myeloid prop.
myeloid <- readRDS('R_objects/20250722_escc_myeloid_flt_km.rds')
# Extract current cell type info
celltypes <- myeloid$annotation_myeloid

# Logical index for the specific cell type
target_idx <- celltypes == "LAMP3+ cDC3"

# Get gene expression for those cells
gene_expr <- FetchData(myeloid, vars = "IL15")[,1]

# Define new labels for the target cell type
new_labels <- celltypes
new_labels[target_idx & gene_expr > 0] <- "IL15+ cDC3"
new_labels[target_idx & gene_expr == 0] <- "IL15- cDC3"

# Add back to metadata
myeloid$annotation_myeloid_v2 <- new_labels

Idents(myeloid) <- 'orig.ident'
table_myeloid <- myeloid@meta.data %>% select(orig.ident, annotation_myeloid_v2) %>% table()
table_myeloid_v2 <- ( table_myeloid / table_myeloid %>% apply(1,sum) ) %>% data.frame()
table_myeloid_v2 <- table_myeloid_v2 %>% filter(annotation_myeloid_v2 == 'IL15+ cDC3')
table_myeloid_v2[is.na(table_myeloid_v2)] <- 0
colnames(table_myeloid_v2)[3] <- 'Freq_myeloid'

#CD4 CTL prop.
cd4t <-readRDS('R_objects/20250729_escc_cd4_celltype_annotated_dl.rds')
Idents(cd4t) <- 'orig.ident'
table_cd4t <- cd4t@meta.data %>% select(orig.ident, Annotation_v2) %>% table()
table_cd4t_v2 <- ( table_cd4t / table_cd4t %>% apply(1,sum) ) %>% data.frame()
table_cd4t_v2 <- table_cd4t_v2 %>% filter(Annotation_v2 == 'CD4 CTL')
table_cd4t_v2[is.na(table_cd4t_v2)] <- 0
colnames(table_cd4t_v2)[3] <- 'Freq_CD4T'

table_myeloid_cd4t_v2 <- merge(table_myeloid_v2, table_cd4t_v2, by = 'orig.ident', all = TRUE)
table_myeloid_cd4t_v2[is.na(table_myeloid_cd4t_v2)] <- 0
table_myeloid_cd4t_v2

cDC3_CTL <- data.frame(
  cDC3 = table_myeloid_cd4t_v2$Freq_myeloid,
  CTL = table_myeloid_cd4t_v2$Freq_CD4T,
  sample = table_myeloid_cd4t_v2$orig.ident)

sample_pathology <- cd4t@meta.data %>% select(orig.ident, pathology) %>% unique()
colnames(sample_pathology)[1] <- 'sample'

cDC3_CTL_v2 <- merge(cDC3_CTL, sample_pathology, by = 'sample')

cDC3_CTL_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',cDC3_CTL_v2$pathology)
cDC3_CTL_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',cDC3_CTL_v2$pathology)
cDC3_CTL_v2$pathology <- factor(cDC3_CTL_v2$pathology, levels = c('Normal','Dysplasia','Microinvasive carcinoma','Macroinvasive carcinoma'))

write.xlsx(cDC3_CTL_v2, '~/Dropbox/escc_figure_share/Fig5/20250816_escc_sc.snrna_CD4CTL_cDC3_corr_km.xlsx')
cDC3_CTL_v2 <- read.xlsx('~/Dropbox/escc_figure_share/Fig5/20250816_escc_sc.snrna_CD4CTL_cDC3_corr_km.xlsx')
library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

#V2
cairo_pdf("~/Dropbox/escc_figure_share/FigS5/20250822/20250823_escc_sc.snrna_CD4CTL_cDC3_corr_km.pdf", 
          width = 6, height = 5, family = "Arial")
ggplot(cDC3_CTL_v2, aes(x = cDC3, y = CTL)) +
  geom_point(aes(shape = pathology, color = pathology), size = 5) +
  scale_color_manual(values = c('#83f52c', '#f3f315', '#ff6600', '#ff0099')) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  geom_smooth(method = 'lm', color ='red', se = FALSE) +
  annotate("text", x = 0.015, y = 0.8, col = "black", size = 7,
           label = paste("r =", signif(cor.test(cDC3_CTL_v2$cDC3, cDC3_CTL_v2$CTL)$estimate,3), ", p =", signif(cor.test(cDC3_CTL_v2$cDC3, cDC3_CTL_v2$CTL)$p.value,3)),
           family = 'Arial')+
  theme_classic(base_family = 'Arial') +
  labs(x = 'IL15+ cDC3 fraction', y = 'CD4 CTL fraction') +
  theme(axis.title.x = element_text(size = 20, family = 'Arial', colour = 'black'),
        axis.text.x = element_text(size = 20, family = 'Arial', colour = 'black'),
        axis.title.y = element_text(size = 20, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 20, family = 'Arial', colour = 'black'),
        legend.position = 'top',
        legend.text = element_text(size = 20, family="Arial"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.title = element_blank()) &
  guides(shape=guide_legend(override.aes = list(size = 5)),
         col = guide_legend(ncol = 2))
dev.off()



##correlation btw CD4 CTL and MHC II + tumor cells####
###Module score####
#epi prop.
epi <- readRDS('R_objects/20250730_escc_squamouscells_only_hj.rds')
epi <- AddModuleScore(epi, features = list(c('HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1')), name='MHC_II')


# Extract current cell type info
celltypes <- epi$Annotation_v2

# Logical index for the specific cell type
tumor_idx <- celltypes == "Tumor"

# Apply threshold to split tumor cells
celltypes[tumor_idx & epi$MHC_II1 > 0]  <- "MHC II+ tumor"
celltypes[tumor_idx & epi$MHC_II1 <= 0] <- "MHC II- tumor"

# Define new labels for the target cell type
epi$Annotation_v3 <- celltypes

Idents(epi) <- 'orig.ident'
table_epi <- epi@meta.data %>% select(orig.ident, Annotation_v3) %>% table()
table_epi_v2 <- ( table_epi / table_epi %>% apply(1,sum) ) %>% data.frame()
table_epi_v2 <- table_epi_v2 %>% filter(Annotation_v3 == 'MHC II+ tumor')
table_epi_v2[is.na(table_epi_v2)] <- 0
colnames(table_epi_v2)[3] <- 'Freq_tumor'

#CD4 CTL prop.
cd4t <-readRDS('R_objects/20250729_escc_cd4_celltype_annotated_dl.rds')
Idents(cd4t) <- 'orig.ident'
table_cd4t <- cd4t@meta.data %>% select(orig.ident, Annotation_v2) %>% table()
table_cd4t_v2 <- ( table_cd4t / table_cd4t %>% apply(1,sum) ) %>% data.frame()
table_cd4t_v2 <- table_cd4t_v2 %>% filter(Annotation_v2 == 'CD4 CTL')
table_cd4t_v2[is.na(table_cd4t_v2)] <- 0
colnames(table_cd4t_v2)[3] <- 'Freq_CD4T'

table_epi_cd4t_v2 <- merge(table_epi_v2, table_cd4t_v2, by = 'orig.ident', all = TRUE)
table_epi_cd4t_v2[is.na(table_epi_cd4t_v2)] <- 0
table_epi_cd4t_v2

tumor_CTL <- data.frame(
  Tumor = table_epi_cd4t_v2$Freq_tumor,
  CTL = table_epi_cd4t_v2$Freq_CD4T,
  sample = table_epi_cd4t_v2$orig.ident)

sample_pathology <- cd4t@meta.data %>% select(orig.ident, pathology) %>% unique()
colnames(sample_pathology)[1] <- 'sample'

tumor_CTL_v2 <- merge(tumor_CTL, sample_pathology, by = 'sample')

tumor_CTL_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',tumor_CTL_v2$pathology)
tumor_CTL_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',tumor_CTL_v2$pathology)
tumor_CTL_v2$pathology <- factor(tumor_CTL_v2$pathology, levels = c('Normal','Dysplasia','Microinvasive carcinoma','Macroinvasive carcinoma'))

write.xlsx(tumor_CTL_v2, '~/Dropbox/escc_figure_share/Fig5/20250817_escc_sc.snrna_CD4CTL_tumor_corr_km.xlsx')
tumor_CTL_v2 <- read.xlsx('~/Dropbox/escc_figure_share/Fig5/20250817_escc_sc.snrna_CD4CTL_tumor_corr_km.xlsx')
library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cairo_pdf("~/Dropbox/escc_figure_share/Fig5/20250816_escc_sc.snrna_CD4CTL_tumor_corr_km.pdf", 
          width = 9, height = 4, family = "Arial")
ggplot(tumor_CTL_v2, aes(x = Tumor, y = CTL)) +
  geom_point(aes(shape = pathology, color = pathology)) +
  scale_color_manual(values = c('#83f52c', '#f3f315', '#ff6600', '#ff0099')) +
  geom_smooth(method = 'lm', color ='red') +
  annotate("text", x = 0.3, y = 0.8, col = "black", size = 7,
           label = paste("r =", signif(cor.test(tumor_CTL_v2$Tumor, tumor_CTL_v2$CTL)$estimate,3), ", p =", signif(cor.test(tumor_CTL_v2$Tumor, tumor_CTL_v2$CTL)$p.value,3)))+
  theme_classic(base_family = 'Arial') +
  labs(x = 'MHC II+ tumor fraction', y = 'CD4 CTL fraction') +
  theme(axis.title.x = element_text(size = 17, family = 'Arial', colour = 'black'),
        axis.text.x = element_text(size = 17, family = 'Arial', colour = 'black'),
        axis.title.y = element_text(size = 17, family = 'Arial', colour = 'black'),
        axis.text.y = element_text(size = 17, family = 'Arial', colour = 'black'),
        legend.position = 'right',
        legend.text = element_text(size = 17, family="Arial"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.title = element_blank()) &
  guides(shape=guide_legend(override.aes = list(size = 5)))
dev.off()




#Fig. S5C####
cd4t <- readRDS('R_objects/20250729_escc_cd4_celltype_annotated_dl.rds')
seurat_obj_subset <- cd4t

## DE test ####
# CD4 CTL vs Others
seurat_obj_subset@meta.data <- seurat_obj_subset@meta.data %>% dplyr::mutate(nebula_anno = case_when(Annotation_v2 == "CD4 CTL" ~ "CD4 CTL",
                                                                                                     T ~ "Others"))

seurat_obj_subset$nebula_anno <- factor(seurat_obj_subset$nebula_anno, levels = rev(c("CD4 CTL", "Others"))) # 중요

seuratdata <- scToNeb(obj = seurat_obj_subset, assay = "RNA", id = "donor", pred = c("nebula_anno", "donor"), offset = "nCount_RNA")
df = model.matrix(~nebula_anno+donor, data = seuratdata$pred)
## include only the first two cell types in the model to avoid separation due to too many binary variables
#data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoCD4 CTL")],offset=seuratdata$offset)
#re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset)
re = nebula(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoCD4 CTL")],offset=seuratdata$offset)

saveRDS(re, "nebula/cd4t/nebula_cd4ctl_km.rds")
re <- readRDS("nebula/cd4t/nebula_cd4ctl_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoCD4 CTL`, `se_nebula_annoCD4 CTL`, `p_nebula_annoCD4 CTL`)
res$Padj <- p.adjust(res$`p_nebula_annoCD4 CTL`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
write.xlsx(res,"nebula/cd4t/20251110_escc_scrna_nebula_cd4ctl_km.xlsx")

# CD4 TNFRSF9+ Treg vs Others
seurat_obj_subset@meta.data <- seurat_obj_subset@meta.data %>% dplyr::mutate(nebula_anno = case_when(Annotation_v2 == "CD4 TNFRSF9+ Treg" ~ "CD4 TNFRSF9+ Treg",
                                                                                                     T ~ "Others"))

seurat_obj_subset$nebula_anno <- factor(seurat_obj_subset$nebula_anno, levels = rev(c("CD4 TNFRSF9+ Treg", "Others"))) # 중요

seuratdata <- scToNeb(obj = seurat_obj_subset, assay = "RNA", id = "donor", pred = c("nebula_anno", "donor"), offset = "nCount_RNA")
df = model.matrix(~nebula_anno+donor, data = seuratdata$pred)
colnames(df)
## include only the first two cell types in the model to avoid separation due to too many binary variables
#data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoCD4 TNFRSF9+ Treg")],offset=seuratdata$offset)
#re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset)
re = nebula(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoCD4 TNFRSF9+ Treg")],offset=seuratdata$offset)

saveRDS(re, "nebula/cd4t/nebula_cd4tnfrsf9postreg_km.rds")
re <- readRDS("nebula/cd4t/nebula_cd4tnfrsf9postreg_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoCD4 TNFRSF9+ Treg`, `se_nebula_annoCD4 TNFRSF9+ Treg`, `p_nebula_annoCD4 TNFRSF9+ Treg`)
res$Padj <- p.adjust(res$`p_nebula_annoCD4 TNFRSF9+ Treg`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
write.xlsx(res,"nebula/cd4t/20251110_escc_scrna_nebula_cd4tnfrsf9postreg_km.xlsx")


# CD4 TNFRSF9- Treg vs Others
seurat_obj_subset@meta.data <- seurat_obj_subset@meta.data %>% dplyr::mutate(nebula_anno = case_when(Annotation_v2 == "CD4 TNFRSF9- Treg" ~ "CD4 TNFRSF9- Treg",
                                                                                                     T ~ "Others"))

seurat_obj_subset$nebula_anno <- factor(seurat_obj_subset$nebula_anno, levels = rev(c("CD4 TNFRSF9- Treg", "Others")))

seuratdata <- scToNeb(obj = seurat_obj_subset, assay = "RNA", id = "donor", pred = c("nebula_anno", "donor"), offset = "nCount_RNA")
df = model.matrix(~nebula_anno+donor, data = seuratdata$pred)
colnames(df)
## include only the first two cell types in the model to avoid separation due to too many binary variables
#data_g = group_cell(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoCD4 TNFRSF9- Treg")],offset=seuratdata$offset)
#re = nebula(data_g$count, data_g$id, pred=data_g$pred, offset=data_g$offset)
re = nebula(count=seuratdata$count,id=seuratdata$id,pred=df[,c("(Intercept)","nebula_annoCD4 TNFRSF9- Treg")],offset=seuratdata$offset)

saveRDS(re, "nebula/cd4t/nebula_cd4tnfrsf9negtreg_km.rds")
re <- readRDS("nebula/cd4t/nebula_cd4tnfrsf9negtreg_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoCD4 TNFRSF9- Treg`, `se_nebula_annoCD4 TNFRSF9- Treg`, `p_nebula_annoCD4 TNFRSF9- Treg`)
res$Padj <- p.adjust(res$`p_nebula_annoCD4 TNFRSF9- Treg`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
write.xlsx(res,"nebula/cd4t/20251110_escc_scrna_nebula_cd4tnfrsf9negtreg_km.xlsx")


##Volcano####
#CD4 TNFRSF9- Treg
re <- readRDS("nebula/cd4t/nebula_cd4tnfrsf9negtreg_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoCD4 TNFRSF9- Treg`, `se_nebula_annoCD4 TNFRSF9- Treg`, `p_nebula_annoCD4 TNFRSF9- Treg`)
res$Padj <- p.adjust(res$`p_nebula_annoCD4 TNFRSF9- Treg`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
pos <- res %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% head(5) %>% pull(gene)
neg <- res %>% filter(p_val_adj < 0.05) %>% arrange(avg_log2FC) %>% head(5) %>% pull(gene)

cd4tnfrsf9negtreg_DEG_volcano <- EnhancedVolcano(res, lab = res$gene, x = "avg_log2FC", y = "p_val_adj",
                                                 title = 'Other CD4 T cells vs. CD4 TNFRSF9- Treg',
                                                 selectLab = c(neg,pos),
                                                 subtitle = NULL,
                                                 labSize = 6.0,
                                                 labCol = 'black',
                                                 pCutoff = 0.05,
                                                 FCcutoff = 1,
                                                 pointSize = 3.0,
                                                 drawConnectors = TRUE,
                                                 legendLabels = c("NS", "avg_Log2FC", "P.adj", "P.adj and avg_Log2FC"),
                                                 legendPosition = "top",
                                                 xlab = "avg_Log2FC",
                                                 ylab = "-Log10(P.adj)") & 
  theme(text = element_text(family = "Arial"))

cd4tnfrsf9negtreg_DEG_volcano$layers[[4]]$aes_params$family <- "Arial" 

cairo_pdf("nebula/cd4t/20251006_escc_nebula_CD4_TNFRSF9negTreg_volcano_km.pdf", width = cm_to_inch(18), height = cm_to_inch(16), family = "Arial")
print(cd4tnfrsf9negtreg_DEG_volcano) & 
  theme(text = element_text(family = "Arial"))
dev.off()


#CD4 TNFRSF9+ Treg
re <- readRDS("nebula/cd4t/nebula_cd4tnfrsf9postreg_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoCD4 TNFRSF9+ Treg`, `se_nebula_annoCD4 TNFRSF9+ Treg`, `p_nebula_annoCD4 TNFRSF9+ Treg`)
res$Padj <- p.adjust(res$`p_nebula_annoCD4 TNFRSF9+ Treg`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
pos <- res %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% head(8) %>% pull(gene)
neg <- res %>% filter(p_val_adj < 0.05) %>% arrange(avg_log2FC) %>% head(8) %>% pull(gene)

cd4tnfrsf9postreg_DEG_volcano <- EnhancedVolcano(res, lab = res$gene, x = "avg_log2FC", y = "p_val_adj",
                                                 title = 'Other CD4 T cells vs. CD4 TNFRSF9+ Treg',
                                                 selectLab = c(neg,pos),
                                                 subtitle = NULL,
                                                 labSize = 6.0,
                                                 labCol = 'black',
                                                 pCutoff = 0.05,
                                                 FCcutoff = 1,
                                                 pointSize = 3.0,
                                                 drawConnectors = TRUE,
                                                 legendLabels = c("NS", "avg_Log2FC", "P.adj", "P.adj and avg_Log2FC"),
                                                 legendPosition = "top",
                                                 xlab = "avg_Log2FC",
                                                 ylab = "-Log10(P.adj)") & 
  theme(text = element_text(family = "Arial"))

cd4tnfrsf9postreg_DEG_volcano$layers[[4]]$aes_params$family <- "Arial" 

cairo_pdf("nebula/cd4t/20251006_escc_nebula_CD4_TNFRSF9posTreg_volcano_km.pdf", width = cm_to_inch(18), height = cm_to_inch(16), family = "Arial")
print(cd4tnfrsf9postreg_DEG_volcano) & 
  theme(text = element_text(family = "Arial"))
dev.off()


#Fig. 5SD####
escc_visium_C2F_latest_dl <- readRDS('~/Dropbox/escc_figure_share/dl/20250812_escc_visium_C2F_latest_dl.rds')
escc_visium_C2F_latest_dl@active.assay <- 'C2L'

for (i in 1:10){
  assign(paste0('p',i), SpatialFeaturePlot(escc_visium_C2F_latest_dl, c('Tumor', 'LA TAM', 'CD4 CTL','CD4 TNFRSF9- Treg','CD4 TNFRSF9+ Treg'),
                                           pt.size.factor = 1.6, images = paste0('s',i),image.alpha = 0,  ncol = 5) &
           theme(text = element_text(family = "Arial")) & NoLegend())
}
p1 <- SpatialFeaturePlot(escc_visium_C2F_latest_dl, c('Tumor', 'LA TAM', 'CD4 CTL','CD4 TNFRSF9- Treg','CD4 TNFRSF9+ Treg'),
                         pt.size.factor = 2.1, images = 's1',image.alpha = 0.1,  ncol = 5) & 
  theme(text = element_text(family = "Arial")) & NoLegend()

p6 <- SpatialFeaturePlot(escc_visium_C2F_latest_dl, c('Tumor', 'LA TAM', 'CD4 CTL','CD4 TNFRSF9- Treg','CD4 TNFRSF9+ Treg'),
                         pt.size.factor = 1.85, images = 's6',image.alpha = 0,  ncol = 5) & 
  theme(text = element_text(family = "Arial")) & NoLegend()

p7 <- SpatialFeaturePlot(escc_visium_C2F_latest_dl, c('Tumor', 'LA TAM', 'CD4 CTL','CD4 TNFRSF9- Treg','CD4 TNFRSF9+ Treg'),
                         pt.size.factor = 1.85, images = 's7',image.alpha = 0,  ncol = 5) & 
  theme(text = element_text(family = "Arial")) & NoLegend()

library(showtext)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cairo_pdf("~/Dropbox/escc_figure_share/FigS5/20250822/20250823_escc_FigS5_visium_CTL_Treg_km.pdf", width = 38, height = 90, family="Arial")
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 1)
dev.off()

cairo_pdf("~/Dropbox/escc_figure_share/FigS5/20250822/20250823_escc_FigS5_visium_CTL_Treg_image_km.pdf", width = 38, height = 90, family="Arial")
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 1)
dev.off()


#Fig. 5SG####
##volcano####
#LA TAM
re <- readRDS("nebula/myeloid/nebula_LAtam_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoLA TAM`, `se_nebula_annoLA TAM`, `p_nebula_annoLA TAM`)
res$Padj <- p.adjust(res$`p_nebula_annoLA TAM`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
pos <- res %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% head(8) %>% pull(gene)
neg <- res %>% filter(p_val_adj < 0.05) %>% arrange(avg_log2FC) %>% head(8) %>% pull(gene)

LA_TAM_DEG_volcano <- EnhancedVolcano(res, lab = res$gene, x = "avg_log2FC", y = "p_val_adj",
                                      title = 'C1QC+ TAM + Int. Mac vs. LA TAM',
                                      selectLab = c(neg,pos),
                                      subtitle = NULL,
                                      labSize = 6.0,
                                      labCol = 'black',
                                      pCutoff = 0.05,
                                      FCcutoff = 1,
                                      pointSize = 3.0,
                                      drawConnectors = TRUE,
                                      legendLabels = c("NS", "avg_Log2FC", "P.adj", "P.adj and avg_Log2FC"),
                                      legendPosition = "top",
                                      xlab = "avg_Log2FC",
                                      ylab = "-Log10(P.adj)") & 
  theme(text = element_text(family = "Arial"))

LA_TAM_DEG_volcano$layers[[4]]$aes_params$family <- "Arial" 

cairo_pdf("nebula/myeloid/20251006_escc_nebula_LA_TAM_volcano_km.pdf", width = cm_to_inch(18), height = cm_to_inch(16), family = "Arial")
print(LA_TAM_DEG_volcano) & 
  theme(text = element_text(family = "Arial"))
dev.off()


#C1QC+ TAM
re <- readRDS("nebula/myeloid/nebula_C1QCposTAM_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoC1QC+ TAM`, `se_nebula_annoC1QC+ TAM`, `p_nebula_annoC1QC+ TAM`)
res$Padj <- p.adjust(res$`p_nebula_annoC1QC+ TAM`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
pos <- res %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% head(5) %>% pull(gene)
neg <- res %>% filter(p_val_adj < 0.05) %>% arrange(avg_log2FC) %>% head(5) %>% pull(gene)

C1QC_TAM_DEG_volcano <- EnhancedVolcano(res, lab = res$gene, x = "avg_log2FC", y = "p_val_adj",
                                        title = 'LA TAM + Int. Mac vs. C1QC+ TAM',
                                        selectLab = c(neg,pos),
                                        subtitle = NULL,
                                        labSize = 6.0,
                                        labCol = 'black',
                                        pCutoff = 0.05,
                                        FCcutoff = 1,
                                        pointSize = 3.0,
                                        drawConnectors = TRUE,
                                        legendLabels = c("NS", "avg_Log2FC", "P.adj", "P.adj and avg_Log2FC"),
                                        legendPosition = "top",
                                        xlab = "avg_Log2FC",
                                        ylab = "-Log10(P.adj)") & 
  theme(text = element_text(family = "Arial"))

C1QC_TAM_DEG_volcano$layers[[4]]$aes_params$family <- "Arial" 

cairo_pdf("nebula/myeloid/20251006_escc_nebula_C1QCpos_TAM_volcano_km.pdf", width = cm_to_inch(18), height = cm_to_inch(16), family = "Arial")
print(C1QC_TAM_DEG_volcano) & 
  theme(text = element_text(family = "Arial"))
dev.off()


#Interstitial macrophages
re <- readRDS("nebula/myeloid/nebula_Intmac_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoInterstitial macrophages`, `se_nebula_annoInterstitial macrophages`, `p_nebula_annoInterstitial macrophages`)
res$Padj <- p.adjust(res$`p_nebula_annoInterstitial macrophages`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")
pos <- res %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% head(5) %>% pull(gene)
neg <- res %>% filter(p_val_adj < 0.05) %>% arrange(avg_log2FC) %>% head(5) %>% pull(gene)

Int_Macro_DEG_volcano <- EnhancedVolcano(res, lab = res$gene, x = "avg_log2FC", y = "p_val_adj",
                                         title = 'LA TAM + C1QC+ TAM vs. Int. Mac',
                                         selectLab = c(neg,pos),
                                         subtitle = NULL,
                                         labSize = 6.0,
                                         labCol = 'black',
                                         pCutoff = 0.05,
                                         FCcutoff = 1,
                                         pointSize = 3.0,
                                         drawConnectors = TRUE,
                                         legendLabels = c("NS", "avg_Log2FC", "P.adj", "P.adj and avg_Log2FC"),
                                         legendPosition = "top",
                                         xlab = "avg_Log2FC",
                                         ylab = "-Log10(P.adj)") & 
  theme(text = element_text(family = "Arial"))

Int_Macro_DEG_volcano$layers[[4]]$aes_params$family <- "Arial" 

cairo_pdf("nebula/myeloid/20251006_escc_nebula_Int_mac_volcano_km.pdf", width = cm_to_inch(18), height = cm_to_inch(16), family = "Arial")
print(Int_Macro_DEG_volcano) & 
  theme(text = element_text(family = "Arial"))
dev.off()



#Fig. S5H####
library(Seurat)
library(dplyr)
library(AUCell)
library(SCENIC)
library(ComplexHeatmap)
library(stringr)
library(showtext)
library(CellChat)
library(colorRamp2)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

# 1. wilcoxon ----

wilcoxon_res <- list()
regulonAUC <- readRDS('scenic/myeloid/int/3.4_regulonAUC.Rds')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

# 2. scaling ----
# sc ------
seurat_obj <- readRDS("R_objects/20250722_escc_myeloid_flt_km.rds")

seurat_obj <- subset(seurat_obj, subset = annotation_myeloid %in% c("LA TAM","C1QC+ TAM","Interstitial macrophages"))

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
  
  for (celltype_tmp in c(c(seurat_obj$annotation_myeloid %>% unique()))) {
    df1 <- df_tmp %>% dplyr::filter(annotation_myeloid == celltype_tmp)
    df2 <- df_tmp %>% dplyr::filter(annotation_myeloid != celltype_tmp)
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

saveRDS(res_df, "scenic/myeloid/downstream/20250810_escc_scenic_tam_sc_test_km.rds")
res_df <- readRDS("scenic/myeloid/downstream/20250810_escc_scenic_tam_sc_test_km.rds")


# sn ------
seurat_obj <- readRDS("R_objects/20250722_escc_myeloid_flt_km.rds")

seurat_obj <- subset(seurat_obj, subset = annotation_myeloid %in% c("LA TAM","C1QC+ TAM","Interstitial macrophages"))

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
  
  for (celltype_tmp in c(c(seurat_obj$annotation_myeloid %>% unique()))) {
    df1 <- df_tmp %>% dplyr::filter(annotation_myeloid == celltype_tmp)
    df2 <- df_tmp %>% dplyr::filter(annotation_myeloid != celltype_tmp)
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

saveRDS(res_df, "scenic/myeloid/downstream/20250810_escc_scenic_tam_sn_test_km.rds")
res_df <- readRDS("scenic/myeloid/downstream/20250810_escc_scenic_tam_sn_test_km.rds")

res_sc_df <- readRDS("scenic/myeloid/downstream/20250810_escc_scenic_tam_sc_test_km.rds")
res_sc_df$platform <- "sc"
res_sn_df <- readRDS("scenic/myeloid/downstream/20250810_escc_scenic_tam_sn_test_km.rds")
res_sn_df$platform <- "sn"

res_df <- rbind(res_sc_df, res_sn_df)
quant <- quantile(res_df$padj, 0.1)

filtered_df <- res_df %>%
  dplyr::group_by(group1, platform) %>%
  dplyr::filter(padj < quantile(padj, 0.15, na.rm = TRUE)) %>%
  ungroup()

res_df <- filtered_df %>% dplyr::mutate(log2FC = log2(group1_mean/group2_mean)) %>% 
  dplyr::filter(group2_mean < group1_mean, padj < 0.05) %>% dplyr::filter(log2FC > 1)

#res_counts <- res_df %>% dplyr::group_by(TF, group1) %>% dplyr::summarise(n=n()) %>% dplyr::arrange(desc(n))
#mytfs <- res_counts %>% dplyr::filter(n == 2) %>% dplyr::pull(TF) %>% unique()
mytfs <- res_df %>% dplyr::pull(TF) %>% unique()
mytfs <- mytfs[-30] #STAT2 (78g) 제거.
#mytfs <- res_df %>% group_by(group1) %>%
slice_min(order_by = padj, n = 5, with_ties = FALSE) %>%
  ungroup() %>% pull(TF)

# 3. merge sc, sn ----

seurat_obj <- readRDS("R_objects/20250722_escc_myeloid_flt_km.rds")

seurat_obj <- subset(seurat_obj, subset = annotation_myeloid %in% c("LA TAM","C1QC+ TAM","Interstitial macrophages"))

seurat_obj <- subset(seurat_obj, subset = platform == "scRNA")
regulonAUC_scaled_sc <- getAUC(regulonAUC)[mytfs, colnames(seurat_obj)]

regulonAUC_scaled_sc <- t(scale(t(regulonAUC_scaled_sc), center = T, scale=T))
# regulonAUC_scaled_sc <- t(apply(regulonAUC_scaled_sc, 1, function(x) {
#   (x - min(x)) / (max(x) - min(x))
# }))
#regulonAUC_scaled_sc <- t(scale(t(getAUC(regulonAUC)[res_df_filtered$TF, colnames(seurat_obj)]), center = T, scale=T))

seurat_obj <- readRDS("R_objects/20250722_escc_myeloid_flt_km.rds")


seurat_obj <- subset(seurat_obj, subset = annotation_myeloid %in% c("LA TAM","C1QC+ TAM","Interstitial macrophages"))

seurat_obj <- subset(seurat_obj, subset = platform == "snRNA")
regulonAUC_scaled_sn <- getAUC(regulonAUC)[mytfs, colnames(seurat_obj)]
regulonAUC_scaled_sn <- t(scale(t(regulonAUC_scaled_sn), center = T, scale=T))
regulonAUC_scaled <- cbind(regulonAUC_scaled_sn, regulonAUC_scaled_sc)

seurat_obj <- readRDS("R_objects/20250722_escc_myeloid_flt_km.rds")

seurat_obj <- subset(seurat_obj, subset = annotation_myeloid %in% c("LA TAM","C1QC+ TAM","Interstitial macrophages"))

celltypes <- seurat_obj@meta.data %>% select(annotation_myeloid)
celltypes$New_Identity <- factor(celltypes$annotation_myeloid, levels = c("Interstitial macrophages", "C1QC+ TAM", "LA TAM"))
selectedResolution <- "New_Identity"
cellsPerTypes <- split(rownames(celltypes), celltypes[,selectedResolution]) 

# 4. scenic heatmap -----

cellorder <- unlist(cellsPerTypes)

cell_type <- data.frame(
  "cell type" = unlist(mapply(rep, names(cellsPerTypes), lapply(cellsPerTypes, length))))

regulonAUC_scaled_order <- regulonAUC_scaled[mytfs, cellorder]
row.names(cell_type) <- colnames(regulonAUC_scaled_order)
palette_length = 100

date <- "20250823"
project_name <- "escc"
cm_to_inch <- function(cm) cm / 2.54
fileidentity <- "tam_scenic"

library(CellChat)
mycol <- scPalette(3)
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
                       annotation_name_gp = gpar(fontsize = 8),
                       annotation_legend_param = list(
                         title_gp = gpar(fontsize = 8, fontfamily = "Arial"),
                         labels_gp = gpar(fontsize = 8, fontfamily = "Arial")
                       ))

col_fun = colorRamp2(c(-1, 0, 1), c("darkblue", "white", "red"))

# png file 
# png파일 만든 이유는 heatmap 부분을 raster 를 주면 상단에 있는 annotation부분까지 raster 적용되어서 색상이 흐릿해지더라고.
# 그래서 png 파일에서 heatmap, dendrogram, annotation부분까지 캡쳐해서 잉크스케이프에 넣고 (선명해)
png(filename = paste0("scenic/myeloid/downstream/",date, "_", project_name, "_", fileidentity, "_km.png"),
    width = 8, height = 3.3, units = "in", res = 600, family = "Arial")
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
                             column_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
                             heatmap_legend_param = list(
                               title = "Score",
                               title_gp = gpar(fontsize = 8, fontfamily = "Arial"),
                               labels_gp = gpar(fontsize = 8, fontfamily = "Arial")
                             ),
                             top_annotation = ha)
draw(ht, use_raster = TRUE)
dev.off()


# 그 다음에 글자 부분은 다시 pdf에서 따왔어! 여기서 raster 를 꼭 지정해줘야 잉크스케이프에서 열릴꺼야.

cairo_pdf(paste0("scenic/myeloid/downstream/",date, "_", project_name, "_", fileidentity, "_km.pdf"),
          width = 8, height = 3.3, family = "Arial")

col_fun <- colorRamp2(c(-1, 0, 1), c("darkblue", "white", "red"))

ht <- Heatmap(
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
  row_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
  column_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
  heatmap_legend_param = list(
    title = "Score",
    title_gp = gpar(fontsize = 8, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 8, fontfamily = "Arial")
  ),
  top_annotation = ha
)

draw(ht)
dev.off()


# 5. ORA ----
library(org.Hs.eg.db)
library(msigdbr)
library(clusterProfiler)

regulons <- readRDS('scenic/myeloid/int/2.6_regulons_asGeneSet.Rds')
#GOBP
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, entrez_gene)
#Hallmark
m_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

hs <- org.Hs.eg.db
row_order_vec <- row_order(ht)
ordered_row_names <- rownames(regulonAUC_scaled_order)[row_order_vec]

tfs <- gsub(ordered_row_names, pattern = " \\([0-9]*g\\)", replacement = "")

# for (i in 1:length(tfs)) {
#   if ("COL1A1" %in% regulons[[tfs[i]]]) {
#     print(tfs[i])
#   }
# }

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

library(openxlsx)
write.xlsx(enricher_res, "scenic/myeloid/downstream/20250812_escc_tam_scenic_targetgenes_gobp_km.xlsx") #v1 top5 each celltype
write.xlsx(enricher_res, "scenic/myeloid/downstream/20250823_escc_tam_scenic_targetgenes_gobp_km.xlsx") #v2 log2FC > 1

enricher_res <- read.xlsx("scenic/myeloid/downstream/20250812_escc_tam_scenic_targetgenes_gobp_km.xlsx") #v1 top5 each celltype
enricher_res <- read.xlsx("scenic/myeloid/downstream/20250823_escc_tam_scenic_targetgenes_gobp_km.xlsx") #v2 log2FC > 1
enricher_res <- enricher_res %>% dplyr::filter(p.adjust < 0.05)
write.xlsx(enricher_res, "scenic/myeloid/downstream/20250812_escc_tam_scenic_targetgenes_gobp_sig_km.xlsx") #v1 top5 each celltype
write.xlsx(enricher_res, "scenic/myeloid/downstream/20250823_escc_tam_scenic_targetgenes_gobp_sig_km.xlsx") #v2 log2FC > 1
enricher_res <- enricher_res %>% dplyr::group_by(TFs) %>% dplyr::arrange(p.adjust, .by_group = T) %>% dplyr::slice(1)
rownames(enricher_res) <- enricher_res$TFs

tf_filtered <- intersect(tfs, enricher_res$TFs)
enricher_res <- enricher_res[tf_filtered,] # sorting top1 results

write.xlsx(enricher_res, "scenic/myeloid/downstream/20250810_escc_tam_scenic_targetgenes_gobp_top1_km.xlsx") #v1 top5 each celltype
write.xlsx(enricher_res, "scenic/myeloid/downstream/20250823_escc_tam_scenic_targetgenes_gobp_top1_km.xlsx") #v2 log2FC > 1
enricher_res <- read.xlsx("scenic/myeloid/downstream/20250810_escc_tam_scenic_targetgenes_gobp_top1_km.xlsx") #v1 top5 each celltype
enricher_res <- read.xlsx("scenic/myeloid/downstream/20250823_escc_tam_scenic_targetgenes_gobp_top1_km.xlsx") #v2 log2FC > 1
# 6. GOBP plot ----
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

rownames(mat) <- sub("^GOBP_", "", rownames(mat))           # 1) GOBP_ 제거
rownames(mat) <- gsub("_", " ", rownames(mat))             # 2) 언더바를 공백으로
rownames(mat) <- str_to_lower(rownames(mat))              # 모두 소문자
rownames(mat) <- str_to_sentence(rownames(mat))
row_lbl <- rownames(mat)

rng <- range(val, na.rm = TRUE)
if (use_log) {
  col_fun <- colorRamp2(c(0, rng[2]), c("white", "red"))
} else {
  col_fun <- colorRamp2(c(rng[1], rng[2]), c("white", "red"))
}


ht <- Heatmap(
  mat,
  name = val_name,
  col  = col_fun,
  na_col = "darkgrey",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_labels = row_lbl,
  row_names_gp = gpar(fontsize = 9, fontfamily = "Arial"),
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(col = "black", fill = NA, lwd = 0.5))
  },
  width = unit(ncol(mat)/4, "cm"), 
  height = unit(nrow(mat)/4, "cm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 9, fontfamily = "Arial"),
    legend_height = unit(1, "cm"),                        
    grid_height = unit(0.2, "cm"),                        
    grid_width = unit(0.1, "cm")                          
  )
)

# save plots
fileidentity <- "scenic_top1_gobp"
cairo_pdf(paste0("scenic/myeloid/downstream/",date, "_", project_name, "_", fileidentity, "_km.pdf"),
          width = cm_to_inch(20), height = cm_to_inch(8), family = "Arial")
draw(
  ht,
  heatmap_legend_side = "left"
)
dev.off()


