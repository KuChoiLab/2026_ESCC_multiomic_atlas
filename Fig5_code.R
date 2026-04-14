# Fig 5A ----
setwd("~/Dropbox/ESCC_snrna/202507")
library(ggalluvial)
library(Seurat)
library(dplyr)
library(cowplot)

##myeloid cells
myeloid <- readRDS('R_objects/20250722_escc_myeloid_flt_km.rds')
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
fib <- readRDS('R_objects/20250728_escc_Fib_fDC_FRC_hj_km.rds')

Idents(fib) <- 'pathology'
table_fib <- fib@meta.data %>% select(pathology, Annotation_v2) %>% table()
table_fib_v2 <- ( table_fib / table_fib %>% apply(1,sum) ) %>% data.frame()
table_fib_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',table_fib_v2$pathology)
table_fib_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',table_fib_v2$pathology)
table_fib_v2$pathology <- factor(table_fib_v2$pathology,  levels = c('Normal', 'Dysplasia', 'Microinvasive carcinoma', 'Macroinvasive carcinoma'))
table_fib_v2$Annotation_v2 %>% table
rm(fib)


##CD8T
cd8t <- readRDS("R_objects/20250728_escc_cd8_annotation_dl.rds")

Idents(cd8t) <- 'pathology'
table_cd8t <- cd8t@meta.data %>% select(pathology, Annotation_v2) %>% table()
table_cd8t_v2 <- ( table_cd8t / table_cd8t %>% apply(1,sum) ) %>% data.frame()
table_cd8t_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',table_cd8t_v2$pathology)
table_cd8t_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',table_cd8t_v2$pathology)
table_cd8t_v2$pathology <- factor(table_cd8t_v2$pathology,  levels = c('Normal', 'Dysplasia', 'Microinvasive carcinoma', 'Macroinvasive carcinoma'))
table_cd8t_v2$Annotation_v2 %>% table
rm(cd8t)


##CD4T
cd4t <- readRDS("R_objects/20250729_escc_cd4_celltype_annotated_dl.rds")

Idents(cd4t) <- 'pathology'
table_cd4t <- cd4t@meta.data %>% select(pathology, Annotation_v2) %>% table()
table_cd4t_v2 <- ( table_cd4t / table_cd4t %>% apply(1,sum) ) %>% data.frame()
table_cd4t_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',table_cd4t_v2$pathology)
table_cd4t_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',table_cd4t_v2$pathology)
table_cd4t_v2$pathology <- factor(table_cd4t_v2$pathology,  levels = c('Normal', 'Dysplasia', 'Microinvasive carcinoma', 'Macroinvasive carcinoma'))
table_cd4t_v2$Annotation_v2 %>% table
rm(cd4t)


##B cells
b_cells <- readRDS("R_objects/20250723_escc_B_GCB_Plasma_flt_km_v2.rds")
b_cells$annotation_B <- gsub('naive','Naive',b_cells$annotation_B)

Idents(b_cells) <- 'pathology'
table_B <- b_cells@meta.data %>% select(pathology, annotation_B) %>% table()
table_B_v2 <- ( table_B / table_B %>% apply(1,sum) ) %>% data.frame()
table_B_v2$pathology <- gsub('Microinvasive','Microinvasive carcinoma',table_B_v2$pathology)
table_B_v2$pathology <- gsub('Macroinvasive','Macroinvasive carcinoma',table_B_v2$pathology)
table_B_v2$pathology <- factor(table_B_v2$pathology,  levels = c('Normal', 'Dysplasia', 'Microinvasive carcinoma', 'Macroinvasive carcinoma'))
table_B_v2$annotation_B %>% table
rm(b_cells)


CD4_Alluvial_legend <- get_legend(ggplot(table_cd4t_v2, aes( x = pathology, y = Freq, alluvium = Annotation_v2)) + 
                                    geom_alluvium(aes(fill = Annotation_v2), colour = "black", alpha = 0.9, decreasing = FALSE) + 
                                    labs(y = "Proportion of Cells", fill = "", colour = "") + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + scale_x_discrete(expand = c(0,0)) + 
                                    theme(text = element_text(family = "Arial"),legend.position = "bottom", legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"), axis.text.y = element_text(colour = "black", size = 10),legend.text = element_text(size = 12),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks = element_line(color = "black", size=0.3),panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+ guides(fill = guide_legend(nrow=3))) 
CD4_Alluvial_Plots <- ggplot(table_cd4t_v2, aes( x = pathology, y = Freq, alluvium = Annotation_v2)) + 
  geom_alluvium(aes(fill = Annotation_v2), colour = "black", alpha = 0.9, decreasing = FALSE) + 
  labs(y = "Proportion of Cells", fill = "", colour = "") + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + scale_x_discrete(expand = c(0,0))+ theme(text = element_text(family = "Arial"), legend.position="none", axis.text.y = element_text(colour = "black", size = 10),legend.text = element_text(size = 12), axis.text.x = element_text(colour = "black", angle = 15, size = 12, vjust = 0.9, hjust = 0.95), axis.title.y =element_blank(),axis.title.x = element_blank(),axis.ticks = element_line(color = "black", size=0.3), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))
CD8_Alluvial_legend <- get_legend(ggplot(table_cd8t_v2, aes( x = pathology, y = Freq, alluvium = Annotation_v2)) + 
                                    geom_alluvium(aes(fill = Annotation_v2), colour = "black", alpha = 0.9, decreasing = FALSE) + 
                                    labs(y = "Proportion of Cells", fill = "", colour = "") + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + scale_x_discrete(expand = c(0,0)) + 
                                    theme(text = element_text(family = "Arial"), legend.position = "bottom", legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"), axis.text.y = element_text(colour = "black", size = 10),legend.text = element_text(size = 12),axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks = element_line(color = "black", size=0.3),panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+ guides(fill = guide_legend(nrow=4))) 
CD8_Alluvial_Plots <- ggplot(table_cd8t_v2, aes( x = pathology, y = Freq, alluvium = Annotation_v2)) + 
  geom_alluvium(aes(fill = Annotation_v2), colour = "black", alpha = 0.9, decreasing = FALSE) + 
  labs(y = "Proportion of Cells", fill = "", colour = "") + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + scale_x_discrete(expand = c(0,0))+ theme(text = element_text(family = "Arial"), legend.position="none", axis.text.y = element_text(colour = "black", size = 10),legend.text = element_text(size = 12), axis.text.x = element_text(colour = "black", angle = 15, size = 12, vjust = 0.9, hjust = 0.95), axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks = element_line(color = "black", size=0.3), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))
B_Alluvial_legend <- get_legend(ggplot(table_B_v2, aes(x = pathology, y = Freq, alluvium = annotation_B)) + 
                                  geom_alluvium(aes(fill = annotation_B), colour = "black", alpha = 0.9, decreasing = FALSE) + 
                                  labs(y = "Proportion of Cells", fill = "", colour = "") + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + scale_x_discrete(expand = c(0,0)) + 
                                  theme(text = element_text(family = "Arial"), legend.position = "bottom", legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"), axis.text.y = element_text(colour = "black", size = 10), legend.text = element_text(size = 12), axis.text.x = element_text(colour = "black", angle = 15, size = 12, vjust = 0.9, hjust = 0.95), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks = element_line(color = "black", size=0.3), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+ guides(fill = guide_legend(nrow=3)))
B_Alluvial_Plots <- ggplot(table_B_v2, aes(x = pathology, y = Freq, alluvium = annotation_B)) + 
  geom_alluvium(aes(fill = annotation_B), colour = "black", alpha = 0.9, decreasing = FALSE) + labs(y = "Proportion of Cells", fill = "", colour = "") + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + scale_x_discrete(expand = c(0,0))+ theme(text = element_text(family = "Arial"), legend.position="none", axis.text.y = element_text(colour = "black", size = 10),legend.text = element_text(size = 12), axis.text.x = element_text(colour = "black", angle = 15, size = 12, vjust = 0.9, hjust = 0.95), axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks = element_line(color = "black", size=0.3), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

Stromal_Alluvial_legend <- get_legend(ggplot(table_fib_v2, aes( x = pathology, y = Freq, alluvium = Annotation_v2)) + 
                                        geom_alluvium(aes(fill = Annotation_v2), colour = "black", alpha = 0.9, decreasing = FALSE) + 
                                        labs(y = "Proportion of Cells", fill = "", colour = "") + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + scale_x_discrete(expand = c(0,0)) + 
                                        theme(text = element_text(family = "Arial"), legend.position = "bottom", legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"), axis.text.y = element_text(colour = "black", size = 10),legend.text = element_text(size = 12),axis.text.x = element_text(colour = "black", angle = 15, size = 10, vjust = 0.9, hjust = 0.95), axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks = element_line(color = "black", size=0.3),panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+ guides(fill = guide_legend(nrow=3))) 
Stromal_Alluvial_Plots <- ggplot(table_fib_v2, aes( x = pathology, y = Freq, alluvium = Annotation_v2)) + 
  geom_alluvium(aes(fill = Annotation_v2), colour = "black", alpha = 0.9, decreasing = FALSE) + 
  labs(y = "Proportion of Cells", fill = "", colour = "") + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + scale_x_discrete(expand = c(0,0))+ theme(text = element_text(family = "Arial"), legend.position="none", axis.text.y = element_text(colour = "black", size = 10),legend.text = element_text(size = 12), axis.text.x = element_text(colour = "black", angle = 15, size = 12, vjust = 0.9, hjust = 0.95), axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks = element_line(color = "black", size=0.3), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))
Myeloid_Alluvial_legend <- get_legend(ggplot(table_myeloid_v2, aes( x = pathology, y = Freq, alluvium = annotation_myeloid)) + 
                                        geom_alluvium(aes(fill = annotation_myeloid), colour = "black", alpha = 0.9, decreasing = FALSE) + labs(y = "Proportion of Cells", fill = "", colour = "") + 
                                        scale_y_continuous(expand = c(0,0), limits = c(0,1)) + scale_x_discrete(expand = c(0,0)) + 
                                        theme(text = element_text(family = "Arial"), legend.position = "bottom", legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"), axis.text.y = element_text(colour = "black", size = 10),legend.text = element_text(size = 12), axis.text.x = element_text(colour = "black", angle = 15, size = 12, vjust = 0.9, hjust = 0.95), axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks = element_line(color = "black", size=0.3), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+ guides(fill = guide_legend(nrow=5))) 
Myeloid_Alluvial_plots <- ggplot(table_myeloid_v2, aes( x = pathology, y = Freq, alluvium = annotation_myeloid)) + 
  geom_alluvium(aes(fill = annotation_myeloid), colour = "black", alpha = 0.9, decreasing = FALSE) + 
  labs(y = "Proportion of Cells", fill = "", colour = "") + scale_y_continuous(expand = c(0,0), limits = c(0,1)) + scale_x_discrete(expand = c(0,0))+ theme(text = element_text(family = "Arial"), legend.position="none", axis.text.y = element_text(colour = "black", size = 10),legend.text = element_text(size = 12), axis.text.x = element_text(colour = "black", angle = 15, size = 12, vjust = 0.9, hjust = 0.95), axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks = element_line(color = "black", size=0.3), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

library(showtext)

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

cairo_pdf("~/Dropbox/escc_figure_share/Fig5/20250822/20250828_escc_Fig5A_abundance_change_alluvial_km.pdf", width = 21, height = 6, family = "Arial")
plot_grid(CD4_Alluvial_Plots/CD4_Alluvial_legend,NULL,
          CD8_Alluvial_Plots/CD8_Alluvial_legend,NULL,
          B_Alluvial_Plots/B_Alluvial_legend,NULL,
          Stromal_Alluvial_Plots/Stromal_Alluvial_legend,NULL,
          Myeloid_Alluvial_plots/Myeloid_Alluvial_legend, ncol=9,
          rel_widths = c(1,0.2, 1, 0.2,1, 0.2,1, 0.2,1))
dev.off()

cairo_pdf("~/Dropbox/escc_figure_share/Fig5/20250822/20250828_escc_Fig5A_abundance_change_alluvial_km_legend.pdf", width = 41, height = 6, family = "Arial") #figure 가져오는 용.
plot_grid(CD4_Alluvial_Plots/CD4_Alluvial_legend,NULL,
          CD8_Alluvial_Plots/CD8_Alluvial_legend,NULL,
          B_Alluvial_Plots/B_Alluvial_legend,NULL,
          Stromal_Alluvial_Plots/Stromal_Alluvial_legend,NULL,
          Myeloid_Alluvial_plots/Myeloid_Alluvial_legend, ncol=9,
          rel_widths = c(1,0.2, 1, 0.2,1, 0.2,1, 0.2,1))
dev.off()


# Fig 5C ----
#nebula
re <- readRDS("nebula/cd4t/nebula_cd4ctl_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoCD4 CTL`, `se_nebula_annoCD4 CTL`, `p_nebula_annoCD4 CTL`)
res$Padj <- p.adjust(res$`p_nebula_annoCD4 CTL`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "p_val_adj")

#plot
cairo_pdf("~/Dropbox/escc_figure_share/Fig5/20250822/20251006_escc_Fig5c_CD4CTL_volcano_km.pdf", width = cm_to_inch(12), height = cm_to_inch(16), family = "Arial")
EnhancedVolcano(res, lab = res$gene, x = "avg_log2FC", y = "p_val_adj",
                title = 'Remaining CD4 T vs. CD4 CTL',
                selectLab = c('CCL5','GZMA','GZMB','GZMK','GZMH','KLRB1','PRF1'),
                labSize = 5.0,
                labCol = 'black',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                drawConnectors = TRUE,
                max.overlaps = Inf,
                arrowheads = TRUE,
                legendLabels = c("NS", "avg_Log2FC", "P.adj", "P.adj and avg_Log2FC"),
                legendPosition = "top",
                xlab = "avg_Log2FC",
                ylab = "-Log10(P.adj)") & theme(text = element_text(family = "Arial"))
dev.off()

#legend
cairo_pdf("~/Dropbox/escc_figure_share/Fig5/20250822/20251006_escc_Fig5c_CD4CTL_volcano_legned_km.pdf", width = cm_to_inch(20), height = cm_to_inch(16), family = "Arial")
EnhancedVolcano(res, lab = res$gene, x = "avg_log2FC", y = "p_val_adj",
                title = 'Remaining CD4 T vs. CD4 CTL',
                selectLab = c('CCL5','GZMA','GZMB','GZMK','GZMH','KLRB1','PRF1'),
                labSize = 5.0,
                labCol = 'black',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                drawConnectors = TRUE,
                max.overlaps = Inf,
                arrowheads = TRUE,
                legendLabels = c("NS", "avg_Log2FC", "P.adj", "P.adj and avg_Log2FC"),
                legendPosition = "top",
                xlab = "avg_Log2FC",
                ylab = "-Log10(P.adj)") & theme(text = element_text(family = "Arial"))
dev.off()


# Fig 5D ----


## GSEA ####

library(fgsea)
library(tibble)

mygsea_res <- list()

# CD4 CTL
celltype <- "CD4 CTL"
re <- readRDS("nebula/cd4t/nebula_cd4ctl_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoCD4 CTL`, `se_nebula_annoCD4 CTL`, `p_nebula_annoCD4 CTL`)
res$Padj <- p.adjust(res$`p_nebula_annoCD4 CTL`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "Padj")

df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)

df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("~/Dropbox/database/c5.go.bp.v2025.1.Hs.symbols.gmt")
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

# CD4 TNFRSF9+ Treg
celltype <- "CD4 TNFRSF9+ Treg"
re <- readRDS("nebula/cd4t/nebula_cd4tnfrsf9postreg_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoCD4 TNFRSF9+ Treg`, `se_nebula_annoCD4 TNFRSF9+ Treg`, `p_nebula_annoCD4 TNFRSF9+ Treg`)
res$Padj <- p.adjust(res$`p_nebula_annoCD4 TNFRSF9+ Treg`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "Padj")

df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)

df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("~/Dropbox/database/c5.go.bp.v2025.1.Hs.symbols.gmt")
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


# CD4 TNFRSF9- Treg
celltype <- "CD4 TNFRSF9- Treg"
re <- readRDS("nebula/cd4t/nebula_cd4tnfrsf9negtreg_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoCD4 TNFRSF9- Treg`, `se_nebula_annoCD4 TNFRSF9- Treg`, `p_nebula_annoCD4 TNFRSF9- Treg`)
res$Padj <- p.adjust(res$`p_nebula_annoCD4 TNFRSF9- Treg`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "Padj")

df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)

df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("~/Dropbox/database/c5.go.bp.v2025.1.Hs.symbols.gmt")
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

mygsea_res

saveRDS(mygsea_res, "nebula/cd4t/20251005_escc_cd4t_gsea_gobp_km.rds")
mygsea_res <- readRDS("nebula/cd4t/20251005_escc_cd4t_gsea_gobp_km.rds")
write.xlsx(mygsea_res, 'nebula/cd4t/20251005_escc_cd4t_gsea_gobp_km.xlsx')


library(stringr)
mygsea_res <- mygsea_res %>%
  mutate(pathway = pathway %>%
           # Remove 'GOPB_' at the start
           str_remove("^GOBP_") %>%
           # Replace '_' with space
           str_replace_all("_", " ") %>%
           # Convert to lowercase
           str_to_lower() %>%
           # Capitalize first letter of each row (string)
           str_replace("^.", toupper))

mygsea_res_cd4ctl_sig <- mygsea_res %>% filter(celltype == 'CD4 CTL', padj < 0.05, NES > 0)
mygsea_res_cd4tregneg_sig <- mygsea_res %>% filter(celltype == 'CD4 TNFRSF9- Treg', padj < 0.05, NES > 0)
mygsea_res_cd4tregpos_sig <- mygsea_res %>% filter(celltype == 'CD4 TNFRSF9+ Treg', padj < 0.05, NES > 0)

specific_pathways <- mygsea_res %>%
  filter(padj < 0.05)
group_by(pathway) %>%
  filter(all(celltype == "CD4 CTL")) %>%
  pull(pathway) %>%
  unique()

# Step 2: Filter based on padj, NES, and those specific pathways
cd4ctl_pathways <- mygsea_res %>%
  filter(celltype == "CD4 CTL",
         padj < 0.05,
         NES > 0,
         pathway %in% specific_pathways)

#dotplot
CTL_gobp_sig_term <- c('Granulocyte migration',
                       'Leukocyte migration',
                       'Antimicrobial humoral immune response mediated by antimicrobial peptide',
                       'Regulation of inflammatory response',
                       'T cell mediated cytotoxicity')
TNFRSF9pos_Treg_gobp_sig_term <- c('Negative regulation of immune response',
                                   'T helper 1 type immune response',
                                   'Negative regulation of cytokine production involved in immune response')
TNFRSF9neg_Treg_gobp_sig_term <- c('Positive regulation of defense response',
                                   'Regulation of adaptive immune response',
                                   'Tumor necrosis factor mediated signaling pathway',
                                   'T cell mediated immunity',
                                   'T cell cytokine production')

all_gobp_sig_term <- unique(c(CTL_gobp_sig_term,TNFRSF9pos_Treg_gobp_sig_term,TNFRSF9neg_Treg_gobp_sig_term))

mygsea_res_plot <- mygsea_res %>% filter(pathway %in% all_gobp_sig_term, padj < 0.05)
mygsea_res_plot$celltype <- factor(mygsea_res_plot$celltype, levels = c('CD4 CTL', 'CD4 TNFRSF9- Treg','CD4 TNFRSF9+ Treg'))
mygsea_res_plot$pathway <- factor(mygsea_res_plot$pathway, levels = unique(c('Antimicrobial humoral immune response mediated by antimicrobial peptide',
                                                                             'Granulocyte migration',
                                                                             'Leukocyte migration',
                                                                             'Negative regulation of immune response',
                                                                             'T cell mediated cytotoxicity',
                                                                             'T helper 1 type immune response',
                                                                             'Regulation of adaptive immune response',
                                                                             'Regulation of inflammatory response',
                                                                             'T cell mediated immunity',
                                                                             'Positive regulation of defense response',
                                                                             'Tumor necrosis factor mediated signaling pathway',
                                                                             'T cell cytokine production',
                                                                             'Negative regulation of cytokine production involved in immune response')))

library(ggplot2)
library(showtext)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

cairo_pdf("nebula/cd4t/20251005_escc_sc.snrna_CD4T_nebula_GSEA_gobp_sig_km.pdf", width = cm_to_inch(32), height = cm_to_inch(21), family = "Arial")
ggplot(mygsea_res_plot, aes(celltype,pathway)) + geom_point(aes(color = NES, size = -log10(padj))) + 
  ylab(NULL) + theme_minimal(base_family = 'Arial') +
  theme(axis.text.x = element_text(colour = "black", angle = 45, size = 20, vjust=1, hjust=1), 
        axis.text.y = element_text(colour = "black", size = 20), plot.title = element_blank(), 
        legend.position = "right", legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black", size=0.5),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.spacing = unit(1, 'cm'), legend.spacing.x = unit(0.5,'cm')) + 
  guides(size = guide_legend(title = expression(-log[10] ~ "(padj)"))) + 
  #scale_size(range = c(2,6), breaks = c(10,20,30)) +
  scale_color_gradientn(colours=c("orange","red")) +
  scale_shape_manual(values = c(13,16))
dev.off()


# Fig 5E ----
setwd('~/Dropbox/ESCC_snrna/202507/')
library(Seurat)
library(dplyr)
library(CellChat)
library(ggplot2)
library(showtext)

myeloid <- readRDS('R_objects/20250722_escc_myeloid_flt_km.rds')

FeaturePlot(myeloid, c('MMP12','FABP4','FABP5','LIPA'), cols = c('grey','red'))
VlnPlot(myeloid, c('MMP12','FABP4','FABP5','LIPA'), group.by = 'annotation_myeloid', ncol = 2)
DimPlot(myeloid, group.by = 'annotation_myeloid')


setwd("~/Dropbox/ESCC_snrna/202507/")

font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cm_to_inch <- function(cm) cm / 2.54

date <- "20250806"
project_name <- "escc"

cellchat_normal <- readRDS(file = "cellchat/cd4t/R_objects/20250806_escc_normal_hj.rds")
cellchat_dysp <- readRDS(file = "cellchat/cd4t/R_objects/20250806_escc_dysplasia_hj.rds")
cellchat_micro <- readRDS(file = "cellchat/cd4t/R_objects/20250806_escc_microinvasive_carcinoma_hj.rds")
cellchat_macro <- readRDS(file = "cellchat/cd4t/R_objects/20250806_escc_macroinvasive_carcinoma_hj.rds")
object.list <- list(Normal = cellchat_normal,
                    Dysplasia = cellchat_dysp,
                    Micro_invasive = cellchat_micro,
                    Macro_invasive = cellchat_macro)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

##comparison between cancer stages####
library(showtext)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cairo_pdf("~/Dropbox/escc_figure_share/Fig5/20250822/20250822_escc_Fig5_Treg_CTL_cellchat_km.pdf", width = 9, height = 9, family="Arial")
netVisual_bubble(cellchat, sources.use = c('C1QC+ TAM','LA TAM','Interstitial macrophages'), 
                 targets.use = c('Tumor'), comparison = c(1,2,3,4)) +
  theme_classic(base_family = 'Arial') +
  theme(axis.text.x = element_text(size = 18, colour = c('#83f52c', '#f3f315', '#ff6600', '#ff0099'), angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18, colour = 'black'),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))
dev.off()


# Fig 5F ----
escc_visium_C2F_latest_dl <- readRDS('~/Dropbox/escc_figure_share/dl/20250812_escc_visium_C2F_latest_dl.rds')
escc_visium_C2F_latest_dl@active.assay <- 'C2L'
p7_main <- SpatialFeaturePlot(escc_visium_C2F_latest_dl, c('Tumor', 'CD4 CTL','CD4 TNFRSF9- Treg','CD4 TNFRSF9+ Treg'),
                              pt.size.factor = 1.85, images = 's7',image.alpha = 0,  ncol = 4) & 
  theme(text = element_text(family = "Arial")) & scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))

library(showtext)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cairo_pdf("~/Dropbox/escc_figure_share/Fig5/20250822/20250823_escc_Fig5_visium_CTL_Treg_km.pdf", width = 30, height = 9, family="Arial")
plot_grid(p7_main, ncol = 1)
dev.off()

# Fig 5G ----
#DEGs for each subtypes were identified in the script 20250813_escc_scrna_myeloid_GSEA_km.R
C1QC_TAM_DEG <- read.xlsx('GSEA/myeloid/20250813_escc_scrna_C1QC_TAM_wilcox_DEG_km.xlsx', rowNames = T)
Int_Macro_DEG <- read.xlsx('GSEA/myeloid/20250813_escc_scrna_Int_Macro_DEG_wilcox_DEG_km.xlsx', rowNames = T)
LA_TAM_DEG <- read.xlsx('GSEA/myeloid/20250813_escc_scrna_LA_TAM_wilcox_DEG_km.xlsx', rowNames = T)

intersect(macro_deg_sig %>% filter(Macrophages.Subpopulations == 'RTM_Int') %>% pull(gene), Int_Macro_DEG %>% arrange(desc(avg_log2FC)) %>% pull(Gene) %>% head(30)) # "F13A1" "FOLR2" "LYVE1" "FOSB" 
intersect(macro_deg_sig %>% filter(Macrophages.Subpopulations == 'Mac_LA') %>% pull(gene), LA_TAM_DEG %>% arrange(desc(avg_log2FC)) %>% pull(Gene) %>% head(40)) # "SPP1" "GPNMB" "APOE" "MMP12"
C1QC_TAM_DEG %>% arrange(desc(avg_log2FC)) %>% pull(Gene) %>% head(40) # "PRKCA" "GAPDH" "AC008691.1" "C1QC"

FeaturePlot(myeloid, c("C1QC","SPP1", "GPNMB", "APOE", "MMP12"), cols = c('grey','red'))

markers <- c("PRKCA", "GAPDH", "AC008691.1", #C1QC TAM
             "SPP1", "GPNMB", "APOE", "MMP12", #LA TAM
             "F13A1", "FOLR2", "LYVE1") #Int mac


tam <- subset(myeloid, annotation_myeloid %in% c('C1QC+ TAM', 'Interstitial macrophages', 'LA TAM'))
tam <- SetIdent(tam, value="annotation_myeloid")
tam$annotation_myeloid <- factor(tam$annotation_myeloid, levels = c('C1QC+ TAM', 'Interstitial macrophages', 'LA TAM'))

cm_to_inch <- function(cm) cm / 2.54
cairo_pdf("annotation/myeloid/20250819_escc_TAM_markers_dotplot_km.pdf", width = cm_to_inch(23), height = cm_to_inch(12), family = "Arial")
DotPlot(object = tam, features = markers) + 
  theme_bw(base_family = "Arial") +
  theme(line = element_line(linewidth = 0.3),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 45, size = 20, vjust = 0.9, hjust = 0.95),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 9), legend.position = "top", 
        legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.1, "cm"), legend.spacing.x = unit(0.02, 'cm'),
        legend.spacing.y = unit(0.1, "cm"),
        #legend.ticks = element_blank(),
        legend.title = element_text(size= 20))+
  scale_size_area(max_size = 7) +
  ylab("") +scale_color_gradientn(colours=c("blue","white", "red"), name = 'Average\nExpression')
dev.off()


# Fig 5H ----
library(fgsea)
library(tibble)

mygsea_res <- list()

# LA TAM
celltype <- "LA TAM"
re <- readRDS("nebula/myeloid/nebula_LAtam_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoLA TAM`, `se_nebula_annoLA TAM`, `p_nebula_annoLA TAM`)
res$Padj <- p.adjust(res$`p_nebula_annoLA TAM`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "Padj")

df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)

df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("~/Dropbox/database/c5.go.bp.v2025.1.Hs.symbols.gmt")
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

# C1QC+ TAM
celltype <- "C1QC+ TAM"
re <- readRDS("nebula/myeloid/nebula_C1QCposTAM_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoC1QC+ TAM`, `se_nebula_annoC1QC+ TAM`, `p_nebula_annoC1QC+ TAM`)
res$Padj <- p.adjust(res$`p_nebula_annoC1QC+ TAM`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "Padj")

df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)

df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("~/Dropbox/database/c5.go.bp.v2025.1.Hs.symbols.gmt")
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


# Interstitial macrophages
celltype <- "Interstitial macrophages"
re <- readRDS("nebula/myeloid/nebula_Intmac_km.rds")
res <- re$summary %>% dplyr::select(gene, `logFC_nebula_annoInterstitial macrophages`, `se_nebula_annoInterstitial macrophages`, `p_nebula_annoInterstitial macrophages`)
res$Padj <- p.adjust(res$`p_nebula_annoInterstitial macrophages`, method = "bonferroni")
colnames(res) <- c("gene", "avg_log2FC", "se", "p_val", "Padj")

df2 <- res %>% dplyr::mutate(statistics_corrected = avg_log2FC/se)

df3 <- df2 %>% dplyr::select(gene, statistics_corrected)
set.seed(1234); ranks <- deframe(df3)
pw <- gmtPathways("~/Dropbox/database/c5.go.bp.v2025.1.Hs.symbols.gmt")
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

mygsea_res

saveRDS(mygsea_res, "nebula/myeloid/20251005_escc_tam_gsea_gobp_km.rds")
mygsea_res <- readRDS("nebula/myeloid/20251005_escc_tam_gsea_gobp_km.rds")
write.xlsx(mygsea_res, 'nebula/myeloid/20251005_escc_tam_gsea_gobp_km.xlsx')


library(stringr)
mygsea_res <- mygsea_res %>%
  mutate(pathway = pathway %>%
           # Remove 'GOPB_' at the start
           str_remove("^GOBP_") %>%
           # Replace '_' with space
           str_replace_all("_", " ") %>%
           # Convert to lowercase
           str_to_lower() %>%
           # Capitalize first letter of each row (string)
           str_replace("^.", toupper))

mygsea_res_latam_sig <- mygsea_res %>% filter(celltype == 'LA TAM', padj < 0.05, NES > 0)
mygsea_res_c1qctam_sig <- mygsea_res %>% filter(celltype == 'C1QC+ TAM', padj < 0.05, NES > 0)
mygsea_res_intmacro_sig <- mygsea_res %>% filter(celltype == 'Interstitial macrophages', padj < 0.05, NES > 0)

#dotplot
latam_sig_term <- c('Lipid catabolic process',
                    'Glycolipid catabolic process',
                    'Phospholipid catabolic process',
                    'Membrane lipid catabolic process')
c1qctam_gobp_sig_term <- c('Cytoplasmic translation',
                           'Ribosomal small subunit biogenesis',
                           'Ribosome biogenesis',
                           'Ribonucleoprotein complex biogenesis',
                           'Peptide antigen assembly with mhc protein complex',
                           'Protein rna complex organization')
intmac_gobp_sig_term <- c('Regulation of alternative mrna splicing via spliceosome',
                          'Myeloid cell differentiation',
                          'Outflow tract morphogenesis',
                          'Intermediate filament organization',
                          'Regulation of cell migration involved in sprouting angiogenesis')

all_gobp_sig_term <- unique(c(latam_sig_term,c1qctam_gobp_sig_term,intmac_gobp_sig_term))

mygsea_res_plot <- mygsea_res %>% filter(pathway %in% all_gobp_sig_term, padj < 0.05)
mygsea_res_plot$celltype <- factor(mygsea_res_plot$celltype, levels = c('LA TAM', 'C1QC+ TAM','Interstitial macrophages'))
mygsea_res_plot$pathway <- factor(mygsea_res_plot$pathway, levels = unique(c(mygsea_res_plot$pathway[mygsea_res_plot$celltype == 'LA TAM'],
                                                                             mygsea_res_plot$pathway[mygsea_res_plot$celltype == 'C1QC+ TAM'],
                                                                             mygsea_res_plot$pathway[mygsea_res_plot$celltype == 'Interstitial macrophages'])))
mygsea_res_plot$pathway <- factor(mygsea_res_plot$pathway, levels = unique(c('Membrane lipid catabolic process',
                                                                             'Phospholipid catabolic process',
                                                                             'Lipid catabolic process',
                                                                             'Glycolipid catabolic process',
                                                                             'Peptide antigen assembly with mhc protein complex',
                                                                             'Cytoplasmic translation',
                                                                             'Ribosomal small subunit biogenesis',
                                                                             'Ribosome biogenesis',
                                                                             'Ribonucleoprotein complex biogenesis',
                                                                             'Protein rna complex organization',
                                                                             'Outflow tract morphogenesis',
                                                                             'Regulation of cell migration involved in sprouting angiogenesis',
                                                                             'Regulation of alternative mrna splicing via spliceosome',
                                                                             'Intermediate filament organization',
                                                                             'Myeloid cell differentiation')))

library(ggplot2)
library(showtext)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)
cm_to_inch <- function(cm) cm / 2.54

cairo_pdf("nebula/myeloid/20251006_escc_sc.snrna_tam_nebula_GSEA_gobp_sig_km.pdf", width = cm_to_inch(21), height = cm_to_inch(17), family = "Arial")
ggplot(mygsea_res_plot, aes(celltype,pathway)) + geom_point(aes(color = NES, size = -log10(padj))) + 
  ylab(NULL) + theme_minimal(base_family = 'Arial') +
  theme(axis.text.x = element_text(colour = "black", angle = 45, size = 15, vjust=1, hjust=1), 
        axis.text.y = element_text(colour = "black", size = 15), plot.title = element_blank(), 
        legend.position = "right", legend.text = element_text(size = 15),
        legend.title = element_text(size = 15), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black", size=0.5),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.spacing = unit(1, 'cm'), legend.spacing.x = unit(0.5,'cm')) + 
  guides(size = guide_legend(title = expression(-log[10] ~ "(padj)"))) + 
  #scale_size(range = c(2,6), breaks = c(10,20,30)) +
  scale_color_gradientn(colours=c("blue","pink","red")) +
  scale_shape_manual(values = c(13,16))
dev.off()


# Fig 5I ----
cellchat_normal <- readRDS(file = "cellchat/cd4t/R_objects/20250806_escc_normal_hj.rds")
cellchat_dysp <- readRDS(file = "cellchat/cd4t/R_objects/20250806_escc_dysplasia_hj.rds")
cellchat_micro <- readRDS(file = "cellchat/cd4t/R_objects/20250806_escc_microinvasive_carcinoma_hj.rds")
cellchat_macro <- readRDS(file = "cellchat/cd4t/R_objects/20250806_escc_macroinvasive_carcinoma_hj.rds")
object.list <- list(Normal = cellchat_normal,
                    Dysplasia = cellchat_dysp,
                    Micro_invasive = cellchat_micro,
                    Macro_invasive = cellchat_macro)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

library(showtext)
font_add("Arial", regular = "/Library/Fonts/Arial Unicode.ttf")
showtext_auto(TRUE)

cairo_pdf("~/Dropbox/escc_figure_share/Fig5/20250822/20250822_escc_Fig5_Treg_CTL_cellchat_km.pdf", width = 9, height = 9, family="Arial")
netVisual_bubble(cellchat, sources.use = c('C1QC+ TAM','LA TAM','Interstitial macrophages'), 
                 targets.use = c('Tumor'), comparison = c(1,2,3,4)) +
  theme_classic(base_family = 'Arial') +
  theme(axis.text.x = element_text(size = 18, colour = c('#83f52c', '#f3f315', '#ff6600', '#ff0099'), angle = 60, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18, colour = 'black'),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))
dev.off()

