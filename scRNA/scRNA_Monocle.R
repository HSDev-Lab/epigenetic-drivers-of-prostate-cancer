library("Seurat")
library("seurat-disk")
library(data.table)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(GenomicRanges)
library('scITD')


##########################
### scRNA Reference 
##########################


library(Seurat)
ref <- readRDS("/Users/tobymilne-clark/Library/Mobile Documents/com~apple~CloudDocs/PhD/Data/Prostate_scRNA/scRNA/scRNAseq.exampleData.20220606.RDS")
ref <- UpdateSeuratObject(ref)

DimPlot(ref, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() +
  ggtitle("scRNA Reference") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("/Users/tobymilne-clark/Documents/PhD_Thesis/Plots/RNA-seq/scRNA_reference_UMAP.svg", dpi = 300, width = 3.5,
       height = 3.5)


markers <- FindAllMarkers(ref, group.by = "celltype", logfc.threshold = 0.5)
markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup() -> top10

DoHeatmap(ref,  features = top10$gene,size = 2) + NoLegend() + 
  theme(
    plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title    = element_text(size = 10),
    axis.text     = element_text(size = 4),
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 6),
    text = element_text(size = 6),
    legend.position = "none"
  )


ggsave("/Users/tobymilne-clark/Documents/PhD_Thesis/Plots/RNA-seq/scRNA_reference_top10_heatmap.svg", dpi = 300, width = 7,
       height = 4)


DotPlot(ref, features = c("EPCAM","KRT5","KLK3","KLK4")) + RotatedAxis()
??Dotplot

DotPlot(ref, split.by = "group" ,features = c("CHGA","CHGB")) + RotatedAxis() 
FeaturePlot(ref, split.by = "group" ,features = c("CHGA","CHGB")) + RotatedAxis() 

DotPlot(ref, split.by = "group" ,features = c("ANPEP", "ABL1", "PSCA", "EFNA1", "HSPB1", "INMT",  "TRIP13"
)) + RotatedAxis() 

"ANPEP", "ABL1", "PSCA", "EFNA1", "HSPB1", "INMT",  "TRIP13"

library(CytoTRACE2)


cyto_data <- CytoTRACE2::loadData_fromSeurat(ref)

# running CytoTRACE 2 main function - cytotrace2 - with default parameters
cytotrace2_result <- CytoTRACE2::cytotrace2(ref,is_seurat = TRUE, species = "human")


cytotrace2_result@plots$CytoTRACE2_UMAP

annotation <- data.frame(phenotype = ref@meta.data$celltype)

rownames(annotation) <- colnames(ref)

plots <- CytoTRACE2::plotData(cytotrace2_result,
                              annotation = annotation, 
                              is_seurat = TRUE)


plots$CytoTRACE2_UMAP

plots$CytoTRACE2_Potency_UMAP

plots$CytoTRACE2_Relative_UMAP

plots$CytoTRACE2_Boxplot_byPheno




library(monocle3)

cds <- as.cell_data_set(ref)
cds <- cluster_cells(cds, resolution=1e-3)

cds <- learn_graph(cds)

## Step 6: Order cells
cds <- order_cells(cds)

rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

plot_cells(cds,color_cells_by = "celltype")

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = TRUE)


ggsave("/Users/tobymilne-clark/Documents/PhD_Thesis/Plots/RNA-seq/celltype_pseudotime.pdf", dpi = 300, width = 4.5,
       height = 3.5)


cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)
cds_graph_test_results %>% 
  arrange(morans_I) %>%
  dplyr::filter(q_value < 0.05) %>%
  dplyr::select(morans_I, q_value) %>%
  head(10)

pr_deg_ids <- row.names(subset(cds_graph_test_results, q_value < 0.05))


deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

deg_cds <- preprocess_cds(cds[deg_ids,])
gene_modules <- find_gene_modules(deg_cds,
                                  resolution=c(10^seq(-6,-1)))
table(gene_modules$module)

cell_groups <- data.frame(cell = row.names(colData(deg_cds)),
                          cell_group = colData(deg_cds)$celltype)
agg_mat <- aggregate_gene_expression(deg_cds,
                                     gene_group_df = gene_modules,
                                     cell_group_df = cell_groups)

row.names(agg_mat) <- paste0("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column",
                   treeheight_row = 0,
                   treeheight_col = 0,
                   clustering_method="ward.D2")

plot_cells(deg_cds,
           genes=gene_modules %>% filter(module %in% c(3,9,6,2)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


plot_cells(deg_cds,
           genes=gene_modules %>% filter(module %in% c(18,13,42,34)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

rowData(deg_cds)$gene_name <- rownames(deg_cds)
rowData(deg_cds)$gene_short_name <- rowData(deg_cds)$gene_name

gm <- gene_modules[which(gene_modules$module %in% c(6, 2, 12,3)),]
plot_cells(deg_cds,
           genes=gm,
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           trajectory_graph_color = "grey60")

neuronal_genes <- c("NRXN1", "NRXN3", "NTRK3", "SLITRK3", "SLITRK4", "GRIN2D", "GRIA3", "SHANK2", "CAMK2A", "CAMK2B", "CAMKK2", "GABRA5", "GABRB1", "GABRG3", "SLC32A1", "SYT7", "ADCY1", "RASGRF1", "CHAT", "KIF17", "PPFIBP2")
AFD_genes <- c("NRXN1", "NRXN3", "NTRK3")

AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$celltype %in% c("LE-KLK3", "LE-KLK4", "BE","OE1","OE2")]

AFD_lineage_cds <- order_cells(AFD_lineage_cds)


plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="celltype",
                         min_expr=0.5)

?plot_genes_in_pseudotime

cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d)

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")

tol_muted_14 <- c(
  "#88CCEE",  # light blue
  "#44AA99",  # teal
  "#117733",  # dark green
  "#332288",  # purple
  "#DDCC77",  # mustard yellow
  "#999933",  # olive
  "#661100",  # brown
  "#CC6677",  # pink
  "#882255",  # wine
  "#AA4499",  # magenta
  "#DDDDDD",  # light grey
  "#000000",  # black
  "#E69F00",  # orange
  "#56B4E9"   # sky blue
)
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="celltype", color_palette = NULL)
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="celltype",color_palette = tol_muted_14)
?plot_cells_3d
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="pseudotime")
?plot_cells_3d




Epithelial <- ref[,  ref$celltype %in% c("LE-KLK3", "LE-KLK4", "BE","OE1","OE2")]
Epithelial.cds <- as.cell_data_set(Epithelial)
Epithelial.cds <- cluster_cells(cds = Epithelial.cds, reduction_method = "UMAP")
Epithelial.cds <- learn_graph(Epithelial.cds, use_partition = TRUE)
Epithelial.cds <- order_cells(Epithelial.cds, reduction_method = "UMAP")
plot_cells(
  cds = Epithelial.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

tumour <- ref[,  ref$group %in% c("tumor")]
tumour.cds <- as.cell_data_set(tumour)
tumour.cds  <- cluster_cells(cds = tumour.cds , reduction_method = "UMAP")
tumour.cds  <- learn_graph(tumour.cds , use_partition = TRUE)
tumour.cds  <- order_cells(tumour.cds , reduction_method = "UMAP")
plot_cells(
  cds = tumour.cds ,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)



c(18,13,42,34)
mod <- gene_modules %>% filter(module == "3") 


?write.csv
write.csv(mod$id, "/Users/tobymilne-clark/Documents/Data/Atlas/scRNA_prostate/scRNA_modules.csv", row.names = FALSE,quote = FALSE)

markers %>% filter(cluster == "Fibroblast") %>% filter(gene %in% mod$id) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC, p_val_adj,cluster) -> mod_markers

logChanges_scRNA <- mod_markers$avg_log2FC
names(logChanges_scRNA) <- rownames(mod_markers)

ranks <- sort(logChanges_scRNA, decreasing = TRUE)


msigdbr_species()
m_df <- msigdbr(species = "Homo sapiens")

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, db_gene_symbol)
head(m_t2g)

em2 <- GSEA(ranks, TERM2GENE = m_t2g, pvalueCutoff = 0.05)

dotplot(em2, showCategory=10,split=".sign") + 
  facet_grid(.~.sign) +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    base_size = 12, base_family = "Arial",
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title    = element_text(size = 12),
    axis.text     = element_text(size = 10),
    legend.title  = element_text(size = 10),
    legend.text   = element_text(size = 9)
  )


