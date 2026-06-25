library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(biomaRt)
library(tidyverse)


Convert("/Users/tobymilne-clark/Documents/Data/Atlas/scRNA_prostate/H5ad/External/5739e080-4dbc-48e2-91dd-feaa1d7ad773.h5ad", dest = "h5seurat", overwrite = TRUE)

source("/Users/tobymilne-clark/Documents/Data/Atlas/scRNA_prostate/H5ad/LoadH5Seurat_v5.R")


obj <- LoadH5Seurat("/Users/tobymilne-clark/Documents/Data/Atlas/scRNA_prostate/H5ad/External/5739e080-4dbc-48e2-91dd-feaa1d7ad773.h5seurat")

saveRDS(obj , file = "/Users/tobymilne-clark/Documents/Data/Atlas/scRNA_prostate/H5ad/External/Apostolov_et_al_2024_h5ad_Seurat_converted_TMC_160426.rds")

neuronal_genes <- c("NRXN1", "NRXN3", "NTRK3", "SLITRK3", "SLITRK4", "GRIN2D", "GRIA3", "SHANK2", 
                    "CAMK2A", "CAMK2B", "CAMKK2", "GABRA5", "GABRB1", "GABRG3", "SLC32A1", "SYT7", 
                    "ADCY1", "RASGRF1", "CHAT", "KIF17", "PPFIBP2")
NEPC_genes <- c("INSM1", "CHGA", "CHGB", "SYP")
rho_genes <- c("ARHGAP25", "ARHGAP26", "ARHGEF7", "BAIAP2", "BCR", "EPHA2","FAM13A", "FARP1", 
               "ITGB1", "MCF2", "PAK6", "PLEKHG3","SRGAP3", "TAOK3", "VAV2", "VAV3")


mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)

neuronal_genes_ENS <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = neuronal_genes,
  mart = mart
)

NEPC_genes_ENS <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = NEPC_genes,
  mart = mart
)

rho_genes_ENS <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = rho_genes,
  mart = mart
)
getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = c("NTF4","NGFR","NTRK2"),
  mart = mart
)

obj <- AddModuleScore(obj, features = list(neuronal_genes_ENS$ensembl_gene_id, NEPC_genes_ENS$ensembl_gene_id, rho_genes_ENS$ensembl_gene_id,"ENSG00000122585","ENSG00000114812","ENSG00000225950"),
                      name=c("neuronal", "NEPC", "rho","NPY","VIPR1","NTF4"))


levels(obj@meta.data$cell_type)

cell_color_code_extended <- c(
  "basal cell of prostate epithelium"           = "#336699",  
  "epithelial cell of prostate"                 = "#66B8D6", 
  "luminal cell of prostate epithelium"         = "#9DCEE2",  
  "club-like cell of the urethral epithelium"   = "#3399CC",  
  "ciliated epithelial cell"                    = "#333366",  
  "prostate neuroendocrine cell"                = "#1F4E79",
  
  "fibroblast of connective tissue of prostate" = "#5CAE63",
  "smooth muscle cell of prostate"              = "#F8A627",
  "vascular associated smooth muscle cell"      = "#F8A627",
  "pericyte"                                    = "#E78F1E",
  
  "endothelial cell"                            = "#986EAC",
  "endothelial cell of lymphatic vessel"        = "#B497C9",
  "vein endothelial cell"                       = "#7A4E99",
  
  "T cell"                                      = "#FEDFD4",
  "regulatory T cell"                           = "#FDD0C2",
  "cytotoxic T cell"                            = "#FCAE91",
  "natural killer cell"                         = "#F8A592",
  
  "B cell"                                      = "#65010C",
  "plasma cell"                                 = "#8B1E2D",
  
  "monocyte"                                    = "#CB1B16",
  "myeloid cell"                                = "#D7301F",
  "macrophage"                                  = "#E34A33",
  "plasmacytoid dendritic cell"                 = "#FB6A4A",
  
  "mast cell"                                   = "#F26A4F",
  
  "glial cell"                                  = "#B8A3C9",
  "myelinating Schwann cell"                    = "#7F5A99",
  
  "unknown"                                     = "#BDBDBD"
)

DimPlot(obj, reduction = "umap", group.by = "cell_type", cols = cell_color_code_extended, label = TRUE) 
FeaturePlot(obj, features = "ENSG00000122585") + DimPlot(obj, reduction = "umap", group.by = "cell_type", cols = cell_color_code_extended, label = TRUE)
FeaturePlot(obj, features = "ENSG00000146469") + DimPlot(obj, reduction = "umap", group.by = "cell_type", cols = cell_color_code_extended, label = TRUE)
FeaturePlot(obj, features = "ENSG00000064300") + DimPlot(obj, reduction = "umap", group.by = "cell_type", cols = cell_color_code_extended, label = TRUE)
FeaturePlot(obj, features = "ENSG00000225950") + DimPlot(obj, reduction = "umap", group.by = "cell_type", cols = cell_color_code_extended, label = TRUE)
FeaturePlot(obj, features = "ENSG00000148053") + DimPlot(obj, reduction = "umap", group.by = "cell_type", cols = cell_color_code_extended, label = TRUE)

obj_df <- obj@meta.data 
colnames(obj_df)

cell_order <- obj_df %>%
  group_by(cell_type) %>%
  summarise(stat_value = median(neuronal1, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(stat_value)) %>%
  pull(cell_type)

highlight_df <- obj_df %>%
  group_by(cell_type) %>%
  summarise(stat_value = median(neuronal1, na.rm = TRUE), .groups = "drop") %>%
  mutate(enriched = stat_value > 0)

obj_df_plot <- obj_df %>%
  left_join(highlight_df, by = "cell_type") %>%
  mutate(cell_type = factor(cell_type, levels = cell_order))

ggplot(
  obj_df_plot,
  aes(x = cell_type, y = neuronal1, fill = enriched)
) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_boxplot(width = 0.7, outlier.size = 0.2, alpha = 0.85) +
  scale_fill_manual(values = c("TRUE" = "#D95F02", "FALSE" = "grey80")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Neuronal module score")

ggplot(
  obj_df,
  aes(x = fct_reorder(celltype_major_v2, neuronal1, .fun = median, .desc = TRUE, na.rm = TRUE),
      y = neuronal1,fill = ISUP_clinical)
) +
  geom_boxplot(width = 0.7, outlier.size = 0.2, alpha = 0.85) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = "Neuronal module score")

ggplot(obj_df, aes(x = cell_type, y = NPY4, fill = ISUP_clinical)) +
  geom_boxplot(width = 0.7, outlier.size = 0.2, alpha = 0.85) +
  scale_fill_manual(values = c("#4C78A8", "#F58518", "#54A24B","purple")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = "NPY module score")

ggplot(obj_df, aes(x = cell_type, y = VIPR15, fill = ISUP_clinical)) +
  geom_boxplot(width = 0.7, outlier.size = 0.2, alpha = 0.85) +
  scale_fill_manual(values = c("#4C78A8", "#F58518", "#54A24B","purple")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = "VIPR1 module score")

ggplot(obj_df, aes(x = cell_type, y = DPP46, fill = ISUP_clinical)) +
  geom_boxplot(width = 0.7, outlier.size = 0.2, alpha = 0.85) +
  scale_fill_manual(values = c("#4C78A8", "#F58518", "#54A24B","purple")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = "DPP4 module score")





curated_genes <- list(
  Inflammation = c("IFNG","STAT1","IRF1","CXCL9","CXCL10","TNF","NFKBIA","IL1B"),
  Cytotoxicity = c("GZMB","PRF1","NKG7","GNLY","CD8A","CD8B"),
  Exhaustion = c("PDCD1","CTLA4","LAG3","TIGIT","HAVCR2")
)

Inflammation_genes_ENS <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = c("IFNG","STAT1","IRF1","CXCL9","CXCL10","TNF","NFKBIA","IL1B"),
  mart = mart
)

Cytotoxicity_genes_ENS <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = c("GZMB","PRF1","NKG7","GNLY","CD8A","CD8B"),
  mart = mart
)

Exhaustion_genes_ENS <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = c("PDCD1","CTLA4","LAG3","TIGIT","HAVCR2"),
  mart = mart
)


obj <- AddModuleScore(obj, features = list(Inflammation_genes_ENS,Cytotoxicity_genes_ENS,Exhaustion_genes_ENS),
                      name=c("Inflammation", "Cytotoxicity", "Exhaustion"))


ggplot(obj_df, aes(x = cell_type, y = Inflammation1, fill = ISUP_clinical)) +
  geom_boxplot(width = 0.7, outlier.size = 0.2, alpha = 0.85) +
  scale_fill_manual(values = c("#4C78A8", "#F58518", "#54A24B","purple")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = "DPP4 module score")

ggplot(obj_df, aes(x = cell_type, y = Cytotoxicity2, fill = ISUP_clinical)) +
  geom_boxplot(width = 0.7, outlier.size = 0.2, alpha = 0.85) +
  scale_fill_manual(values = c("#4C78A8", "#F58518", "#54A24B","purple")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = "DPP4 module score")

ggplot(obj_df, aes(x = cell_type, y = Exhaustion3, fill = ISUP_clinical)) +
  geom_boxplot(width = 0.7, outlier.size = 0.2, alpha = 0.85) +
  scale_fill_manual(values = c("#4C78A8", "#F58518", "#54A24B","purple")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = "DPP4 module score")







