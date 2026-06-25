####################
##### Giotto #######
####################
library(tidyverse)
library(Seurat)
library(Giotto)
# if data.table corrupt run .rs.restartR()

# ============================================================
# Define samples to loop through
# ============================================================
sample_ids <- paste0("NEADT", 9:11)

# Base paths - update these as needed
base_data_path    <- "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics"
base_results_path <- "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Processed_data"

# ============================================================
# load scRNA-seq reference data
# ============================================================


ref <- readRDS("/Users/tobymilne-clark/Library/Mobile Documents/com~apple~CloudDocs/PhD/Data/Prostate_scRNA/scRNA/scRNAseq.exampleData.20220606.RDS")
ref <- UpdateSeuratObject(ref)

markers <- FindAllMarkers(ref,
                          group.by       = "celltype",
                          logfc.threshold = 0.5)

ref_matrix <- Seurat::GetAssayData(ref, assay = "RNA", layer = "data")


# ============================================================
# Make sign matrices for PAGE enrichment
# ============================================================

neuronal_genes <- c("NRXN1", "NRXN3", "NTRK3", "SLITRK3", "SLITRK4", "GRIN2D", "GRIA3", "SHANK2", "CAMK2A", "CAMK2B", "CAMKK2", "GABRA5", "GABRB1", "GABRG3", "SLC32A1", "SYT7", "ADCY1", "RASGRF1", "CHAT", "KIF17", "PPFIBP2")
Rho_GTPase_genes <- c(  "ARHGAP25", "ARHGAP26", "ARHGEF7", "BAIAP2", "BCR", "EPHA2",   "FAM13A", "FARP1", "ITGB1", "MCF2", "PAK6", "PLEKHG3",   "SRGAP3", "TAOK3", "VAV2", "VAV3")

expressed_DMRs <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Metadata/DMRs_annotated_expressed.csv")

expressed_DMRs <- expressed_DMRs %>% mutate(direction = if_else(difference > 0, "Hyper", "Hypo")) 

genes_list3 <- expressed_DMRs %>%
  group_by(direction) %>%
  summarise(genes = list(annot.symbol), .groups = "drop") %>%
  deframe()

Hallmark_AR = m_t2g_2 %>% filter(gs_name == "HALLMARK_ANDROGEN_RESPONSE") %>% pull(ensembl_gene)
Prolaris = m_t2g_2 %>% filter(gs_name == "Prolaris") %>% pull(ensembl_gene)


Meth_genes <- list(neuronal = neuronal_genes,
                   Rho_GTPase = Rho_GTPase_genes,
                   Hallmark_AR = Table_Genes$SYMBOL[match(Hallmark_AR,Table_Genes$ENSEMBL)],
                   Prolaris = Table_Genes$SYMBOL[match(Prolaris,Table_Genes$ENSEMBL)]) 

Meth_genes <- c(Meth_genes ,genes_list3)

sign_matrix_meth <- makeSignMatrixPAGE(sign_names = names(Meth_genes),
                                       sign_list = Meth_genes )
cell_types_Meth <- colnames(sign_matrix_meth)



basal_kp <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/Results_and_scripts/PPInetworks_Basal/data/key_proteins.txt") 
luminal_kp <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/Results_and_scripts/PPInetworks_Luminal/data/key_proteins.txt")
Immune_kp <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/Results_and_scripts/PPInetworks_Immune/data/key_proteins.txt")
OE_kp <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/Results_and_scripts/PPInetworks_OE/data/key_proteins.txt")
fibroblast_kp <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/Results_and_scripts/PPInetworks_Fibroblast/data/key_proteins.txt")
Endothelial_kp <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/Results_and_scripts/PPInetworks_Endothelial/data/key_proteins.txt")

key_proteins_cells <- list(basal_kp = basal_kp$Gene,
                           luminal_kp = luminal_kp$Gene,
                           Immune_kp = Immune_kp$Gene,
                           OE_kp = OE_kp$Gene,
                           fibroblast_kp = fibroblast_kp$Gene,
                           Endothelial_kp = Endothelial_kp$Gene)


sign_matrix_KP <- makeSignMatrixPAGE(sign_names = names(key_proteins_cells),
                                     sign_list = key_proteins_cells )

cell_types_KP <- colnames(sign_matrix_KP)


genes_basal <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/cell_atlas/Genes/Genes_Basal_cells.csv", col_names = FALSE) 
genes_luminal <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/cell_atlas/Genes/Genes_Luminal_cells.csv", col_names = FALSE)
genes_immune <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/cell_atlas/Genes/Genes_Immune_cells.csv", col_names = FALSE)
genes_OE <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/cell_atlas/Genes/Genes_Other_Epi_cells.csv", col_names = FALSE)
genes_fibroblast <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/cell_atlas/Genes/Genes_Fibroblast_cells.csv", col_names = FALSE)
genes_Endothelial <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/cell_atlas/Genes/Genes_Endothelial_cells.csv", col_names = FALSE)

genes_cell_atlas <- list(basal = genes_basal$X1,
                         luminal = genes_luminal$X1,
                         immune = genes_immune$X1,
                         OE = genes_OE$X1,
                         fibroblast = genes_fibroblast$X1,
                         Endothelial = genes_Endothelial$X1)



sign_matrix_Atlas<- makeSignMatrixPAGE(sign_names = names(genes_cell_atlas),
                                       sign_list = genes_cell_atlas )
cell_types_Atlas <- colnames(sign_matrix_Atlas)


# ============================================================
# Main loop
# ============================================================
for (sample_id in sample_ids) {
  
  message("\n========================================")
  message("Processing sample: ", sample_id)
  message("========================================\n")
  
  # ----------------------------------------------------------
  # 1. Set up paths and output folder
  # ----------------------------------------------------------
  data_path      <- file.path(base_data_path, sample_id)
  results_folder <- file.path(base_results_path, paste0("results_Giotto_", sample_id))
  
  # Create output directory if it doesn't exist
  if (!dir.exists(results_folder)) {
    dir.create(results_folder, recursive = TRUE)
    message("Created output folder: ", results_folder)
  }
  
  instructions <- createGiottoInstructions(
    save_dir   = results_folder,
    save_plot  = TRUE,
    show_plot  = FALSE,
    return_plot = FALSE,
    plot_format = "pdf",
    dpi        = 600
  )
  
  # ----------------------------------------------------------
  # 2. Load Visium data
  # ----------------------------------------------------------
  visium_obj <- createGiottoVisiumObject(
    visium_dir       = data_path,
    expr_data        = "filter",
    png_name         = "tissue_lowres_image.png",
    gene_column_index = 2,
    instructions     = instructions
  )
  
  ## Check metadata
  pDataDT(visium_obj)
  
  ## Show available image names
  showGiottoImageNames(visium_obj)
  
  ## Show aligned image
#  spatPlot(gobject    = visium_obj,
#           cell_color = "in_tissue",
#           show_image = TRUE,
#           point_alpha = 0.7)
  
  # ----------------------------------------------------------
  # 3. Subset to tissue-covered spots
  # ----------------------------------------------------------
  metadata          <- pDataDT(visium_obj)
  in_tissue_barcodes <- metadata[in_tissue == 1]$cell_ID
  
  visium_obj <- subsetGiotto(visium_obj,
                             cell_ids = in_tissue_barcodes)
  
  # ----------------------------------------------------------
  # 4. Filter
  # ----------------------------------------------------------
  visium_obj <- filterGiotto(
    gobject               = visium_obj,
    expression_threshold  = 1,
    feat_det_in_min_cells = 50,
    min_det_feats_per_cell = 1000,
    expression_values     = "raw",
    verbose               = TRUE
  )
  
  # ----------------------------------------------------------
  # 5. Normalize
  # ----------------------------------------------------------
  visium_obj <- normalizeGiotto(gobject     = visium_obj,
                                scalefactor = 6000,
                                verbose     = TRUE)
  
  ## Add gene & cell statistics
  visium_obj <- addStatistics(gobject = visium_obj)
  
  ## Visualize
#  spatPlot2D(gobject    = visium_obj,
#             show_image = TRUE,
#             point_alpha = 0.7)
  
#  spatPlot2D(gobject        = visium_obj,
#             show_image     = TRUE,
#             point_alpha    = 0.7,
#             cell_color     = "nr_feats",
#             color_as_factor = FALSE)
  
  # ----------------------------------------------------------
  # 6. Dimensionality reduction
  # ----------------------------------------------------------
  visium_obj <- calculateHVF(gobject   = visium_obj,
                             save_plot = TRUE)
  
  visium_obj <- runPCA(gobject = visium_obj)
  
  screePlot(visium_obj, ncp = 30)
  
  plotPCA(gobject           = visium_obj,
         dim1_to_use = 1, dim2_to_use = 2)
  
  visium_obj <- runUMAP(gobject = visium_obj, 
                        dimensions_to_use = 1:10)
  
  plotUMAP(gobject = visium_obj)
  
  visium_obj <- runtSNE(visium_obj,
                        dimensions_to_use = 1:10)
  
  plotTSNE(gobject = visium_obj)
  
  # ----------------------------------------------------------
  # 7. Clustering
  # ----------------------------------------------------------
  ## sNN network (default)
  visium_obj <- createNearestNetwork(gobject           = visium_obj,
                                     dimensions_to_use = 1:10,
                                     k                 = 15)
  
  ## Leiden clustering
  visium_obj <- doLeidenCluster(gobject      = visium_obj,
                                resolution   = 0.2,
                                n_iterations = 1000)
  
  plotUMAP(gobject        = visium_obj,
           cell_color     = "leiden_clus",
           show_NN_network = TRUE,
           point_size     = 2.5)
  
  # ----------------------------------------------------------
  # 8. Marker genes from reference
  # ----------------------------------------------------------

  markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::filter(gene %in% visium_obj@feat_ID$rna) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup() -> top10
  

  # ----------------------------------------------------------
  # 9. DWLS deconvolution
  # ----------------------------------------------------------
  sign_matrix <- makeSignMatrixDWLSfromMatrix(
    matrix    = ref_matrix,
    cell_type = droplevels(ref$celltype),
    sign_gene = top10$gene
  )
  
  visium_obj <- runDWLSDeconv(visium_obj,
                              sign_matrix = sign_matrix)
  
  spatDeconvPlot(visium_obj,
                 show_image = FALSE,
                 radius     = 100,
                 cell_color_code = c(
                   "BE" = "#333366",
                   "OE1" = "#336699",
                   "OE2" = "#3399CC",
                   "LE-KLK3" = "#D7EDF4",
                   "LE-KLK4" = "#9DCEE2",
                   "Fibroblast" = "#5CAE63",
                   "Endothelial" = "#986EAC",
                   "Tcell" = "#FEDFD4",
                   "Bcell" = "#65010C",
                   "MNP" = "#CB1B16",
                   "Mastcell" = "#F26A4F",
                   "NKcell" = "#F8A592",
                   "Sperm" = "#DEDEDE"
                 ))

    # ----------------------------------------------------------
  # 10. PAGE enrichment - Methylation
  # ----------------------------------------------------------
  visium_obj <- runPAGEEnrich(gobject          = visium_obj,
                              sign_matrix      = sign_matrix_meth,
                              min_overlap_genes = 1,
                              name             = "Meth_enrich")
  
  cell_types_Meth <- colnames(sign_matrix_meth)
  
  spatCellPlot2D(gobject                = visium_obj,
                 spat_enr_names         = "Meth_enrich",
                 cell_annotation_values = cell_types_Meth[1:6],
                 cow_n_col              = 2,
                 coord_fix_ratio        = 1,
                 point_size             = 1,
                 gradient_midpoint      = 0,
                 gradient_style         = "divergent",
                 cell_color_gradient    = c("darkblue", "blue","lightblue", "white", "yellow", "orange", "darkred"),
                 show_legend            = TRUE)
  
  # ----------------------------------------------------------
  # 11. PAGE enrichment - Atlas
  # ----------------------------------------------------------
  visium_obj <- runPAGEEnrich(gobject          = visium_obj,
                              sign_matrix      = sign_matrix_Atlas,
                              min_overlap_genes = 1,
                              name             = "Atlas_enrich")
  
  cell_types_Atlas <- colnames(sign_matrix_Atlas)
  
  spatCellPlot2D(gobject                = visium_obj,
                 spat_enr_names         = "Atlas_enrich",
                 cell_annotation_values = cell_types_Atlas[1:6],
                 cow_n_col              = 2,
                 coord_fix_ratio        = 1,
                 point_size             = 1,
                 gradient_midpoint      = 0,
                 gradient_style         = "divergent",
                 cell_color_gradient    = c("darkblue", "blue","lightblue", "white", "yellow", "orange", "darkred"),
                 show_legend            = TRUE)
  
  # ----------------------------------------------------------
  # 12. PAGE enrichment - Key Proteins
  # ----------------------------------------------------------
  visium_obj <- runPAGEEnrich(gobject          = visium_obj,
                              sign_matrix      = sign_matrix_KP,
                              min_overlap_genes = 1,
                              name             = "Key_protein_enrich")
  
  cell_types_KP <- colnames(sign_matrix_KP)
  
  spatCellPlot2D(gobject                = visium_obj,
                 spat_enr_names         = "Key_protein_enrich",
                 cell_annotation_values = cell_types_KP[1:6],
                 cow_n_col              = 2,
                 coord_fix_ratio        = 1,
                 point_size             = 1,
                 gradient_midpoint      = 0,
                 gradient_style         = "divergent",
                 cell_color_gradient    = c("darkblue", "blue","lightblue", "white", "yellow", "orange", "darkred"),
                 show_legend            = TRUE)

    # ----------------------------------------------------------
  # 13. Optionally save the Giotto object for this sample
  # ----------------------------------------------------------
  saveRDS(visium_obj,
          file = file.path(results_folder, paste0(sample_id, "_giotto_object.rds")))
  
  message("Finished sample: ", sample_id, "\n")
  
} # end for loop

message("All samples processed!")