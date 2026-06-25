library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library("infercnv")
library(SpaCET)

localdir <- "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/VisiumHD"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))

# Setting default assay changes between 8um and 16um binning
Assays(object)
DefaultAssay(object) <- "Spatial.008um"

vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot


# normalize both 8um and 16um bins
DefaultAssay(object) <- "Spatial.008um"
object <- NormalizeData(object)

DefaultAssay(object) <- "Spatial.016um"
object <- NormalizeData(object)

# switch spatial resolution to 16um from 8um
DefaultAssay(object) <- "Spatial.016um"
p1 <- SpatialFeaturePlot(object, features = "KLK3") + ggtitle("KLK3 expression (16um)")

# switch back to 8um
DefaultAssay(object) <- "Spatial.008um"
p2 <- SpatialFeaturePlot(object, features = "KLK3") + ggtitle("KLK3 expression (8um)")

KLK3_resolution <- p1 | p2


# note that data is already normalized
DefaultAssay(object) <- "Spatial.008um"
object <- FindVariableFeatures(object)
object <- ScaleData(object)
# we select 50,0000 cells and create a new 'sketch' assay
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch",
  features = VariableFeatures(object)
)


# switch analysis to sketched cells
DefaultAssay(object) <- "sketch"

# perform clustering workflow
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
#object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 3)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 2)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)")

options(future.globals.maxSize = 8000 * 1024^2)
object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + NoLegend()

# switch to full dataset
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label = F) + ggtitle("Projected clustering (full dataset)") + NoLegend()

p1 | p2

SpatialDimPlot(object, label = T, repel = T, label.size = 4)





object <- add_gene_signatures(object,genes_list = genes_list)


FeaturePlot(object, features = c("tumor","luminal epithelium","basal epithelium","sensescent epithelium")) +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

FeaturePlot(object, features = c("Neuronal","fibroblast","muscle","capillary venous endothelium")) +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "tumor") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "luminal epithelium") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "basal epithelium") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "Macrophages") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "Monocytes") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "T cells") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "Mast cells") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")



SpatialFeaturePlot(object, features = "Neuronal") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

object <- add_gene_signatures(object,genes_list = genes_list2)

SpatialFeaturePlot(object, features = "Luminal") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "Immune") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "Fibroblast") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "Endothelial") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "Other_Epi") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")


SpatialFeaturePlot(object, features = c("Fibroblast","Endothelial")) +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

object <- add_gene_signatures(object,genes_list = genes_list3)


SpatialFeaturePlot(object, features = "Hypo") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "Hyper") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")





# set up reference
ref <- readRDS("/Users/tobymilne-clark/Library/Mobile Documents/com~apple~CloudDocs/PhD/Data/Prostate_scRNA/scRNA/scRNAseq.exampleData.20220606.RDS")
ref <- UpdateSeuratObject(ref)

# extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
cluster <- droplevels(as.factor(ref$celltype))
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- spacexr::Reference(round(counts,0), cluster, round(nUMI,0))

# set up query with the RCTD function SpatialRNA
counts_hd <- object[["sketch"]]$counts
object_cells_hd <- colnames(object[["sketch"]])
coords <- GetTissueCoordinates(object)[object_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))


RCTD <- spacexr::create.RCTD(query, reference)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
object <- AddMetaData(object, metadata = RCTD@results$results_df)

p1 <- SpatialDimPlot(object, group.by = "first_type")
p2 <- SpatialDimPlot(object, group.by = "second_type")
p1 | p2

object <- AddModuleScore(
  object = object,
  features = list(neuronal_genes),
  name = "NeuronalScore",
  ctrl = 50         # number of control genes per bin (default)
)


SpatialFeaturePlot(object, features = "NeuronalScore1") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")


collagen_genes <- c("ADAMTS2","COL5A1","COL9A3","COL11A1","COL18A1","COL21A1","P4HB","PXDN")



key_protiens <- read_csv("/Users/tobymilne-clark/Downloads/key_proteins 1.txt")



object <- AddModuleScore(
  object = object,
  features = list(key_protiens$Gene),
  name = "ProtienScore",
  ctrl = 50         # number of control genes per bin (default)
)

SpatialFeaturePlot(object, features = "ProtienScore1") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")


VEGF_keyprotein_list <- c("KDR", "RAC1", "CALM1", "NOS3", "PRKACG", "CAV1", "PRKACB", "HRAS",
                          "PRKACA", "FYN", "ACTG1", "ACTB", "AKT1", "MTOR", "PTK2", "RHOA",
                          "SRC", "CTNNB1", "KRAS", "MAPK14", "PIK3CA", "VEGFA", "CDC42", "HSP90AA1")


immune_keyprotein_list <- c("IL10", "HSPA8", "MMP2", "MYC", "AKT1", "IL4", "NANOG", "CCND1", "FOS", "IL1B", "IL6", "CCL2", "HGF", "FN1", "FOXO1", "CXCL8", "MMP9", "POMC",
                            "TP53", "SOX2", "MUC1", "TNF", "STAT3", "GATA3", "VEGFA", "ITGB1",  "PTGS2", "FGF2", "HSP90AA1")



object <- AddModuleScore(
  object = object,
  features = list(VEGF_keyprotein_list, immune_keyprotein_list),
  name = "HensonScore",
  ctrl = 50         # number of control genes per bin (default)
)

SpatialFeaturePlot(object, features = "HensonScore1") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")

SpatialFeaturePlot(object, features = "HensonScore2") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")


Prolaris <- read_excel("/Users/tobymilne-clark/Documents/Documents/CANCAP03/RNA/Gene_lists/Prolaris_genes.xlsx")


object <- AddModuleScore(
  object = object,
  features = list(Prolaris$`Gene name`),
  name = "ProlarisScore",
  ctrl = 50         # number of control genes per bin (default)
)

SpatialFeaturePlot(object, features = "ProlarisScore1") +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "right")


#### SpaCET ####

# set the path to the in-house breast cancer ST data. 
# User can set the paths to their own data.
visiumPath <- file.path("/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/VisiumHD/binned_outputs/square_016um")

# load ST data to create an SpaCET object.
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)

# filter out spots with less than 100 expressed genes
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes=100)

# plot the QC metrics
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "QualityControl", 
  spatialFeatures=c("UMI","Gene"),
  imageBg = TRUE,
  imageSize = "CompleteImage"
)
# This takes some time
SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType="PRAD", coreNo=6)

saveRDS(SpaCET_obj, file = "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Processed_data/VisiumHD_SpaCET_deconv.RDS")
# load the SpaCET object
SpaCET_obj <- readRDS("/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Processed_data/VisiumHD_SpaCET_deconv.RDS")
SpaCET_obj@results$deconvolution$propMat[1:13,1:6]

SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures=c("Malignant","Macrophage")
)

SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures="All", 
  sameScaleForFraction = TRUE,
  pointSize = 0.1, 
  nrow=5
)

# calculate the cell-cell colocalization.
SpaCET_obj <- SpaCET.CCI.colocalization(SpaCET_obj)

# visualize the cell-cell colocalization.
SpaCET.visualize.colocalization(SpaCET_obj)

# calculate the L-R network score across ST spots.
SpaCET_obj <- SpaCET.CCI.LRNetworkScore(SpaCET_obj,coreNo=6)

saveRDS(SpaCET_obj, file = "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Processed_data/VisiumHD_SpaCET_deconv.RDS")

# visualize the L-R network score.
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "LRNetworkScore", 
  spatialFeatures=c("Network_Score","Network_Score_pv")
)

# Ligand-Receptor analysis for a co-localized cell-type pair
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("CAF","Macrophage"))
## [1] "CAF and Macrophage M2 have potential intercellular interaction in the current tissue."

# Visualize the interaction analysis of a co-localized cell-type pair.
SpaCET.visualize.cellTypePair(SpaCET_obj, cellTypePair=c("CAF","Macrophage"))


# Identify the Tumor-Stroma Interface
SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj)

saveRDS(SpaCET_obj, file = "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Processed_data/VisiumHD_SpaCET_deconv.RDS")

SpaCET.visualize.spatialFeature(SpaCET_obj, 
                                spatialType = "Interface", 
                                spatialFeature = "Interface",
                                pointSize = 0.1)


SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures=c("Malignant","Macrophage"),
  pointSize = 0.1
)

SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures=c("Malignant","CAF"),
  pointSize = 0.1
)

SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "LRNetworkScore", 
  spatialFeatures=c("Network_Score","Network_Score_pv"),
  pointSize = 0.1
)

# Combine the interface and interaction spots
SpaCET_obj <- SpaCET.combine.interface(SpaCET_obj, cellTypePair=c("CAF","Macrophage"))

# Visualize the Interface. The spatialFeature should be formated as "Interface&celltype1_celltype2". Celltype 1 and 2 is in the order of alphabet.
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeature = "Interface&CAF_Macrophage", pointSize = 0.1)
SpaCET.distance.to.interface(SpaCET_obj, cellTypePair=c("CAF", "Macrophage"))


SpaCET_obj <- SpaCET.deconvolution.malignant(SpaCET_obj, coreNo=6)



####################
##### Giotto #######
####################


library(Giotto)
library(usethis)
usethis::edit_r_environ()
# if data.table corropt run .rs.restartR()
# 1. set working directory
results_folder <- "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Processed_data/results_Giotto_VisiumHD"

instructions <- createGiottoInstructions(
  save_dir = results_folder,
  save_plot = TRUE,
  show_plot = FALSE,
  return_plot = FALSE,
  plot_format = "pdf",
  dpi = 600
)

data_path <- "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/VisiumHD"

visiumHD_dir <- file.path(data_path, "binned_outputs/square_002um")

expression_path <- file.path(
  data_path,
  "/binned_outputs/square_002um/raw_feature_bc_matrix"
)


# 1 for ensembl gene ID, 2 for gene name
expr_results <- get10Xmatrix(
  path_to_data = expression_path,
  gene_column_index = 2
)

tissue_positions_path <- file.path(
  data_path,
  "/binned_outputs/square_002um/spatial/tissue_positions.parquet"
)

tissue_positions <- data.table::as.data.table(arrow::read_parquet(tissue_positions_path))

# convert expression matrix to minimal data.frame or data.table object
matrix_tile_dt <- data.table::as.data.table(Matrix::summary(expr_results))
genes <- expr_results@Dimnames[[1]]
samples <- expr_results@Dimnames[[2]]
matrix_tile_dt[, gene := genes[i]]
matrix_tile_dt[, pixel := samples[j]]

# merge data.table matrix and spatial coordinates to create input for Giotto Polygons
expr_pos_data <- data.table::merge.data.table(matrix_tile_dt,
                                              tissue_positions,
                                              by.x = "pixel",
                                              by.y = "barcode"
)

expr_pos_data <- expr_pos_data[, .(pixel, pxl_row_in_fullres, pxl_col_in_fullres, gene, x)]
colnames(expr_pos_data) <- c("pixel", "x", "y", "gene", "count")

# colnames(expr_pos_data) <- c("pixel", "y", "x", "gene", "count")

giotto_points <- createGiottoPoints(x = expr_pos_data[, .(x, y, gene, pixel, count)])

# create giotto polygons, here we create hexagons
hexbin400 <- tessellate(
  extent = ext(giotto_points),
  shape = "hexagon",
  shape_size = 400,
  name = "hex400"
)

plot(hexbin400)


# gpoints provides spatial gene expression information
# gpolygons provides spatial unit information (here = hexagon tiles)
visiumHD <- createGiottoObjectSubcellular(
  gpoints = list("rna" = giotto_points),
  gpolygons = list("hex400" = hexbin400),
  instructions = instructions
)

# create spatial centroids for each spatial unit (hexagon)
visiumHD <- addSpatialCentroidLocations(
  gobject = visiumHD,
  poly_info = "hex400"
)

rm(expr_pos_data, matrix_tile_dt, tissue_positions, expr_results, hexbin400)

feature_data <- fDataDT(visiumHD)

spatInSituPlotPoints(visiumHD,
                     show_image = FALSE,
                     feats = list("rna" = feature_data$feat_ID[10:20]),
                     show_legend = TRUE,
                     spat_unit = "hex400",
                     point_size = 0.25,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = "hex400",
                     polygon_bg_color = NA,
                     polygon_color = "white",
                     polygon_line_size = 0.1,
                     expand_counts = TRUE,
                     count_info_column = "count",
                     jitter = c(25, 25)
)


# calculate overlap between points and polygons
visiumHD <- calculateOverlap(visiumHD,
                             spatial_info = "hex400",
                             feat_info = "rna"
)

showGiottoSpatialInfo(visiumHD)


# convert overlap results to bin by gene matrix
visiumHD <- overlapToMatrix(visiumHD,
                            poly_info = "hex400",
                            feat_info = "rna",
                            name = "raw"
)

# this action will automatically create an active spatial unit, e.g. hexbin 400
activeSpatUnit(visiumHD)


hexbin100 <- tessellate(
  extent = ext(visiumHD),
  shape = "hexagon",
  shape_size = 100,
  name = "hex100"
)

visiumHD <- setPolygonInfo(
  gobject = visiumHD,
  x = hexbin100,
  name = "hex100",
  initialize = TRUE
)

visiumHD <- addSpatialCentroidLocations(
  gobject = visiumHD,
  poly_info = "hex100"
)

activeSpatUnit(visiumHD) <- "hex100"

spatInSituPlotPoints(visiumHD,
                     show_image = FALSE,
                     feats = list("rna" = feature_data$feat_ID[1:20]),
                     show_legend = TRUE,
                     spat_unit = "hex100",
                     point_size = 0.1,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = "hex100",
                     polygon_bg_color = NA,
                     polygon_color = "white",
                     polygon_line_size = 0.2,
                     expand_counts = TRUE,
                     count_info_column = "count",
                     jitter = c(25, 25)
)


visiumHD <- calculateOverlap(visiumHD,
                             spatial_info = "hex100",
                             feat_info = "rna"
)

visiumHD <- overlapToMatrix(visiumHD,
                            poly_info = "hex100",
                            feat_info = "rna",
                            name = "raw"
)

visiumHD <- filterGiotto(visiumHD,
                         spat_unit = "hex100",
                         expression_threshold = 1,
                         feat_det_in_min_cells = 10,
                         min_det_feats_per_cell = 10
)

visiumHD <- normalizeGiotto(visiumHD,
                            spat_unit = "hex100",
                            scalefactor = 1000,
                            verbose = TRUE
)

visiumHD <- addStatistics(visiumHD,
                          spat_unit = "hex100"
)


visiumHD <- calculateHVF(visiumHD,
                         spat_unit = "hex100",
                         zscore_threshold = 1
)

visiumHD <- runPCA(visiumHD,
                   spat_unit = "hex100",
                   expression_values = "normalized",
                   feats_to_use = "hvf"
)

plotPCA(visiumHD,
        spat_unit = "hex100"
)


visiumHD <- runUMAP(visiumHD,
                    spat_unit = "hex100",
                    dimensions_to_use = 1:14,
                    n_threads = 10
)

# plot UMAP, coloring cells/points based on nr_feats
plotUMAP(
  gobject = visiumHD,
  spat_unit = "hex100",
  point_size = 2
)


# sNN network (default)
visiumHD <- createNearestNetwork(visiumHD,
                                 spat_unit = "hex100",
                                 dimensions_to_use = 1:14,
                                 k = 5
)

## leiden clustering ####
visiumHD <- doLeidenClusterIgraph(visiumHD,
                                  spat_unit = "hex100",
                                  resolution = 0.2,
                                  n_iterations = 1000
)

plotUMAP(
  gobject = visiumHD,
  spat_unit = "hex100",
  cell_color = "leiden_clus",
  point_size = 1.5,
  show_NN_network = T,
  edge_alpha = 0.05
)

spatInSituPlotPoints(visiumHD,
                     show_image = FALSE,
                     feats = NULL,
                     show_legend = FALSE,
                     spat_unit = "hex100",
                     point_size = 0.5,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = "hex100",
                     polygon_fill_as_factor = TRUE,
                     polygon_fill = "leiden_clus",
                     polygon_color = "black",
                     polygon_line_size = 0.3
)

####### Deconvolution ######
library(Seurat)
ref <- readRDS("/Users/tobymilne-clark/Library/Mobile Documents/com~apple~CloudDocs/PhD/Data/Prostate_scRNA/scRNA/scRNAseq.exampleData.20220606.RDS")
ref <- UpdateSeuratObject(ref)

Giotto_Prostate_reference <- seuratToGiottoV5(
  ref,
  spatial_assay = "RNA",
  dim_reduction = c("pca", "umap"),
  subcellular_assay = "RNA",
  sp_network = NULL,
  nn_network = NULL,
  verbose = TRUE
)

Giotto_Prostate_reference <- normalizeGiotto(Giotto_Prostate_reference)

markers_scran <- findMarkers_one_vs_all(gobject = Giotto_Prostate_reference, 
                                        method = "scran",
                                        expression_values = "normalized",
                                        cluster_column = "celltype", 
                                        min_feats = 10,
                                        pval = 0.01,
                                        logFC = 0.5)

top_markers <- markers_scran[, head(.SD, 25), by = "cluster"]$feats


markers <- FindAllMarkers(ref, group.by = "celltype", logfc.threshold = 0.5)
markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(gene %in% visiumHD@feat_metadata$hex100$rna$feat_ID) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup() -> top10

DoHeatmap(ref,slot = "data",  features = top10$gene) + NoLegend()


sign_matrix <- makeSignMatrixDWLSfromMatrix(
  matrix = getExpression(Giotto_Prostate_reference,
                         values = "normalized",
                         output = "matrix"),
  cell_type = droplevels(pDataDT(Giotto_Prostate_reference)$celltype),
  sign_gene = top10$gene)


visiumHD <- runDWLSDeconv(visiumHD,
                   sign_matrix = sign_matrix
)

top10$gene %in% visiumHD@feat_metadata$hex100$rna$feat_ID

spatDeconvPlot(visiumHD, 
               show_image = FALSE,
               radius = 50
)


spatial_markers <- findMarkers_one_vs_all(visiumHD,
                       cluster_column = "leiden_clus", min_feats = 20
)

celltypes <- levels(factor(top10$cluster)) 

sign_list <- list()

for (i in 1:length(celltypes)){
  sign_list[[i]] = top10[which(top10$cluster == celltypes[i]),]$gene
}

sign_matrix <- makeSignMatrixPAGE(sign_names = celltypes,
                                  sign_list = sign_list)


visiumHD <- runPAGEEnrich(gobject = visiumHD, 
                   sign_matrix = sign_matrix,
                   min_overlap_genes = 1)

cell_types_PAGE <- colnames(sign_matrix)
library(viridis)
plotMetaDataCellsHeatmap(gobject = visiumHD,
                         metadata_cols = "leiden_clus",
                         value_cols = cell_types_PAGE,
                         spat_enr_names = "PAGE",
                         x_text_size = 8,
                         y_text_size = 8,
                         x_text_angle = 90)

spatCellPlot2D(gobject = visiumHD,
               spat_enr_names = "PAGE",
               cell_annotation_values = cell_types_PAGE[1:4],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               show_legend = TRUE)

spatCellPlot2D(gobject = visiumHD,
               spat_enr_names = "PAGE",
               cell_annotation_values = cell_types_PAGE[5:8],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               show_legend = TRUE)

spatCellPlot2D(gobject = visiumHD,
               spat_enr_names = "PAGE",
               cell_annotation_values = cell_types_PAGE[9:12],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               show_legend = TRUE)


visiumHD_subset <- subsetGiottoLocs(
  gobject = visiumHD,
  x_min = 12000,
  x_max = 18000,
  y_min = 10000,
  y_max = 15000
)

spatInSituPlotPoints(visiumHD_subset,
                     spat_unit = "hex100",
                     show_image = FALSE,
                     feats = NULL,
                     show_legend = FALSE,
                     point_size = 0.5,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = "hex100",
                     polygon_fill_as_factor = TRUE,
                     polygon_fill = "leiden_clus",
                     polygon_color = "black",
                     polygon_line_size = 0.3
)

spatCellPlot2D(gobject = visiumHD_subset,
               spat_enr_names = "PAGE",
               cell_annotation_values = cell_types_PAGE[1:4],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               show_legend = TRUE)

spatCellPlot2D(gobject = visiumHD_subset,
               spat_enr_names = "PAGE",
               cell_annotation_values = cell_types_PAGE[5:8],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               show_legend = TRUE)

spatCellPlot2D(gobject = visiumHD_subset,
               spat_enr_names = "PAGE",
               cell_annotation_values = cell_types_PAGE[9:12],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               show_legend = TRUE)


top10 %>% dplyr::filter(cluster == "OE2") %>% dplyr::select(gene) -> OE2_genes

OE2_genes$gene

#Adrenergic fibers - Slc16a2, Slc18a2, Dbh
#Cholinergic fibers - slc18a3, chat
#Immune - RUNX1, RUNX3, FOXP3 



dataPath <- file.path(system.file(package = "SpaCET"), "extdata/GeneSets/")
CCS <- SpaCET::read.gmt(paste0(dataPath,"/","CancerCellState",".gmt"))


sign_matrix <- makeSignMatrixPAGE(sign_names = names(CCS),
                                  sign_list = CCS )


visiumHD <- runPAGEEnrich(gobject = visiumHD, 
                          sign_matrix = sign_matrix,
                          min_overlap_genes = 1,
                          name = "CCS_enrich")

cell_types_CCS<- colnames(sign_matrix)

plotMetaDataCellsHeatmap(gobject = visiumHD,
                         metadata_cols = "cell_types",
                         value_cols = cell_types_CCS,
                         spat_enr_names = "CCS_enrich",
                         x_text_size = 8,
                         y_text_size = 8)


spatCellPlot2D(gobject = visiumHD,
               spat_enr_names = "CCS_enrich",
               cell_annotation_values = cell_types_CCS[1:4],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               show_legend = TRUE)

spatCellPlot2D(gobject = visiumHD,
               spat_enr_names = "CCS_enrich",
               cell_annotation_values = cell_types_CCS[5:8],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               show_legend = TRUE)

spatCellPlot2D(gobject = visiumHD,
               spat_enr_names = "CCS_enrich",
               cell_annotation_values = cell_types_CCS[9:12],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               show_legend = TRUE)

spatCellPlot2D(gobject = visiumHD,
               spat_enr_names = "CCS_enrich",
               cell_annotation_values = cell_types_CCS[13:16],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               show_legend = TRUE)




neuronal_genes <- c("NRXN1", "NRXN3", "NTRK3", "SLITRK3", "SLITRK4", "GRIN2D", "GRIA3", "SHANK2", "CAMK2A", "CAMK2B", "CAMKK2", "GABRA5", "GABRB1", "GABRG3", "SLC32A1", "SYT7", "ADCY1", "RASGRF1", "CHAT", "KIF17", "PPFIBP2")
Rho_GTPase_genes <- c(  "ARHGAP25", "ARHGAP26", "ARHGEF7", "BAIAP2", "BCR", "EPHA2",   "FAM13A", "FARP1", "ITGB1", "MCF2", "PAK6", "PLEKHG3",   "SRGAP3", "TAOK3", "VAV2", "VAV3")
MeCP2 <- c("HDAC1", "HIPK2", "TBL1XR1", "TNRC6C")
notch <- c("HDAC1", "HDAC4", "HDAC9", "MAML2", "NCOR2", "RBPJ", "TBL1XR1")

Meth_genes <- list(neuronal_genes,
                   Rho_GTPase_genes) 

names(Meth_genes) <- c("neuronal",
                       "Rho_GTPase")

Meth_genes <- c(Meth_genes ,genes_list3)
sign_matrix_meth <- makeSignMatrixPAGE(sign_names = names(Meth_genes),
                                  sign_list = Meth_genes )


visiumHD <- runPAGEEnrich(gobject = visiumHD, 
                          sign_matrix = sign_matrix_meth,
                          min_overlap_genes = 1,
                          name = "Meth_enrich")

cell_types_Meth <- colnames(sign_matrix_meth)

spatCellPlot2D(gobject = visiumHD,
               spat_enr_names = "Meth_enrich",
               cell_annotation_values = cell_types_Meth[1:4],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_style = "sequential",
               cell_color_gradient = c("darkblue","blue" ,"white", "yellow", "orange", "darkred"),
               show_legend = TRUE)


plotMetaDataCellsHeatmap(gobject = visiumHD,
                         metadata_cols = "cell_types",
                         value_cols = cell_types_Meth,
                         spat_enr_names = "Meth_enrich",
                         x_text_size = 8,
                         y_text_size = 8)

plotMetaDataCellsHeatmap(gobject = visiumHD,
                         metadata_cols = "cell_types",
                         value_cols = celltypes,
                         spat_enr_names = "PAGE",
                         x_text_size = 8,
                         y_text_size = 8)

library(tidyverse)
expressed_DMRs <- read_csv("/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Metadata/DMRs_annotated_expressed.csv")

expressed_DMRs <- expressed_DMRs %>% mutate(direction = if_else(difference > 0, "Hyper", "Hypo")) 

genes_list3 <- expressed_DMRs %>%
  group_by(direction) %>%
  summarise(genes = list(annot.symbol), .groups = "drop") %>%
  deframe()


sign_matrix_meth <- makeSignMatrixPAGE(sign_names = names(genes_list3),
                                       sign_list = genes_list3 )

visiumHD <- runPAGEEnrich(gobject = visiumHD, 
                          sign_matrix = sign_matrix_meth,
                          min_overlap_genes = 1,
                          name = "Meth_enrich")

cell_types_Meth <- colnames(sign_matrix_meth)

plotMetaDataCellsHeatmap(gobject = visiumHD,
                         metadata_cols = "cell_types",
                         value_cols = cell_types_Meth,
                         spat_enr_names = "Meth_enrich",
                         x_text_size = 8,
                         y_text_size = 8)

spatCellPlot2D(gobject = visiumHD,
               spat_enr_names = "Meth_enrich",
               cell_annotation_values = cell_types_Meth[1:2],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_midpoint = 0,
               gradient_style = "sequential",
               cell_color_gradient = c("darkblue","blue" ,"white", "yellow", "orange", "darkred"),
               show_legend = TRUE)

?spatCellPlot2D

#### subset down to main celltypes and annotate


spatPlot(visiumHD,
         group_by = 'leiden_clus', group_by_subset = c(1,3,4,6,7,9,10,11,14,25),
         cell_color = 'leiden_clus', cow_n_col = 4,
         point_size = 1, axis_text = 6, axis_title = 6, legend_text = 6)


metadata = pDataDT(visiumHD)
subset_cell_IDs = metadata[leiden_clus %in% c(1,3,6,7,9,10,11,14)]$cell_ID
visiumHD_small = subsetGiotto(visiumHD, cell_ids = subset_cell_IDs)

spatDimPlot(visiumHD_small,
            cell_color = 'leiden_clus', 
            plot_alignment = 'horizontal', spat_point_size = 1)


spatInSituPlotPoints(visiumHD_small,
                     show_image = FALSE,
                     feats = NULL,
                     show_legend = FALSE,
                     spat_unit = "hex100",
                     point_size = 0.5,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = "hex100",
                     polygon_fill_as_factor = TRUE,
                     polygon_fill = "leiden_clus",
                     polygon_color = "black",
                     polygon_line_size = 0.3
)

pDataDT(visiumHD_small)
names = c(1,3,6,7,9,10,11,14)
feats = c('Epithelial',"Stroma","Non tumour epithelial","Tumour proximinal fibroblast","Immune","Endothelial","Tumour", "Epithelial CHGA+CHGB+")

annot_dt = data.table::data.table('names' = names, 'feats' = feats)
annot_dt[, annotname := paste0(feats)]
cell_annot = annot_dt$annotname;names(cell_annot) = annot_dt$names
visiumHD_small = annotateGiotto(visiumHD_small, annotation_vector = cell_annot, cluster_column = 'leiden_clus')

pDataDT(visiumHD_small)
spatInSituPlotPoints(visiumHD_small,
                     show_image = FALSE,
                     feats = NULL,
                     show_legend = TRUE,
                     spat_unit = "hex100",
                     point_size = 0.5,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = "hex100",
                     polygon_fill_as_factor = TRUE,
                     polygon_fill = "cell_types",
                     polygon_color = "black",
                     polygon_line_size = 0.3
)

#### cell-cell interactions

metadata = pDataDT(visiumHD_small)
uniq_clusters = length(unique(metadata$leiden_clus))

clusters_cell_types = paste0('cell ', 1:uniq_clusters)
names(clusters_cell_types) = 1:uniq_clusters

visiumHD = annotateGiotto(gobject = visiumHD, 
                          annotation_vector = clusters_cell_types, 
                          cluster_column = 'leiden_clus', 
                          name = 'cell_types')

# create network (required for binSpect methods)
visiumHD_small = createSpatialNetwork(gobject = visiumHD_small, minimum_k = 2)

# identify genes with a spatial coherent expression profile
km_spatialgenes = binSpect(visiumHD_small, bin_method = 'kmeans')

set.seed(seed = 2841)
cell_proximities = cellProximityEnrichment(gobject = visiumHD_small,
                                           cluster_column = 'cell_types',
                                           spatial_network_name = 'Delaunay_network',
                                           adjust_method = 'fdr',
                                           number_of_simulations = 1000)
# visualize results

# barplot
cellProximityBarplot(gobject = visiumHD_small, 
                     CPscore = cell_proximities, 
                     min_orig_ints = 3, min_sim_ints = 3)

# heatmap
cellProximityHeatmap(gobject = visiumHD_small, 
                     CPscore = cell_proximities, 
                     order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), 
                     color_names = c('blue', 'white', 'red'))

# network
cellProximityNetwork(gobject = visiumHD_small, 
                     CPscore = cell_proximities, 
                     remove_self_edges = T, only_show_enrichment_edges = T)


# network with self-edges
cellProximityNetwork(gobject = visiumHD_small, 
                     CPscore = cell_proximities,
                     remove_self_edges = F, self_loop_strength = 0.3,
                     only_show_enrichment_edges = F,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_weight_range_depletion = c(1, 2),
                     edge_weight_range_enrichment = c(2,5))


unique(cell_proximities)

# Option 1
spec_interaction = "Epithelial CHGA+CHGB+--Tumour proximinal fibroblast"
cellProximitySpatPlot2D(gobject = visiumHD_small,
                        interaction_name = spec_interaction,
                        show_network = T,
                        cluster_column = 'cell_types',
                        cell_color = 'cell_types',
                        cell_color_code = c('Epithelial CHGA+CHGB+' = 'lightblue', 'Tumour proximinal fibroblast' = 'red'),
                        point_size_select = 4, point_size_other = 2)


# Option 2: create additional metadata
spec_interaction = "Epithelial CHGA+CHGB+--Tumour proximinal fibroblast"

visiumHD_small = addCellIntMetadata(visiumHD_small,
                                      spatial_network = 'Delaunay_network',
                                      cluster_column = 'cell_types',
                                      cell_interaction = spec_interaction,
                                      name = 'Epi_TPF_interactions')
spatPlot(visiumHD_small, cell_color = 'Epi_TPF_interactions', legend_symbol_size = 3, point_size = 1,
         select_cell_groups =  c('other_Epithelial CHGA+CHGB+', 'other_Tumour proximinal fibroblast', 'select_Epithelial CHGA+CHGB+', 'select_Tumour proximinal fibroblast'))



spec_interaction = "Epithelial CHGA+CHGB+--Tumour"

visiumHD_small = addCellIntMetadata(visiumHD_small,
                                    spatial_network = 'Delaunay_network',
                                    cluster_column = 'cell_types',
                                    cell_interaction = spec_interaction,
                                    name = 'Epi_Tumour_interactions')
spatPlot(visiumHD_small, cell_color = 'Epi_Tumour_interactions', legend_symbol_size = 3, point_size = 1,
         select_cell_groups =  c('other_Epithelial CHGA+CHGB+', 'other_Tumour', 'select_Epithelial CHGA+CHGB+', 'select_Tumour'))





spatial_markers <- findMarkers_one_vs_all(visiumHD_small,
                                          cluster_column = "cell_types", min_feats = 20
)
top_markers <- spatial_markers[, head(.SD, 25), by = "cluster"]



spatial_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_min(ranking, n = 25) %>%
  dplyr::ungroup() -> top25_spatial

spatial_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_head(n = 25) %>%
  dplyr::ungroup() -> top25_spatial
slice_min()



plotHeatmap(
  visiumHD_small,
  feats =  neuronal_genes,
  cluster_column = "cell_types",
  legend_nrows = 3
)



library("fgsea")
library("msigdbr")
gene_sets_df <- msigdbr(species = 'Homo sapiens', category = 'H')

gene_sets <- gene_sets_df %>%
  split(x = .$gene_symbol, f = .$gs_name)

tumour_df <- spatial_markers %>% dplyr::filter(cluster == "Tumour")
rankings <- sign(tumour_df$logFC)*(-log10(tumour_df$p.value)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- tumour_df$feats # genes as names#

rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking

GSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

head(GSEAres[order(pval), ],10)

GSEA_tumour <- GSEAres %>% slice_min(padj,n = 15) %>% 
  mutate(direction = ifelse(NES > 0, "Positive", "Negative")) %>%
  ggplot() + geom_bar(mapping = aes(x = reorder(pathway, -NES), y = NES, fill = direction), stat = "identity") +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position="none")

ggsave(GSEA_tumour, filename = "/Users/tobymilne-clark/Documents/Data/Plots/GSEA_tumour.png", units = "in", width = 7, height = 7)




tpf_df <- spatial_markers %>% dplyr::filter(cluster == "Tumour proximinal fibroblast")
rankings <- sign(tpf_df$logFC)*(-log10(tpf_df$p.value)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- tpf_df$feats # genes as names#

rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking

GSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

head(GSEAres[order(pval), ],10)

GSEA_tpf <- GSEAres %>% slice_min(padj,n = 10) %>% 
  mutate(direction = ifelse(NES > 0, "Positive", "Negative")) %>%
  ggplot() + geom_bar(mapping = aes(x = reorder(pathway, -NES), y = NES, fill = direction), stat = "identity") +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position="none")

ggsave(GSEA_tpf, filename = "/Users/tobymilne-clark/Documents/Data/Plots/GSEA_TPF.png", units = "in", width = 7, height = 7)



Epi_df <- spatial_markers %>% dplyr::filter(cluster == "Epithelial CHGA+CHGB+")
rankings <- sign(Epi_df$logFC)*(-log10(Epi_df$p.value)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- Epi_df$feats # genes as names#

rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking

GSEAres <- fgsea(pathways = gene_sets, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

head(GSEAres[order(pval), ],10)

GSEA_Epi <- GSEAres %>% slice_min(padj,n = 10) %>% 
  mutate(direction = ifelse(NES > 0, "Positive", "Negative")) %>%
  ggplot() + geom_bar(mapping = aes(x = reorder(pathway, -NES), y = NES, fill = direction), stat = "identity") +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position="none")

ggsave(GSEA_Epi, filename = "/Users/tobymilne-clark/Documents/Data/Plots/GSEA_EPI_CDHA.png", units = "in", width = 7, height = 7)

























LR_data = data.table::fread(system.file("extdata", "human_ligand_receptors.txt", package = 'Giotto'))







featData <- fDataDT(visiumHD_small,
                    spat_unit = "hex100"
)

hvf_genes <- featData[hvf == "yes"]$feat_ID

visiumHD_small <- createSpatialNetwork(visiumHD_small,
                                 name = "kNN_network",
                                 spat_unit = "hex100",
                                 method = "kNN",
                                 k = 8
)

ranktest <- binSpect(visiumHD_small,
                     spat_unit = "hex100",
                     subset_feats = hvf_genes,
                     bin_method = "rank",
                     calc_hub = FALSE,
                     do_fisher_test = TRUE,
                     spatial_network_name = "kNN_network"
)


set0 <- ranktest[high_expr < 50][1:2]$feats
set1 <- ranktest[high_expr > 50 & high_expr < 100][1:2]$feats
set2 <- ranktest[high_expr > 100 & high_expr < 200][1:2]$feats
set3 <- ranktest[high_expr > 200 & high_expr < 400][1:2]$feats
set4 <- ranktest[high_expr > 400 & high_expr < 1000][1:2]$feats
set5 <- ranktest[high_expr > 1000][1:2]$feats

spatFeatPlot2D(visiumHD_small,
               spat_unit = "hex100",
               expression_values = "scaled",
               feats = c(set0, set1, set2),
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               cow_n_col = 2, point_size = 1
)

spatFeatPlot2D(visiumHD_small,
               spat_unit = "hex100",
               expression_values = "scaled",
               feats = c(set3, set4, set5),
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               cow_n_col = 2, point_size = 1
)


ext_spatial_genes <- ranktest[adj.p.value < 0.001]$feats

spat_cor_netw_DT <- detectSpatialCorFeats(visiumHD_small,
                                          spat_unit = "hex100",
                                          method = "network",
                                          spatial_network_name = "kNN_network",
                                          subset_feats = ext_spatial_genes
)


# cluster spatial genes
spat_cor_netw_DT <- clusterSpatialCorFeats(spat_cor_netw_DT,
                                           name = "spat_netw_clus",
                                           k = 16
)

# visualize clusters
heatmSpatialCorFeats(visiumHD_small,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = "spat_netw_clus",
                     heatmap_legend_param = list(title = NULL)
)




spatFeatPlot2D(visiumHD_small,
               spat_unit = "hex100",
               expression_values = "scaled",
               feats = c("SYP", "CHGA", "CHGB","ENO2"),
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               cow_n_col = 2, point_size = 0.5
)


# create metagene enrichment score for clusters
cluster_genes_DT <- showSpatialCorFeats(spat_cor_netw_DT,
                                        use_clus_name = "spat_netw_clus",
                                        show_top_feats = 1
)

cluster_genes <- cluster_genes_DT$clus
names(cluster_genes) <- cluster_genes_DT$feat_ID

visiumHD_small <- createMetafeats(visiumHD_small,
                            spat_unit = "hex100",
                            expression_values = "normalized",
                            feat_clusters = cluster_genes,
                            name = "cluster_metagene"
)

showGiottoSpatEnrichments(visiumHD_small)



spatCellPlot(visiumHD_small,
             spat_unit = "hex100",
             spat_enr_names = "cluster_metagene",
             gradient_style = "sequential",
             cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
             cell_annotation_values = as.character(c(1:4)),
             point_size = 1, cow_n_col = 2
)
spatCellPlot(visiumHD_small,
             spat_unit = "hex100",
             spat_enr_names = "cluster_metagene",
             gradient_style = "sequential",
             cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
             cell_annotation_values = as.character(c(5:8)),
             point_size = 1, cow_n_col = 2
)
spatCellPlot(visiumHD_small,
             spat_unit = "hex100",
             spat_enr_names = "cluster_metagene",
             gradient_style = "sequential",
             cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
             cell_annotation_values = as.character(c(9:12)),
             point_size = 1, cow_n_col = 2
)
spatCellPlot(visiumHD_small,
             spat_unit = "hex100",
             spat_enr_names = "cluster_metagene",
             gradient_style = "sequential",
             cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
             cell_annotation_values = as.character(c(13:16)),
             point_size = 1, cow_n_col = 2
)


balanced_genes <- getBalancedSpatCoexpressionFeats(
  spatCorObject = spat_cor_netw_DT,
  maximum = 5
)
selected_feats <- names(balanced_genes)

# give genes from same cluster same color
distinct_colors <- getDistinctColors(n = 20)
names(distinct_colors) <- 1:20

my_colors <- distinct_colors[balanced_genes]
names(my_colors) <- names(balanced_genes)

spatInSituPlotPoints(visiumHD_small,
                     spat_unit = "hex100",
                     show_image = FALSE,
                     feats = list("rna" = selected_feats),
                     feats_color_code = my_colors,
                     show_legend = FALSE,
                     point_size = 0.20,
                     show_polygon = FALSE,
                     use_overlap = FALSE,
                     polygon_feat_type = "hex100",
                     polygon_bg_color = NA,
                     polygon_color = "white",
                     polygon_line_size = 0.01,
                     expand_counts = TRUE,
                     count_info_column = "count",
                     jitter = c(25, 25)
)




















#######
# subset to smaller area

visiumHD_subset <- subsetGiottoLocs(
  gobject = visiumHD,
  x_min = 12000,
  x_max = 18000,
  y_min = 10000,
  y_max = 15000
)


spatFeatPlot2D(visiumHD_subset,
               spat_unit = "hex100",
               expression_values = "scaled",
               feats = c("SYP", "CHGA", "CHGB","ENO2"),
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               cow_n_col = 2, point_size = 1
)
spatFeatPlot2D(visiumHD,
               spat_unit = "hex100",
               expression_values = "scaled",
               feats = c("NRXN1","NRXN3"),
               gradient_style = "sequential",
               cell_color_gradient = c("blue", "white", "yellow", "orange", "red", "darkred"),
               cow_n_col = 2, point_size = 1
)



write_csv(pDataDT(visiumHD_small),
          "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Metadata/visiumHD_small_hex100_cell_metadata.csv"
)
write_csv(fDataDT(visiumHD_small),
          "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Metadata/visiumHD_small_hex100_feature_featdata.csv"
)

write_csv(spatial_markers,
          "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Metadata/visiumHD_small_hex100_spatial_markers.csv"
)





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

visiumHD <- runPAGEEnrich(gobject = visiumHD, 
                          sign_matrix = sign_matrix_KP,
                          min_overlap_genes = 1,
                          name = "Key_protein_enrich")

cell_types_KP <- colnames(sign_matrix_KP)

plotMetaDataCellsHeatmap(gobject = visiumHD,
                         metadata_cols = "cell_types",
                         value_cols = cell_types_KP,
                         spat_enr_names = "Key_protein_enrich",
                         x_text_size = 8,
                         y_text_size = 8)

spatCellPlot2D(gobject = visiumHD,
               spat_enr_names = "Key_protein_enrich",
               cell_annotation_values = cell_types_KP[1:6],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_midpoint = 0,
               gradient_style = "sequential",
               cell_color_gradient = c("darkblue","blue" ,"white", "yellow", "orange", "darkred"),
               show_legend = TRUE)



Tumour_KP <- read_csv("/Users/tobymilne-clark/Documents/Data/ICGC/PPI/key_proteins.txt")

key_proteins_tumour <- list(tumour = Tumour_KP$Gene)

sign_matrix_KP_tumour <- makeSignMatrixPAGE(sign_names = names(key_proteins_tumour),
                                       sign_list = key_proteins_tumour )


visiumHD <- runPAGEEnrich(gobject = visiumHD, 
                          sign_matrix = sign_matrix_KP_tumour,
                          min_overlap_genes = 1,
                          name = "Key_protein_tum_enrich")

cell_types_KP <- colnames(sign_matrix_KP)

plotMetaDataCellsHeatmap(gobject = visiumHD,
                         metadata_cols = "cell_types",
                         value_cols = cell_types_KP,
                         spat_enr_names = "Key_protein_enrich",
                         x_text_size = 8,
                         y_text_size = 8)

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

visiumHD <- runPAGEEnrich(gobject = visiumHD, 
                          sign_matrix = sign_matrix_Atlas,
                          min_overlap_genes = 1,
                          name = "Atlas_enrich")

cell_types_Atlas <- colnames(sign_matrix_Atlas)

plotMetaDataCellsHeatmap(gobject = visiumHD,
                         metadata_cols = "cell_types",
                         value_cols = cell_types_Atlas,
                         spat_enr_names = "Atlas_enrich",
                         x_text_size = 8,
                         y_text_size = 8)

spatCellPlot2D(gobject = visiumHD,
               spat_enr_names = "Atlas_enrich",
               cell_annotation_values = cell_types_Atlas[1:6],
               cow_n_col = 2,
               coord_fix_ratio = 1, 
               point_size = 1, 
               gradient_midpoint = 0,
               gradient_style = "sequential",
               cell_color_gradient = c("darkblue","blue" ,"white", "yellow", "orange", "darkred"),
               show_legend = TRUE)


saveGiotto(gobject = visiumHD_small,
           dir = "/Users/tobymilne-clark/Documents/Data/Atlas/Spatial_Transcriptomics/Processed_data/",
           foldername = "VisiumHD_small_hex100",
           method = "qs")












