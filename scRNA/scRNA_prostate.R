library(Seurat)
library(SeuratObject)
library(Matrix)
library(tidyverse)
library(reticulate)
library("SeuratDisk")
use_condaenv("r-reticulate")


scRNA <- readRDS("/Users/tobymilne-clark/Library/Mobile Documents/com~apple~CloudDocs/PhD/Data/Prostate_scRNA/scRNA/scRNAseq.exampleData.20220606.RDS")

scRNA <- UpdateSeuratObject(scRNA)

UMAP <- DimPlot(scRNA,label = T) + NoLegend()
ggplot2::ggsave(filename = "/Users/tobymilne-clark/Library/Mobile Documents/com~apple~CloudDocs/PhD/Data/Prostate_scRNA/umap_scRNAseq.tiff", plot = UMAP) # can add additional parameters similar to png

prostate.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
prostate.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

prostate.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top50


write_csv(top25,"/Users/tm739/Dropbox (Cambridge University)/Mac/Documents/MPhil/Data/Top25_scRNA_genes.csv")
plot <- DoHeatmap(object = scRNA, features = top50$gene) + NoLegend()
ggplot2::ggsave(filename = "/Users/tobymilne-clark/Library/Mobile Documents/com~apple~CloudDocs/PhD/Data/Prostate_scRNA/Top10_markers_scRNAseq.tiff", plot = plot) # can add additional parameters similar to png

library("scater")

scRNA_with_TSNE <- RunTSNE(
  scRNA,
  assay = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  reduction.key = "tSNE_",
)
TSNE <- TSNEPlot(scRNA_with_TSNE,label = T) + NoLegend()
ggplot2::ggsave(filename = "/Users/tobymilne-clark/Library/Mobile Documents/com~apple~CloudDocs/PhD/Data/Prostate_scRNA/tsne_scRNAseq.tiff", plot = TSNE) # can add additional parameters similar to png



adata <- convertFormat(scRNA, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
adata
??convertFormat

write_h5ad(adata, "/Users/tobymilne-clark/Documents/Data/Atlas/scRNA_Prostate_Charlie_ANNData.h5ad")

SaveH5Seurat(scRNA, filename = "/Users/tobymilne-clark/Documents/Data/Atlas/scRNA_Prostate_Charlie.h5Seurat")
Convert(
  "/Users/tobymilne-clark/Documents/Data/Atlas/scRNA_Prostate_Charlie.h5Seurat",
  dest = "h5ad"
)

# CellPhoneDB Protein Protein interaction networks for scRNAseq
setwd("/Users/tobymilne-clark/Documents/Data/Atlas/scRNA_prostate")

# data preparation
Idents(scRNA) = scRNA$celltype
# Save normalised counts - NOT scaled!
writeMM(scRNA@assays$RNA@data, file = 'prostate_counts_mtx/matrix.mtx')
# save gene and cell names
write(x = rownames(scRNA@assays$RNA@data), file = "prostate_counts_mtx/features.tsv")
write(x = colnames(scRNA@assays$RNA@data), file = "prostate_counts_mtx/barcodes.tsv")

scRNA@meta.data$Cell = rownames(scRNA@meta.data)
df = scRNA@meta.data[, c('Cell', 'celltype')]
write.table(df, file ='prostate_meta.tsv', sep = '\t', quote = F, row.names = F)

## Extract DEGs for each cell_type
DEGs <- FindAllMarkers(scRNA, 
                       test.use = 'LR', 
                       verbose = F, 
                       only.pos = T, 
                       random.seed = 1, 
                       logfc.threshold = 0.2, 
                       min.pct = 0.1, 
                       return.thresh = 0.05)

fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC > 0.1)

# 1st column = cluster; 2nd column = gene 
fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 
write.table(fDEGs, file ='prostate_DEGs.tsv', sep = '\t', quote = F, row.names = F)


DoHeatmap(object = scRNA, features = fDEGs$gene) + NoLegend()



DEGs_all <- FindAllMarkers(scRNA, 
                       test.use = 'LR', 
                       verbose = F, 
                       only.pos = F, 
                       random.seed = 1, 
                       logfc.threshold = 0.2, 
                       min.pct = 0.1, 
                       return.thresh = 0.05)
??findMarkers
fDEGs_all = subset(DEGs_all, p_val_adj < 0.05 & (avg_log2FC > 0.1 | avg_log2FC < -0.1))


DEGs_all_wilcoxon <- FindAllMarkers(scRNA, 
                           test.use = 'wilcox', 
                           verbose = F, 
                           only.pos = F, 
                           random.seed = 1, 
                           logfc.threshold = 0.2, 
                           min.pct = 0.1, 
                           return.thresh = 0.05)

fDEGs_all_wilcoxon = subset(DEGs_all_wilcoxon, p_val_adj < 0.05 & (avg_log2FC > 0.1 | avg_log2FC < -0.1))

fDEGs_all_wilcoxon %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top50

DoHeatmap(object = scRNA, features = top50$gene) + NoLegend()



write.table(fDEGs_all_wilcoxon, file ='/Users/tobymilne-clark/Documents/Data/Atlas/scRNA_prostate/prostate_DEGs_all_wilcox_seurat.tsv', sep = '\t', quote = F, row.names = F)



DEGs_all %>% 
  mutate(colour = if_else(avg_log2FC > 0.1 & p_val_adj < 0.05, "up_regualted",if_else(avg_log2FC < -0.1 & p_val_adj < 0.05,"down_regulated","non_DEG"))) %>% 
  mutate(log10pval = if_else(log10(p_val_adj) < -10, -10, log10(p_val_adj))) %>% 
  ggplot() + 
  geom_point(aes(x = avg_log2FC, y = -log10pval, color = colour),alpha = 0.4)  + 
  scale_colour_manual(values=c("dodgerblue","grey","lightcoral")) + theme_minimal() +
  labs(x = "log2FC", y = "-log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  geom_vline(xintercept = -0.1, linetype = "dashed") +
  xlim (-15,15)
  
unique(fDEGs_all$gene)

library("Nebulosa")

plot_density(scRNA, "AR")




dat <- data.frame(scRNA@meta.data)
dat$UMAP1 <- Embeddings(scRNA, "umap")[,1]
dat$UMAP2 <- Embeddings(scRNA, "umap")[,2]

dat_bg <- dat[,-(which(colnames(dat)=="group"))]


density_plot <- ggplot(dat, aes(x=UMAP1, y=UMAP2)) +
  stat_density_2d(geom="raster", aes(fill=stat(ndensity)), contour=F) + #ndensity calculates the normalized density for each sample--otherwise density would be affected by the number of cells for each sample, which is variable
  geom_point(data=dat_bg, shape=16, size=0.1, alpha=0.2, color="white") +
  scale_fill_gradientn(colours=viridisLite::mako(100), name="Density") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~group, ncol=8) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, color="black"),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"))


ggplot2::ggsave(filename = "/Users/tobymilne-clark/Library/Mobile Documents/com~apple~CloudDocs/PhD/Data/Prostate_scRNA/density_scRNAseq.tiff", plot = density_plot) # can add additional parameters similar to png



p3 <- plot_density(scRNA, c("KLK3", "KLK4"), joint = TRUE)
p3 + plot_layout(ncol = 1)



plot_density(scRNA, "AR")


# CellChat
library(CellChat) 

data.input <- scRNA[["RNA"]]$data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(scRNA)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellChat <- createCellChat(object = scRNA, group.by = "ident", assay = "RNA")

CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)
cellChat@DB <- CellChatDB.use

devtools::install_github('immunogenomics/presto')
cellchat <- subsetData(cellChat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#options(future.globals.maxSize = 8000 * 1024^2)
cellchat <- computeCommunProb(cellchat, type = "triMean")

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

?netVisual_aggregate
pathways.show <- c("TGFb") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

netAnalysis_contribution(cellchat, signaling = pathways.show)

netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:12), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:12), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 12, targets.use = c(1:11), remove.isolate = FALSE)

netVisual_chord_gene(cellchat, sources.use = 5, targets.use = c(1:12), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat, sources.use = 6, targets.use = c(1:12), lab.cex = 0.5,legend.pos.y = 30)

netVisual_chord_gene(cellchat, sources.use = c(5,6), targets.use = c(1:12), slot.name = "netP", legend.pos.x = 10)


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

netAnalysis_signalingRole_scatter(cellchat)

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

library(ggalluvial)

netAnalysis_river(cellchat, pattern = "outgoing")


cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)


?netVisual_embedding


?netVisual_diffInteraction


gg2 <- compareInteractions(cellchat, show.legend = F, group.facet = "group", group.facet.levels = c("tumour","benign"), measure = "weight")

netVisual_diffInteraction(cellchat, weight.scale = T)




devtools::install_github("saeyslab/nichenetr")



devtools::install_github('cole-trapnell-lab/monocle3')




library("monocle")


#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(scRNA@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = scRNA@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

expressed_genes <- row.names(subset(fData(monocle_cds),
                                    num_cells_expressed >= 10))

plot_cell_clusters(monocle_cds, 1, 2, color = "group")


diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],
                                      fullModelFormulaStr = "~group")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

monocle_cds_tum <- setOrderingFilter(monocle_cds, ordering_genes)
plot_ordering_genes(monocle_cds_tum)

monocle_cds_tum <- reduceDimension(monocle_cds_tum, max_components = 2,
                            method = 'DDRTree')

orderCells_fix <- function (cds, root_state = NULL, num_paths = NULL, reverse = NULL) 
{
  if (class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  if (any(c(length(cds@reducedDimS) == 0, length(cds@reducedDimK) == 
            0))) {
    stop("Error: dimension reduction didn't prodvide correct results. Please check your reduceDimension() step and ensure correct dimension reduction are performed before calling this function.")
  }
  root_cell <- monocle:::select_root_cell(cds, root_state, reverse)
  cds@auxOrderingData <- new.env(hash = TRUE)
  if (cds@dim_reduce_type == "ICA") {
    if (is.null(num_paths)) {
      num_paths = 1
    }
    adjusted_S <- t(cds@reducedDimS)
    dp <- as.matrix(dist(adjusted_S))
    cellPairwiseDistances(cds) <- as.matrix(dist(adjusted_S))
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    minSpanningTree(cds) <- dp_mst
    next_node <<- 0
    res <- pq_helper(dp_mst, use_weights = FALSE, root_node = root_cell)
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    order_list <- monocle:::extract_good_branched_ordering(res$subtree, 
                                                 res$root, cellPairwiseDistances(cds), num_paths, 
                                                 FALSE)
    cc_ordering <- order_list$ordering_df
    row.names(cc_ordering) <- cc_ordering$sample_name
    minSpanningTree(cds) <- as.undirected(order_list$cell_ordering_tree)
    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), 
    ]$pseudo_time
    pData(cds)$State <- cc_ordering[row.names(pData(cds)), 
    ]$cell_state
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 
                                                        2)]$name
    minSpanningTree(cds) <- dp_mst
    cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree <- as.undirected(order_list$cell_ordering_tree)
  }
  else if (cds@dim_reduce_type == "DDRTree") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- monocle:::extract_ddrtree_ordering(cds, root_cell)
    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), 
    ]$pseudo_time
    K_old <- monocle:::reducedDimK(cds)
    old_dp <- monocle:::cellPairwiseDistances(cds)
    old_mst <- monocle:::minSpanningTree(cds)
    old_A <- monocle:::reducedDimA(cds)
    old_W <- monocle:::reducedDimW(cds)
    cds <- project2MST_fix(cds, monocle:::project_point_to_line_segment)
    monocle:::minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree
    root_cell_idx <- which(V(old_mst)$name == root_cell, 
                           arr.ind = T)
    cells_mapped_to_graph_root <- which(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == 
                                          root_cell_idx)
    if (length(cells_mapped_to_graph_root) == 0) {
      cells_mapped_to_graph_root <- root_cell_idx
    }
    cells_mapped_to_graph_root <- V(minSpanningTree(cds))[cells_mapped_to_graph_root]$name
    tip_leaves <- names(which(degree(minSpanningTree(cds)) == 
                                1))
    root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% 
                                              tip_leaves][1]
    if (is.na(root_cell)) {
      root_cell <- monocle:::select_root_cell(cds, root_state, reverse)
    }
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    cc_ordering_new_pseudotime <- monocle:::extract_ddrtree_ordering(cds, 
                                                           root_cell)
    pData(cds)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(cds)), 
    ]$pseudo_time
    if (is.null(root_state) == TRUE) {
      closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
      pData(cds)$State <- cc_ordering[closest_vertex[, 
                                                     1], ]$cell_state
    }
    monocle:::reducedDimK(cds) <- K_old
    monocle:::cellPairwiseDistances(cds) <- old_dp
    monocle:::minSpanningTree(cds) <- old_mst
    monocle:::reducedDimA(cds) <- old_A
    monocle:::reducedDimW(cds) <- old_W
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 
                                                        2)]$name
  }
  else if (cds@dim_reduce_type == "SimplePPT") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- monocle:::extract_ddrtree_ordering(cds, root_cell)
    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), 
    ]$pseudo_time
    pData(cds)$State <- cc_ordering[row.names(pData(cds)), 
    ]$cell_state
    mst_branch_nodes <- V(minSpanningTree(cds))[which(degree(minSpanningTree(cds)) > 
                                                        2)]$name
  }
  cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- mst_branch_nodes
  cds
}


project2MST_fix <- function (cds, Projection_Method) 
{
  dp_mst <- monocle:::minSpanningTree(cds)
  Z <- reducedDimS(cds)
  Y <- reducedDimK(cds)
  cds <- monocle:::findNearestPointOnMST(cds)
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- as.matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  tip_leaves <- names(which(degree(dp_mst) == 1))
  if (!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  }
  else {
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
    for (i in 1:length(closest_vertex)) {
      neighbors <- names(V(dp_mst)[suppressWarnings(.nei(closest_vertex_names[i], 
                                                        mode = "all"))])
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]
      for (neighbor in neighbors) {
        if (closest_vertex_names[i] %in% tip_leaves) {
          tmp <- monocle:::projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], 
                                            neighbor)])
        }
        else {
          tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], 
                                              neighbor)])
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, dist(rbind(Z_i, tmp)))
      }
      if (!is(projection, "matrix")) 
        projection <- as.matrix(projection)
      P[, i] <- projection[which(distance == min(distance))[1], 
      ]
    }
  }
  colnames(P) <- colnames(Z)
  dp <- as.matrix(dist(t(P)))
  min_dist = min(dp[dp != 0])
  dp <- dp + min_dist
  diag(dp) <- 0
  monocle:::cellPairwiseDistances(cds) <- dp
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
  cds
}

monocle_cds_tum <- orderCells_fix(monocle_cds_tum)

plot_cell_trajectory(monocle_cds_tum, color_by = "celltype")

plot_cell_trajectory(monocle_cds_tum, color_by = "State") 
plot_cell_trajectory(monocle_cds_tum, color_by = "State") +
  facet_wrap(~State, nrow = 1)

table(pData(monocle_cds_tum)$State, pData(monocle_cds_tum)$group)[,"normal"]

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$group)[,"normal"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
monocle_cds_tum <- orderCells_fix(monocle_cds_tum, root_state = GM_state(monocle_cds_tum))
plot_cell_trajectory(monocle_cds_tum, color_by = "Pseudotime")


diff_test_res_pseudo <- differentialGeneTest(monocle_cds_tum[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(monocle_cds_tum[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)


BEAM_res <- BEAM(monocle_cds_tum, branch_point = 4, cores = 4)
