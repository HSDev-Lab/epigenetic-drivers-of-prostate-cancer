#Relevant Fig1 panels

library(data.table)
library(ComplexHeatmap)
library(viridis)
library(ggrepel)
library(annotatr)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(tidyverse)


## ---- c: methylation heatmap -------------------------------------------------

meth_df  <- fread("data/ICGC_recommended_heatmap_methylation_beta.tsv")
metadata <- fread("data/ICGC_recommended_heatmap_methylation_metadata.tsv")

meth_df <- dplyr::mutate(meth_df, Sample_ID = paste0(Patient, "_", sample_type))

# rows = DMRs, columns = samples
meth_matrix <- meth_df %>%
  dplyr::select(-Patient, -sample_type) %>%
  column_to_rownames("Sample_ID") %>%
  as.matrix() %>%
  t()

annotation_df <- meth_df %>%
  dplyr::select(Patient, sample_type, Sample_ID) %>%
  dplyr::filter(Patient != "CRUK_PC_0050") %>%
  dplyr::left_join(metadata, by = c("Patient" = "ICGC_ID_short")) %>%
  dplyr::mutate(
    mevotype   = ifelse(sample_type != "TUM", "benign", mevotype),
    SPOP       = ifelse(sample_type != "TUM", "benign", SPOP),
    ETS_status = ifelse(sample_type != "TUM", "benign", ETS_status)
  ) %>%
  column_to_rownames("Sample_ID") %>%
  dplyr::select(sample_type, mevotype, SPOP, ETS_status) %>%
  na.omit()

sample_order <- annotation_df %>%
  dplyr::mutate(mevotype = factor(mevotype)) %>%
  dplyr::arrange(mevotype) %>%
  rownames()

col_ha <- HeatmapAnnotation(
  df  = annotation_df[sample_order, ],
  col = list(sample_type = c("TUM" = "#F15A2B", "BEN" = "#27A9E1", "BPH" = "#92ac90"))
)

ht_dmr <- Heatmap(
  meth_matrix[, sample_order],
  name = "Methylation",
  top_annotation = col_ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = plasma(100),
  use_raster = FALSE
)
ht_dmr
pdf("Fig1c_DMR_methylation_heatmap.pdf", width = 10, height = 8)
draw(ht_dmr)
dev.off()


## ---- e: SPOP promoter methylation -------------------------------------------

spop_meth <- fread("data/SPOP_methylation_by_mevo.tsv")
icgc_ppcg <- fread("data/ICGC_PPCG_ID.tsv") %>%
  dplyr::mutate(Sanger_ID_clean = sub("^(PD\\d+).*", "\\1", Sanger_ID))
ppcg_evotype <- fread("data/PPCG_evotype.csv")

spop_meth <- spop_meth %>%
  dplyr::left_join(icgc_ppcg) %>%
  dplyr::left_join(ppcg_evotype, by = c("Sanger_ID_clean" = "PPCG_ID")) %>%
  dplyr::mutate(evotype = ifelse(is.na(evotype), NA,
                                 ifelse(evotype == "1", "alternative", "canonical")))

spop_plot <- ggplot(spop_meth, aes(x = sample_type, y = methylation, group = ICGC_ID_short,
                                   color = mevotype, shape = evotype)) +
  geom_point(size = 3) +
  geom_line(alpha = 0.7) +
  scale_color_manual(values = c("mevotype_1" = "#986EAC", "mevotype_2" = "#5CAE63")) +
  scale_shape_manual(values = c("alternative" = 16, "canonical" = 1), na.value = 2) +
  scale_x_discrete(expand = c(0.1, 0)) +
  labs(x = "Sample type", y = "Methylation (%)", color = "Mevotype") +
  theme_minimal(base_size = 16) +
  theme(axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  geom_text_repel(data = dplyr::filter(spop_meth, sample_type == "BEN"),
                  aes(label = ICGC_ID_short), size = 4, nudge_x = 0.2,
                  direction = "y", segment.size = 0.2)

ggsave("Fig1e_SPOP_methylation_by_mevotype.pdf", spop_plot, width = 5, height = 7)


## ---- f, g: pathway enrichment of DMR-associated genes -----------------------

annotations <- build_annotations(genome = "hg38",
                                 annotations = c("hg38_basicgenes", "hg38_cpgs"))

annotate_genes <- function(regions) {
  gr <- makeGRangesFromDataFrame(regions[, c("chr", "start", "end")])
  df <- as.data.frame(annotate_regions(gr, annotations = annotations,
                                       ignore.strand = TRUE, quiet = TRUE))
  unique(na.omit(df$annot.symbol))
}

dmrs <- fread("data/ICGC_TMC_DMRs_all.tsv")
epic_panel <- fread("data/EPIC_panel_hg38.bed")
colnames(epic_panel) <- c("chr", "start", "end")

dmr_genes  <- annotate_genes(dmrs)
epic_genes <- annotate_genes(epic_panel)

# gene symbols -> Entrez; the EPIC array genes are the background
to_entrez <- function(genes) {
  clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)$ENTREZID
}

dmr_entrez      <- to_entrez(dmr_genes)
universe_entrez <- to_entrez(epic_genes)

kegg <- enrichKEGG(gene = dmr_entrez, organism = "hsa", pvalueCutoff = 0.05,
                   universe = universe_entrez)
reactome <- enrichPathway(gene = dmr_entrez, organism = "human", pvalueCutoff = 0.05,
                          readable = TRUE, universe = universe_entrez)

ggsave("Fig1f_KEGG_DMR_genes.pdf", dotplot(kegg, showCategory = 8), width = 6, height = 7)
ggsave("Fig1g_Reactome_DMR_genes.pdf", barplot(reactome, showCategory = 8), width = 6, height = 7)
