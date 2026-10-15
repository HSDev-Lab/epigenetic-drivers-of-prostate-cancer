## Relevant Extended Data Fig 3 panels

library(data.table)
library(VennDiagram)
library(pheatmap)
library(viridis)
library(annotatr)
library(tidyverse)


## ---- a: overlap of differentially methylated cytosines ----------------------

as_dmc <- function(file) {
  x <- fread(file)
  unique(paste0(x$V1, ":", x$V4, "-", x$V5))
}

dmc_list <- list(current_work = as_dmc("data/Supp_Table1_ICGC_DMR_CpG.tsv"),
                 Metilene     = as_dmc("data/metilene_DMR_CpG.tsv"),
                 methylKit    = as_dmc("data/methylKit_DMR_CpG.tsv"))

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
venn_plot <- venn.diagram(dmc_list, filename = NULL,
                          fill = c("#E41A1C", "#377EB8", "#4DAF4A"), alpha = 0.5,
                          cex = 2, cat.cex = 2, cat.pos = c(0, 0, 180),
                          cat.dist = c(0.05, 0.05, 0.05))

pdf("EFig3a_DMC_overlap.pdf", width = 6, height = 6)
grid::grid.draw(venn_plot)
dev.off()


## ---- b: TCGA-PRAD heatmaps --------------------------------------------------

clin <- fread("data/TCGA_PRAD_clindata.tsv")

# mean beta per DMR and sample
dmr_means <- function(file) {
  fread(file) %>%
    group_by(chr, start, end) %>%
    summarise(across(starts_with("TCGA"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    mutate(region = paste(chr, start, end, sep = ".")) %>%
    column_to_rownames("region") %>%
    select(-chr, -start, -end) %>%
    as.matrix()
}

tcga_heatmap <- function(mat, filename) {
  annotation_col <- data.frame(sample_type = clin$sample_type[match(colnames(mat), clin$barcode)],
                               row.names = colnames(mat))
  pheatmap(mat,
           cluster_rows = TRUE, cluster_cols = TRUE,
           annotation_col = annotation_col,
           annotation_colors = list(sample_type = c("Solid Tissue Normal" = "#27A9E1",
                                                    "Primary Tumor" = "#F15A2B",
                                                    "Metastatic" = "#8C564B")),
           show_rownames = FALSE, show_colnames = FALSE,
           color = plasma(100), legend = FALSE,
           filename = filename, width = 5, height = 5)
}

tcga_heatmap(dmr_means("data/PRAD.meth.mtx.hyperDMR.tsv"), "EFig3b_TCGA_hypermethylated_DMRs.pdf")
tcga_heatmap(dmr_means("data/PRAD.meth.mtx.hypoDMR.tsv"), "EFig3b_TCGA_hypomethylated_DMRs.pdf")


## ---- c, d: genomic annotation of the DMRs -----------------------------------

annotations <- build_annotations(genome = "hg38",
                                 annotations = c("hg38_basicgenes", "hg38_cpgs"))

annotate_panel <- function(regions) {
  gr <- makeGRangesFromDataFrame(regions[, c("chr", "start", "end")])
  as.data.frame(annotate_regions(gr, annotations = annotations, ignore.strand = TRUE, quiet = TRUE))
}

dmr_annot  <- annotate_panel(fread("data/ICGC_TMC_DMRs_all.tsv"))
epic_annot <- annotate_panel(setNames(fread("data/EPIC_panel_hg38.bed")[, 1:3], c("chr", "start", "end")))

# percentage of annotations in each category, within CpG and within gene annotations
annotation_percent <- function(prefix, labels) {
  count_by_type <- function(annot) {
    annot %>%
      filter(startsWith(annot.type, prefix)) %>%
      count(annot.type)
  }
  full_join(count_by_type(dmr_annot), count_by_type(epic_annot),
            by = "annot.type", suffix = c("_dmr", "_epic")) %>%
    mutate(`ICGC DMRs` = 100 * n_dmr / sum(n_dmr),
           `EPIC panel` = 100 * n_epic / sum(n_epic),
           annot.type = factor(annot.type, levels = names(labels), labels = labels)) %>%
    select(annot.type, `EPIC panel`, `ICGC DMRs`) %>%
    pivot_longer(-annot.type, names_to = "Panel", values_to = "Percent") %>%
    mutate(Panel = factor(Panel, levels = c("EPIC panel", "ICGC DMRs")))
}

annotation_barplot <- function(df, xlab) {
  ggplot(df, aes(x = annot.type, y = Percent, fill = Panel)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("EPIC panel" = "#1f77b4", "ICGC DMRs" = "#ff7f0e"), name = "Genomic Region") +
    labs(x = xlab, y = "Percentage (%)") +
    theme_minimal(base_size = 24) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

cpg_labels <- c(hg38_cpg_inter = "CpG_inter", hg38_cpg_islands = "CpG_islands",
                hg38_cpg_shelves = "CpG_shelves", hg38_cpg_shores = "CpG_shores")
gene_labels <- c(hg38_genes_1to5kb = "1to5kb", hg38_genes_3UTRs = "3UTRs", hg38_genes_5UTRs = "5UTRs",
                 hg38_genes_exons = "exons", hg38_genes_introns = "introns", hg38_genes_promoters = "promoters")

ggsave("EFig3c_CpG_annotation.pdf",
       annotation_barplot(annotation_percent("hg38_cpg", cpg_labels), "CpG Annotation"),
       width = 8, height = 7)
ggsave("EFig3d_gene_annotation.pdf",
       annotation_barplot(annotation_percent("hg38_genes", gene_labels), "Gene Annotation"),
       width = 8, height = 7)
