## #Relevant Fig2 panels

library(data.table)
library(tidyverse)
library(ggrepel)
library(annotatr)

dmr <- fread("data/ICGC_DMRs.tsv")
deg <- fread("data/ICGC_RNAseq_DEGs.csv")

# one DMR and one transcript per gene: the largest effect
dmr_unique <- dmr %>%
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  slice_max(abs(difference), n = 1, with_ties = FALSE) %>%
  ungroup()

deg_unique <- deg %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  slice_max(abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

dmr_deg <- inner_join(dmr_unique, deg_unique, by = c("gene" = "SYMBOL")) %>%
  filter(abs(difference) > 0.2 & abs(log2FoldChange) > 0.5) %>%
  mutate(difference = difference * 100)

# promoter or not
annotations <- build_annotations(
  genome = "hg38",
  annotations = c("hg38_genes_promoters", "hg38_genes_5UTRs", "hg38_genes_exons",
                  "hg38_genes_introns", "hg38_genes_3UTRs"))

dmr_deg_gr <- makeGRangesFromDataFrame(dmr_deg, seqnames.field = "chr",
                                       start.field = "start", end.field = "end",
                                       keep.extra.columns = TRUE)

drivers <- annotate_regions(dmr_deg_gr, annotations = annotations,
                            ignore.strand = TRUE, quiet = TRUE) %>%
  as.data.frame() %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(driver_score = -difference * log2FoldChange,
         promoter_status = ifelse(annot.type == "hg38_genes_promoters", "Promoter", "Non-promoter"),
         driver = ifelse(driver_score < 0, "Non-Canonical", "Canonical"))

write.table(drivers, "DMR_DEG_genes_annotated.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

label_genes <- bind_rows(top_n(drivers, 10, driver_score),
                         top_n(drivers, -10, driver_score))

promoter_cols <- c("Promoter" = "#986EAC", "Non-promoter" = "#5CAE63")


## ---- c: methylation difference vs expression --------------------------------

quadrant_labels <- data.frame(
  difference = c(25, -25, 25, -25),
  log2FoldChange = c(-2.8, 4.4, 4.4, -2.8),
  label = c("Repression", "Activation", "Non Canonical", "Non Canonical"))

fig2c <- ggplot(drivers, aes(x = difference, y = log2FoldChange)) +
  geom_point(aes(color = promoter_status, size = -log10(pvalue)), alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(data = quadrant_labels, aes(label = label), size = 6, fontface = "italic") +
  geom_text_repel(data = label_genes, aes(label = annot.symbol), size = 4.5, max.overlaps = 20) +
  scale_color_manual(values = promoter_cols) +
  labs(x = "Average Methylation Difference", y = "RNA log2 Fold Change") +
  theme_minimal(base_size = 16)
fig2c

ggsave("Fig2c_methylation_vs_expression.pdf", fig2c, width = 8, height = 8)


## ---- d: genes with the most extreme driver scores ---------------------------

fig2d <- label_genes %>%
  filter(!is.na(annot.symbol)) %>%
  ggplot(aes(x = reorder(annot.symbol, driver_score), y = driver_score, fill = promoter_status)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = promoter_cols) +
  labs(x = "Gene Symbol", y = "Driver Score") +
  theme_minimal(base_size = 16)

ggsave("Fig2d_driver_scores.pdf", fig2d, width = 8, height = 8)


## ---- numbers quoted in the text ---------------------------------------------

# canonical = hypermethylated and down, or hypomethylated and up
table(drivers$driver)

# methylation direction vs expression direction
meth_vs_rna <- table(RNA = ifelse(drivers$log2FoldChange > 0, "upreg", "downreg"),
                     methylation = ifelse(drivers$difference > 0, "hyper", "hypo"))
meth_vs_rna <- meth_vs_rna[c("upreg", "downreg"), c("hypo", "hyper")]
meth_vs_rna
fisher.test(meth_vs_rna)

# canonical vs non-canonical, promoter vs non-promoter
promoter_vs_driver <- table(drivers$promoter_status, drivers$driver)[c("Promoter", "Non-promoter"), ]
promoter_vs_driver
fisher.test(promoter_vs_driver)
