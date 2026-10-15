## Relevant Extended Data Fig 1 panels

library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggpubr)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(tidyverse)


## ---- a, b: mutations and copy number by mevotype ----------------------------

final_metadata <- fread("data/final_ICGC_metadata.tsv")

# columns 1-3 and 60 are sample ids and mevotype, 69-70 and 72-107 SPOP, ETS and mutation calls
mutation_df <- final_metadata %>%
  dplyr::select(c(1:3, 60, 69:70, 72:107)) %>%
  dplyr::filter(colour == "TUM") %>%
  na.omit() %>%
  dplyr::mutate(ETS_num = ifelse(ETS_status.y == "ETS+", 1, 0),
                mevotype = ifelse(mevotype.x == "MEvotype1", "mevotype_1", "mevotype_2")) %>%
  dplyr::select(-ETS_status.y, -mevotype.x) %>%
  dplyr::select(sampleIDs, ICGC_ID_short, mevotype, colour, SPOP, ETS_num, everything())

# columns 108 onwards are the copy number calls
cnv_df <- final_metadata %>%
  dplyr::select(c(1:3, 60, 108:ncol(final_metadata))) %>%
  dplyr::filter(colour == "TUM") %>%
  na.omit() %>%
  dplyr::mutate(mevotype = ifelse(mevotype.x == "MEvotype1", "mevotype_1", "mevotype_2")) %>%
  dplyr::select(-mevotype.x) %>%
  dplyr::select(sampleIDs, ICGC_ID_short, mevotype, colour, everything())

to_binary <- function(df) {
  m <- (df[, 5:ncol(df)] > 0) * 1
  rownames(m) <- df$ICGC_ID_short
  as.data.frame(m)
}
mutation_bin <- to_binary(mutation_df)
cnv_bin <- to_binary(cnv_df)

# Fisher's exact test of each feature against mevotype
fisher_by_mevotype <- function(bin_df, mevotype) {
  bin_df$mevotype <- mevotype
  lapply(setdiff(colnames(bin_df), "mevotype"), function(feature) {
    tab <- table(bin_df[[feature]], bin_df$mevotype)
    if (!all(dim(tab) == c(2, 2))) return(NULL)
    test <- fisher.test(tab)
    tibble(feature = feature, p.value = test$p.value, odds.ratio = test$estimate)
  }) %>% bind_rows()
}

# right-hand annotation: log2 odds ratio and significance
or_annotation <- function(fisher_res) {
  log_or <- log2(fisher_res$odds.ratio)
  finite <- log_or[is.finite(log_or)]
  rowAnnotation(
    OR = anno_simple(log_or,
                     col = colorRamp2(c(min(finite), 0, max(finite)), c("#986EAC", "white", "#5CAE63")),
                     width = unit(6, "mm")),
    pval = anno_text(ifelse(fisher_res$p.value < 0.05, "s", "ns"),
                     gp = gpar(fontsize = 8), width = unit(4, "mm")),
    annotation_name_side = "top")
}

or_legend <- Legend(
  title = "OR",
  col_fun = colorRamp2(c(-1, 0, 1), c("#986EAC", "white", "#5CAE63")),
  at = c(-1, 0, 1),
  labels = c("Enriched in mevotype 1", "Neutral", "Enriched in mevotype 2"))

# a: mutations
fisher_mut <- fisher_by_mevotype(mutation_bin, mutation_df$mevotype)

ht_mut <- Heatmap(
  t(mutation_bin[, fisher_mut$feature]),
  name = "Mutation",
  col = c("0" = "#EDEDED", "1" = "#025196"),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = FALSE,
  column_split = mutation_df$mevotype,
  column_gap = unit(2, "mm"),
  right_annotation = or_annotation(fisher_mut))

pdf("EFig1a_mutation_heatmap.pdf", width = 10, height = 5)
draw(ht_mut, annotation_legend_list = or_legend)
dev.off()

# b: copy number, restricted to the loci with a known cytoband and evotype
cnv_loci <- fread("data/mevotype_evotype_cnv_correl.csv")

cnv_sub <- cnv_bin[, cnv_loci$locus]
colnames(cnv_sub) <- cnv_loci$cytoband[match(colnames(cnv_sub), cnv_loci$locus)]

fisher_cnv <- fisher_by_mevotype(cnv_sub, cnv_df$mevotype)
cnv_plot <- cnv_sub[, fisher_cnv$feature]

evotype_of <- setNames(cnv_loci$evotype, cnv_loci$cytoband)
label_colour <- ifelse(evotype_of[colnames(cnv_plot)] == "alternative", "firebrick",
                       ifelse(evotype_of[colnames(cnv_plot)] == "canonical", "steelblue", "black"))
label_wrapped <- sapply(colnames(cnv_plot), function(x) paste(strwrap(x, width = 22), collapse = "\n"))

ht_cnv <- Heatmap(
  t(cnv_plot),
  name = "CNV",
  col = c("0" = "#EDEDED", "1" = "#025196"),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = FALSE,
  row_labels = label_wrapped,
  row_names_gp = gpar(col = label_colour),
  column_split = cnv_df$mevotype,
  column_gap = unit(2, "mm"),
  right_annotation = or_annotation(fisher_cnv))

evotype_legend <- Legend(
  title = "Evotype",
  labels = c("canonical", "alternative", "ns"),
  type = "text",
  legend_gp = gpar(col = NA),
  labels_gp = gpar(col = c("steelblue", "firebrick", "black"), fontsize = 10))

pdf("EFig1b_cnv_heatmap.pdf", width = 10, height = 5)
draw(ht_cnv, annotation_legend_list = list(or_legend, evotype_legend))
dev.off()


## ---- c-f: purity and Gleason grade ------------------------------------------

meth_hyper <- fread("data/ave_beta_at_TMChyperDMRs.purity.tsv")
meth_hypo  <- fread("data/ave_beta_at_TMChypoDMRs.purity.tsv")

# "3+4" -> 7
gleason_sum <- function(x) sapply(strsplit(x, "\\+"), function(g) sum(as.numeric(g)))
meth_hyper$Gleason_Sum <- gleason_sum(meth_hyper$Gleason_Sum)
meth_hypo$Gleason_Sum  <- gleason_sum(meth_hypo$Gleason_Sum)

purity_scatter <- function(df, direction) {
  ggplot(df, aes(x = purity, y = mean_meth)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(size = 6) +
    labs(x = "Sample purity",
         y = paste("Mean methylation across", direction, "DMRs"),
         title = paste("Methylation vs Purity at", direction, "DMRs")) +
    theme_minimal(base_size = 14) +
    theme(axis.text = element_text(size = 16))
}

ggsave("EFig1c_purity_hypermethylated.pdf", purity_scatter(meth_hyper, "Hypermethylated"), width = 6, height = 6)
ggsave("EFig1d_purity_hypomethylated.pdf", purity_scatter(meth_hypo, "Hypomethylated"), width = 6, height = 6)

grade_comparisons <- list(c("6", "7"), c("6", "8"), c("6", "9"), c("7", "8"), c("7", "9"), c("8", "9"))

gleason_boxplot <- function(df, direction) {
  df %>%
    dplyr::filter(!is.na(Gleason_Sum)) %>%
    ggplot(aes(x = factor(Gleason_Sum), y = mean_meth, fill = factor(Gleason_Sum))) +
    geom_boxplot() +
    stat_compare_means(comparisons = grade_comparisons, method = "wilcox.test") +
    labs(x = "Gleason Grade", y = paste0("Methylation at ", direction, " DMRs (%)")) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none")
}

ggsave("EFig1e_gleason_hypermethylated.pdf", gleason_boxplot(meth_hyper, "Hypermethylated"), width = 5, height = 5)
ggsave("EFig1f_gleason_hypomethylated.pdf", gleason_boxplot(meth_hypo, "Hypomethylated"), width = 5, height = 5)


## ---- g: DMRs that best separate the mevotypes -------------------------------

mevotype_dmrs <- fread("data/top_mevotype_DMRs.tsv")

top_dmrs <- mevotype_dmrs %>%
  dplyr::arrange(dplyr::desc(abs(diff))) %>%
  dplyr::slice_head(n = 50)

# keep the strongest DMR per gene
top_dmrs_unique <- top_dmrs %>%
  dplyr::group_by(gene) %>%
  dplyr::slice_max(abs(diff), n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(abs(diff)))

fig_g <- ggplot(top_dmrs_unique, aes(x = diff, y = reorder(gene, diff))) +
  geom_col(aes(fill = type, alpha = -log10(p_adj))) +
  scale_fill_manual(values = c("hyper" = "#FBA627", "hypo" = "#803E97")) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  labs(x = "More mevotype_1-like (negative) or mevotype_2-like (positive)",
       y = "DMR", title = "Top differential DMRs") +
  theme_classic(base_size = 14)

ggsave("EFig1g_top_mevotype_DMRs.pdf", fig_g, width = 6, height = 12)

# pathway enrichment of these genes (no shared pathway, quoted in the text)
top_entrez <- clusterProfiler::bitr(unique(top_dmrs_unique$gene), fromType = "SYMBOL",
                                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
write.csv(as.data.frame(enrichKEGG(top_entrez, organism = "hsa", pvalueCutoff = 0.05)),
          "EFig1g_top_DMR_genes_KEGG.csv", row.names = FALSE)
write.csv(as.data.frame(enrichPathway(top_entrez, organism = "human", pvalueCutoff = 0.05, readable = TRUE)),
          "EFig1g_top_DMR_genes_Reactome.csv", row.names = FALSE)
