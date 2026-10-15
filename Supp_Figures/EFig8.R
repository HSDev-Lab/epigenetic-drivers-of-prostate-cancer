##Relevant Extended Data Fig 8 panels

library(data.table)
library(readxl)
library(patchwork)
library(gridExtra)
library(tidyverse)

neuronal_genes <- c("NRXN1", "NRXN3", "NTRK3", "SLITRK3", "SLITRK4", "GRIN2D", "GRIA3", "SHANK2",
                    "CAMK2A", "CAMK2B", "CAMKK2", "GABRA5", "GABRB1", "GABRG3", "SLC32A1", "SYT7",
                    "ADCY1", "RASGRF1", "CHAT", "KIF17", "PPFIBP2")

supp_table <- "data/260515_Milne-Clark_Lee Yu_SuppTable.xlsx"


## ---- a: gene-level summary heatmap ------------------------------------------

st1 <- read_excel(supp_table, sheet = "ST1_ICGC_all_DMR")
st3 <- read_excel(supp_table, sheet = "ST3_ICGC_RNAseq_DEGs")

# average methylation of each neuronal gene's DMRs per tumour patient
meth_by_patient <- fread("data/ave_beta_per_DMR_mevo_purity.tsv") %>%
  filter(sample_type == "TUM") %>%
  mutate(chr = str_extract(region, "^[^.]+"),
         start = as.integer(str_extract(region, "(?<=\\.)[0-9]+(?=\\.)")),
         end = as.integer(str_extract(region, "[0-9]+$"))) %>%
  left_join(select(st1, chr, start, end, gene), by = c("chr", "start", "end")) %>%
  filter(gene %in% neuronal_genes) %>%
  group_by(Patient, gene) %>%
  summarise(mean_beta = mean(beta, na.rm = TRUE), .groups = "drop")

# average expression per tumour patient
counts <- fread("data/icgc_counts_keep.txt") %>% column_to_rownames("V1") %>% as.matrix()
counts <- counts[, !duplicated(colnames(counts))]

expr_by_patient <- as.data.frame(counts[rownames(counts) %in% neuronal_genes, , drop = FALSE]) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  mutate(Patient = str_extract(sample, "CRUK_PC_\\d+"),
         sample_type = str_extract(sample, "TUM$|BEN$")) %>%
  filter(sample_type == "TUM") %>%
  group_by(Patient, gene) %>%
  summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")

gene_corr <- inner_join(meth_by_patient, expr_by_patient, by = c("Patient", "gene")) %>%
  group_by(gene) %>%
  summarise(correlation = cor(mean_beta, mean_expr, use = "complete.obs"), .groups = "drop")

# tumour vs benign: methylation difference (ST1) and RNA log2 fold change (ST3)
meth_diff <- st1 %>%
  filter(gene %in% neuronal_genes) %>%
  group_by(gene) %>%
  summarise(meth_difference = mean(difference, na.rm = TRUE), .groups = "drop")

rna_lfc <- st3 %>%
  filter(SYMBOL %in% neuronal_genes) %>%
  select(gene = SYMBOL, rna_log2FC = log2FoldChange)

gene_summary <- meth_diff %>%
  full_join(rna_lfc, by = "gene") %>%
  full_join(gene_corr, by = "gene") %>%
  mutate(hyper_hypo = ifelse(meth_difference > 0, "Hyper", "Hypo"),
         up_down = ifelse(rna_log2FC > 0, "Up", "Down"),
         driver_score = -meth_difference * rna_log2FC)

gene_order <- gene_summary %>%
  arrange(desc(up_down == "Up"), correlation) %>%
  pull(gene)
gene_summary$gene <- factor(gene_summary$gene, levels = rev(gene_order))

tile_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(), plot.margin = margin(2, 2, 2, 2),
        axis.text.x = element_blank(), axis.text.y = element_blank())

tile_column <- function(fill_var, xlab = NULL, labelled = FALSE) {
  p <- ggplot(gene_summary, aes(x = 1, y = gene, fill = .data[[fill_var]])) +
    geom_tile(color = "white", linewidth = 0.6) +
    labs(x = xlab, y = NULL) +
    tile_theme
  if (labelled) p <- p + geom_text(aes(label = round(.data[[fill_var]], 2)), size = 3)
  p
}

p_meth <- tile_column("hyper_hypo") +
  scale_fill_manual(values = c("Hyper" = "#E8A33D", "Hypo" = "#4C6EF5"), na.value = "grey92", name = "Methylation") +
  theme(axis.text.y = element_text(size = 10))

p_rna <- tile_column("up_down") +
  scale_fill_manual(values = c("Up" = "#D55E00", "Down" = "#0072B2"), na.value = "grey92", name = "RNA")

p_corr <- tile_column("correlation", "Correlation", labelled = TRUE) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
                       na.value = "grey92", name = "Meth-Expr\ncorrelation")

p_driver <- tile_column("driver_score", "Driver score", labelled = TRUE) +
  scale_fill_gradient2(low = "#986EAC", mid = "white", high = "#5CAE63", midpoint = 0,
                       na.value = "grey92", name = "Driver\nscore")

fig8a <- (p_meth | p_rna | p_corr | p_driver) +
  plot_layout(widths = rep(0.5, 4), guides = "collect") +
  plot_annotation(title = "Neuronal gene panel: methylation, expression, correlation, and driver status") &
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8))

ggsave("EFig8a_neuronal_gene_panel.pdf", fig8a, width = 11, height = 7)


## ---- b: detection and concordance in the single-cell datasets ---------------

# % of panel genes detected in >5% of cells (genes up or down in bulk RNA-seq), and % of genes whose
# tumour-vs-normal direction in single cells agrees with bulk (paired datasets only)
detection_table <- tribble(
  ~Dataset,                ~`Up genes detected`, ~`Down genes detected`, ~`Up genes concordant`,    ~`Down genes concordant`,
  "Tuong et al. [69]",     "75%",                "27%",                  "76%",                      "52.6%",
  "Apostolov et al. [66]", "88%",                "73%",                  "NA (not paired sample)",   "NA (not paired sample)",
  "Hirz et al. [68]",      "86%",                "50%",                  "92.3%",                    "8.3%",
  "Zaidi et al. [67]",     "75%",                "45%",                  "NA (not paired sample)",   "NA (not paired sample)")

table_grob <- tableGrob(detection_table, rows = NULL,
                        theme = ttheme_minimal(base_size = 11,
                                               core = list(fg_params = list(hjust = 0, x = 0.05)),
                                               colhead = list(fg_params = list(fontface = "bold"))))

pdf("EFig8b_detection_concordance_table.pdf", width = 10, height = 3)
grid::grid.draw(table_grob)
dev.off()
