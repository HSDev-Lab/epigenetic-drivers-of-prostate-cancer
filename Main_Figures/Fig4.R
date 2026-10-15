## Relevant Fig 4 panels

library(data.table)
library(Seurat)
library(ggpubr)
library(org.Hs.eg.db)
library(tidyverse)

neuronal_genes <- c("NRXN1", "NRXN3", "NTRK3", "SLITRK3", "SLITRK4", "GRIN2D", "GRIA3", "SHANK2",
                    "CAMK2A", "CAMK2B", "CAMKK2", "GABRA5", "GABRB1", "GABRG3", "SLC32A1", "SYT7",
                    "ADCY1", "RASGRF1", "CHAT", "KIF17", "PPFIBP2")
nepc_genes <- c("SYP", "CHGA", "CHGB", "ENO2")


## ---- b: Apostolov et al. ----------------------------------------------------

apo <- readRDS("data/Apostolov_et_al_2024_h5ad_Seurat_converted_TMC_160426.rds")

# genes are Ensembl ids; rename to symbols, keeping the Ensembl id when there is no symbol
expr <- GetAssayData(apo, assay = "RNA", layer = "data")
symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(expr), keytype = "ENSEMBL",
                                 column = "SYMBOL", multiVals = "first")
new_names <- ifelse(is.na(symbols), rownames(expr), symbols)
keep <- !duplicated(new_names)
expr <- expr[keep, ]
rownames(expr) <- new_names[keep]

apo_symbols <- CreateSeuratObject(counts = expr)
apo_symbols <- SetAssayData(apo_symbols, layer = "data", new.data = expr)
apo_symbols$cell_type <- apo$cell_type
apo_symbols$sample_id <- apo$sample_id
apo_symbols <- AddModuleScore(apo_symbols, features = list(neuronal_genes), name = "neuronal_score")

apo_scores <- apo_symbols@meta.data %>%
  group_by(sample_id, cell_type) %>%
  summarise(neuronal_score = mean(neuronal_score1, na.rm = TRUE), .groups = "drop")

rm(apo, apo_symbols, expr); gc()

cell_type_order <- apo_scores %>%
  group_by(cell_type) %>%
  summarise(med = median(neuronal_score, na.rm = TRUE)) %>%
  arrange(desc(med)) %>%
  pull(cell_type)

neuronal_high <- c("glial cell", "prostate neuroendocrine cell",
                   "luminal cell of prostate epithelium", "ciliated epithelial cell")

apo_scores <- apo_scores %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_order),
         highlight = ifelse(cell_type %in% neuronal_high, "Neuronal-high", "Other"))

fig4b <- ggplot(apo_scores, aes(x = cell_type, y = neuronal_score, fill = highlight)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, linewidth = 0.3) +
  scale_x_discrete(labels = function(x) str_wrap(x, 12)) +
  scale_fill_manual(values = c("Neuronal-high" = "#9DCEE2", "Other" = "grey80")) +
  labs(x = NULL, y = "Neuronal score", title = "Neuronal gene expression by cell type") +
  theme_bw(base_size = 5.5) +
  theme(axis.text.x = element_text(size = 4.2, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 5),
        plot.title = element_text(size = 6.5),
        panel.grid = element_blank(), legend.position = "none")

ggsave("Fig4b_Apostolov_neuronal_score.pdf", fig4b, width = 205, height = 60, units = "mm")


## ---- c: Hirz et al. ---------------------------------------------------------

hirz <- readRDS("data/Hirz_scRNA_annotated_HLY.rds")

hirz_scores <- hirz@meta.data %>%
  group_by(sample, cell1, type) %>%
  summarise(neuronal_score = mean(neuronal1, na.rm = TRUE), .groups = "drop") %>%
  mutate(type = factor(type, levels = c("normal", "adj_ben", "tum")))

rm(hirz); gc()

sample_type_colours <- c("normal" = "#A8D8F0", "adj_ben" = "#27A9E1", "tum" = "#F15A2B")

fig4c <- ggplot(hirz_scores, aes(x = cell1, y = neuronal_score, fill = type)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.4, size = 0.4) +
  scale_fill_manual(values = sample_type_colours,
                    labels = c("Normal", "Adjacent Benign", "Tumour"), name = "Sample type") +
  labs(x = NULL, y = "Neuronal score", title = "Neuronal gene expression by cell type and sample type") +
  theme_bw(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
        axis.text.y = element_text(size = 5),
        plot.title = element_text(size = 6.5),
        legend.text = element_text(size = 5), legend.title = element_text(size = 5.5),
        legend.key.size = unit(0.3, "cm"))

ggsave("Fig4c_Hirz_neuronal_score.pdf", fig4c, width = 205, height = 60, units = "mm")


## ---- d, e: Zaidi et al. -----------------------------------------------------

zaidi <- readRDS("data/Zaidi_PBMC_combi_scRNA_annotated_HLY.rds")
zaidi <- subset(zaidi, subset = !is.na(coarse_ano))
zaidi <- AddModuleScore(zaidi, features = list(neuronal_genes), name = "neuronal_score")
zaidi <- AddModuleScore(zaidi, features = list(nepc_genes), name = "nepc_score")

zaidi_scores <- zaidi@meta.data %>%
  group_by(patient, coarse_ano) %>%
  summarise(neuronal_score = mean(neuronal_score1, na.rm = TRUE),
            nepc_score = mean(nepc_score1, na.rm = TRUE),
            subtype = first(subtype), .groups = "drop") %>%
  mutate(subtype = factor(subtype, levels = c("CSPC", "CRPC", "NEPC")))

rm(zaidi); gc()

zaidi_boxplot <- function(score, title, ylab) {
  ggplot(zaidi_scores, aes(x = subtype, y = .data[[score]], fill = subtype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.6) +
    stat_compare_means(comparisons = list(c("CSPC", "CRPC"), c("CRPC", "NEPC"), c("CSPC", "NEPC")),
                       method = "wilcox.test", size = 2, tip.length = 0.01, bracket.size = 0.2) +
    facet_wrap(~coarse_ano, nrow = 1, scales = "free_y") +
    labs(title = title, x = "Subtype", y = ylab) +
    theme_bw(base_size = 6) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          strip.text = element_text(size = 5.5),
          panel.spacing = unit(2, "pt"),
          plot.title = element_text(size = 6.5))
}

ggsave("Fig4d_Zaidi_neuronal_score.pdf",
       zaidi_boxplot("neuronal_score", "Neuronal gene expression by cell compartment and subtype", "Neuronal score"),
       width = 85, height = 85, units = "mm")
ggsave("Fig4e_Zaidi_NEPC_score.pdf",
       zaidi_boxplot("nepc_score", "NEPC gene expression by cell compartment and subtype", "NEPC score"),
       width = 85, height = 85, units = "mm")
