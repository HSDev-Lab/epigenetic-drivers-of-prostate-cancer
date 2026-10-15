## Relevant Extended Data Fig 7 panels

library(data.table)
library(Seurat)
library(scater)
library(ggpubr)
library(tidyverse)

neuronal_genes <- c("NRXN1", "NRXN3", "NTRK3", "SLITRK3", "SLITRK4", "GRIN2D", "GRIA3", "SHANK2",
                    "CAMK2A", "CAMK2B", "CAMKK2", "GABRA5", "GABRB1", "GABRG3", "SLC32A1", "SYT7",
                    "ADCY1", "RASGRF1", "CHAT", "KIF17", "PPFIBP2")
nepc_genes <- c("SYP", "CHGA", "CHGB", "ENO2")

sce <- readRDS("data/prostate_portal_300921.RDS")


## ---- b: neuronal score per cell type, tumour vs normal ----------------------

sce_seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")
sce_seurat <- AddModuleScore(sce_seurat, features = list(neuronal_genes), name = "neuronal_score")

cell_types <- c("LE-KLK3", "LE-KLK4", "Fibroblast", "Endothelial", "MNP", "Mast cell",
                "NK cell", "T cell", "B cell", "BE", "HE", "CE")

patient_scores <- sce_seurat@meta.data %>%
  filter(celltype %in% cell_types) %>%
  group_by(patient, celltype, group) %>%
  summarise(neuronal_score = mean(neuronal_score1, na.rm = TRUE), .groups = "drop") %>%
  mutate(group = factor(group, levels = c("normal", "tumor")),
         celltype = factor(celltype, levels = cell_types))

rm(sce_seurat); gc()

# paired Wilcoxon test per cell type, patients with both samples only
paired_test <- function(ct) {
  wide <- patient_scores %>%
    filter(celltype == ct) %>%
    pivot_wider(names_from = group, values_from = neuronal_score) %>%
    filter(!is.na(normal), !is.na(tumor))
  if (nrow(wide) < 3) return(tibble(celltype = ct, n = nrow(wide), p_value = NA))
  tibble(celltype = ct, n = nrow(wide), p_value = wilcox.test(wide$tumor, wide$normal, paired = TRUE)$p.value)
}

stars <- function(p) {
  case_when(is.na(p) ~ "NA", p < 0.0001 ~ "****", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}

test_labels <- map_dfr(cell_types, paired_test) %>%
  mutate(label = stars(p_value)) %>%
  left_join(patient_scores %>% group_by(celltype) %>%
              summarise(y_pos = max(neuronal_score, na.rm = TRUE) * 1.1) %>%
              mutate(celltype = as.character(celltype)),
            by = "celltype")

fig7b <- ggplot(patient_scores, aes(x = celltype, y = neuronal_score, fill = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              alpha = 0.4, size = 0.8) +
  geom_text(data = test_labels, aes(x = celltype, y = y_pos, label = label), inherit.aes = FALSE, size = 4) +
  scale_fill_manual(values = c("normal" = "#27A9E1", "tumor" = "#F15A2B"), name = "Sample type") +
  labs(x = NULL, y = "Neuronal score",
       title = "Neuronal gene expression, patient-level, by cell type") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("EFig7b_neuronal_score_tumour_vs_normal.pdf", fig7b, width = 13, height = 6)


## ---- c: NEPC score on the UMAP of tumour cells ------------------------------

sce_tumour <- sce[, sce$group == "tumor"]
nepc_present <- nepc_genes[nepc_genes %in% rownames(sce_tumour)]
sce_tumour$NEPC_genes <- colMeans(logcounts(sce_tumour)[nepc_present, , drop = FALSE])

fig7c <- plotReducedDim(sce_tumour, dimred = "UMAP", colour_by = "NEPC_genes") +
  ggtitle("UMAP NEPC genes") +
  theme(plot.title = element_text(size = 18), axis.title = element_text(size = 16),
        axis.text = element_text(size = 16), legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))

ggsave("EFig7c_NEPC_genes_UMAP.pdf", fig7c, width = 5, height = 5)

rm(sce, sce_tumour); gc()


## ---- d: neuronal DMR methylation across disease stages ----------------------

mets <- fread("data/Feng_mCRPC_betas_at_DMRs.tsv") %>% select(-V1)
loc  <- fread("data/meth_summary_perDMR.wide.tsv") %>%
  mutate(DMR = paste0("chr", chr, ".", start, ".", end)) %>%
  select(DMR, everything())

dmr_annot <- fread("data/all_patients_ave_beta_hiQC.summary.tsv") %>%
  mutate(type = if_else(difference < 0, "hypo", "hyper")) %>%
  select(chr, start, end, type, gene)

mets_unique <- mets %>%
  distinct(DMR, .keep_all = TRUE) %>%
  mutate(chr = paste0("chr", chr)) %>%
  select(DMR, chr, start, end, starts_with("beta"))

loc_tum_ben <- loc %>% select(DMR, chr, start, end, matches("TUM$|BEN$"))

merged <- mets_unique %>%
  left_join(loc_tum_ben, by = c("DMR", "chr", "start", "end")) %>%
  left_join(dmr_annot, by = c("chr", "start", "end")) %>%
  rename(DMR_type = type) %>%
  as.data.frame() %>%
  mutate(across(starts_with("beta_"), ~ .x * 100))

# metastatic = beta_ columns, localised = *_TUM, benign = *_BEN
neuronal_stage <- merged %>%
  pivot_longer(cols = -c(DMR, chr, start, end, gene, DMR_type), names_to = "sample", values_to = "beta") %>%
  mutate(sample_type = case_when(grepl("^beta_", sample) ~ "MET",
                                 grepl("_TUM$", sample) ~ "LOC",
                                 grepl("_BEN$", sample) ~ "BEN")) %>%
  filter(gene %in% neuronal_genes, sample_type %in% c("BEN", "LOC", "MET")) %>%
  group_by(DMR, gene, DMR_type, sample_type) %>%
  summarise(mean_beta = mean(beta, na.rm = TRUE), .groups = "drop") %>%
  mutate(sample_type = factor(sample_type, levels = c("BEN", "LOC", "MET")))

fig7d <- ggplot(neuronal_stage, aes(x = sample_type, y = mean_beta, fill = sample_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  facet_wrap(~DMR_type, scales = "free_y") +
  stat_compare_means(comparisons = list(c("BEN", "LOC"), c("LOC", "MET"), c("BEN", "MET")),
                     method = "wilcox.test") +
  scale_fill_manual(values = c("BEN" = "#27A9E1", "LOC" = "#F15A2B", "MET" = "#800000")) +
  labs(x = NULL, y = "Mean methylation (beta)",
       title = "Neuronal DMR methylation: benign vs localised vs metastatic disease") +
  theme_classic(base_size = 14)

ggsave("EFig7d_neuronal_methylation_by_stage.pdf", fig7d, width = 10, height = 6)
