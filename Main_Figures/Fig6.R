## Relevant Fig 6 panels

library(data.table)
library(scater)
library(ggpubr)
library(ggrepel)
library(tidyverse)

cell_levels <- c("buffycoat", "immune", "endo", "fibroblast", "smoothmuscle", "basal", "luminal")
cell_colours <- c("buffycoat" = "#FEDFD4", "immune" = "#CB1B16", "endo" = "#986EAC",
                  "fibroblast" = "#5CAE63", "smoothmuscle" = "#F8A627",
                  "basal" = "#292BC1", "luminal" = "#9DCEE2")


## ---- a, b: bulk tissue ------------------------------------------------------

recommended <- fread("data/recommend_samples_renamed.csv")
icgc_metadata <- fread("data/ICGC_metadata_filter.tsv")

icgc <- fread("data/ICGCtissue_epidish_w_buffycoat_TMC_deconvolve.tsv") %>%
  filter(sample %in% recommended$sampleIDs) %>%
  rowwise() %>%
  mutate(patient = paste0(strsplit(sample, "_")[[1]][1:3], collapse = "_"),
         sample_type = strsplit(sample, "_")[[1]][6]) %>%
  ungroup() %>%
  filter(patient != "CRUK_PC_0050") %>%
  left_join(icgc_metadata, by = c("patient" = "ICGC_ID_short")) %>%
  mutate(gleason_group = case_when(Gleason_biopsy == "3+3" ~ "Gleason 6",
                                   Gleason_biopsy %in% c("3+4", "4+3") ~ "Gleason 7",
                                   TRUE ~ "Gleason 8up"))

tissue_long <- icgc %>%
  filter(sample_type %in% c("BEN", "NRM", "TUM")) %>%
  pivot_longer(fibroblast:buffycoat, names_to = "cell_type", values_to = "fraction") %>%
  mutate(cell_type = factor(cell_type, levels = cell_levels),
         sample_type = factor(sample_type, levels = c("NRM", "BEN", "TUM")))

fig6a <- ggplot(tissue_long, aes(x = sample_type, y = fraction, fill = cell_type)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = cell_colours) +
  labs(title = "ICGC Tissue Deconvolution", x = "Cell Type", y = "Fraction") +
  theme_minimal() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("Fig6a_tissue_composition.pdf", fig6a, width = 7, height = 7)

# tumours only, patients with at most two tumour samples
tumour_long <- icgc %>%
  filter(sample_type == "TUM") %>%
  group_by(patient) %>%
  filter(n() <= 2) %>%
  ungroup() %>%
  pivot_longer(fibroblast:buffycoat, names_to = "cell_type", values_to = "fraction") %>%
  mutate(cell_type = factor(cell_type, levels = cell_levels))

gleason_comparisons <- list(c("Gleason 6", "Gleason 7"), c("Gleason 6", "Gleason 8up"), c("Gleason 7", "Gleason 8up"))

fig6b <- ggplot(tumour_long, aes(x = cell_type, y = fraction, fill = gleason_group)) +
  geom_boxplot() +
  stat_compare_means(comparisons = gleason_comparisons, method = "wilcox.test",
                     label = "p.signif", aes(group = gleason_group)) +
  labs(title = "ICGC Tissue Deconvolution by Gleason Group", x = "Gleason Group", y = "Fraction") +
  theme_minimal() +
  theme(text = element_text(size = 18), legend.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        legend.position = c(0.95, 0.95), legend.justification = c("right", "top"))

ggsave("Fig6b_tumour_composition_by_grade.pdf", fig6b, width = 8, height = 5)


## ---- c: plasma matched to the cell-sorted patients --------------------------

matched_plasma <- fread("data/Atlasplasma_epidish_w_buffycoat_TMC_deconvolve.tsv") %>%
  pivot_longer(fibroblast:buffycoat, names_to = "cell_type", values_to = "fraction")

fig6c <- ggplot(matched_plasma, aes(x = sample, y = fraction, fill = cell_type)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = cell_colours) +
  labs(title = "Atlas Plasma Deconvolution", x = "Sample", y = "Fraction") +
  theme_minimal() +
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("Fig6c_matched_plasma_composition.pdf", fig6c, width = 7, height = 7)


## ---- f, g: plasma cohort ----------------------------------------------------

stage_info <- fread("data/twist_canary_stage_tumfrac.csv")

plasma <- fread("data/canary_epidish_w_buffycoat_TMC_deconvolve.tsv") %>%
  left_join(stage_info, by = c("sample" = "CMDL_ID")) %>%
  mutate(Stage = ifelse(is.na(Stage), "ben", Stage),
         Stage = case_when(Stage %in% c("CPG1", "CPG5") ~ "localised",
                           Stage == "met" ~ "metastatic",
                           Stage == "ben" ~ "benign",
                           TRUE ~ Stage))

plasma_long <- plasma %>%
  pivot_longer(fibroblast:buffycoat, names_to = "cell_type", values_to = "fraction") %>%
  mutate(cell_type = factor(cell_type, levels = cell_levels))

fig6f <- ggplot(plasma_long, aes(x = sample, y = fraction, fill = cell_type)) +
  geom_col() +
  facet_grid(~Stage, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cell_colours) +
  labs(title = "Cell type fractions per sample grouped by Stage", x = "Sample", y = "Cell Fraction", fill = "Cell Type") +
  theme_minimal() +
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("Fig6f_plasma_deconvolution.pdf", fig6f, width = 13, height = 7)

stage_comparisons <- list(c("benign", "localised"), c("benign", "metastatic"), c("localised", "metastatic"))

fig6g <- ggplot(plasma_long, aes(x = Stage, y = fraction, fill = Stage)) +
  geom_boxplot() +
  facet_wrap(~cell_type, scales = "free_y", nrow = 7) +
  stat_compare_means(comparisons = stage_comparisons, method = "wilcox.test",
                     label = "p.signif", aes(group = Stage)) +
  coord_cartesian(clip = "off") +
  labs(title = "Plasma cell fractions", x = "Stage", y = "Cell Fraction") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "none")

ggsave("Fig6g_plasma_fractions_by_stage.pdf", fig6g, width = 5, height = 16)


## ---- d: tumour-vs-benign against immune-vs-buffy coat methylation -----------

n_immune <- 4
n_blood <- 4

quadrant_df <- fread("data/ICGC_cellsort_plasma_BC_betas.tsv") %>%
  mutate(Tumor_diff = diff,
         Immune_diff = Immune_ave_beta - Blood_ave_beta,
         # Welch t-test of immune vs buffy coat from the group means and sds
         Immune_se = sqrt(Immune_beta_sd^2 / n_immune + Blood_beta_sd^2 / n_blood),
         Immune_t = Immune_diff / Immune_se,
         df = (Immune_beta_sd^2 / n_immune + Blood_beta_sd^2 / n_blood)^2 /
           ((Immune_beta_sd^2 / n_immune)^2 / (n_immune - 1) + (Blood_beta_sd^2 / n_blood)^2 / (n_blood - 1)),
         Immune_pval = 2 * pt(-abs(Immune_t), df),
         Quadrant = case_when(Tumor_diff > 0 & Immune_diff > 0 ~ "Q1",
                              Tumor_diff < 0 & Immune_diff < 0 ~ "Q4",
                              TRUE ~ "Other"))

# the ten most significant genes in each concordant quadrant
highlight_df <- quadrant_df %>%
  filter(!is.na(gene), Quadrant %in% c("Q1", "Q4")) %>%
  group_by(Quadrant) %>%
  arrange(Immune_pval) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  group_by(gene) %>%
  slice_min(Immune_pval, n = 1) %>%
  ungroup()

fig6d <- ggplot(quadrant_df, aes(x = Tumor_diff, y = Immune_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(data = filter(quadrant_df, Immune_pval > 0.05), color = "#DEDEDE", alpha = 0.4, size = 3) +
  geom_point(data = filter(quadrant_df, Immune_pval <= 0.05), aes(color = -log10(Immune_pval)),
             alpha = 0.8, size = 3) +
  scale_color_gradient(low = "#377eb9", high = "#e2524f") +
  geom_text_repel(data = highlight_df, aes(label = gene), size = 5, fontface = "bold", color = "black",
                  max.overlaps = 30, box.padding = 0.5, point.padding = 0.5) +
  coord_cartesian(xlim = range(quadrant_df$Tumor_diff) * 1.1, ylim = range(quadrant_df$Immune_diff) * 1.1,
                  expand = FALSE) +
  labs(x = "Tumor vs Benign Beta Diff", y = "Immune vs Buffy Beta Diff", color = "-log10(p-value)") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold"))

ggsave("Fig6d_tumour_immune_quadrants.png", fig6d, width = 6, height = 6, dpi = 600)


## ---- e: top concordant genes in the Tuong et al. single-cell data -----------

immune_quadrant_genes <- c("ZNF343", "NCOR2", "TTC28", "FGD1", "KIF21B", "SPRY4", "SLC35G6", "RNF123",
                           "CDC14A", "TBCD", "LOC105372639", "RIMBP2", "LINC01749", "MGAT2", "AC110614.1",
                           "CCND2", "DPP6")

# genes in the list that are not among the highlighted genes in d
setdiff(immune_quadrant_genes, highlight_df$gene)

sce <- readRDS("data/prostate_portal_300921.RDS")
present <- immune_quadrant_genes[immune_quadrant_genes %in% rownames(sce)]

signature_scores <- tibble(celltype = sce$celltype,
                           group = ifelse(sce$group == "tumor", "tumour", "normal"),
                           score = colMeans(logcounts(sce)[present, , drop = FALSE])) %>%
  filter(celltype != "Sperm")

rm(sce); gc()

fig6e <- ggplot(signature_scores, aes(x = celltype, y = score, fill = group)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.9), outlier.shape = NA) +
  stat_compare_means(aes(group = group), method = "wilcox.test", label = "p.signif",
                     label.y = max(signature_scores$score, na.rm = TRUE), size = 5) +
  labs(title = "Comparison of tumour_immune_genes between Tumour and Normal",
       x = "Cell type", y = "tumour_immune_genes") +
  theme_minimal() +
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Fig6e_immune_genes_tumour_vs_normal.pdf", fig6e, width = 8, height = 4)
