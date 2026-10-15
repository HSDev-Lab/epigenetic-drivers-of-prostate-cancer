## Relevant Extended Data Fig 5 panels

library(data.table)
library(ComplexHeatmap)
library(viridis)
library(tidyverse)

reference_colours <- c("Blood-T" = "#FEDFD4", "Blood-B" = "#65010C", "Blood-Mono+Macro" = "#CB1B16",
                       "Blood-Granul" = "#F26A4F", "Blood-NK" = "#F8A592", "Endothelium" = "#986EAC",
                       "Heart-Fibro" = "#5CAE63", "Smooth-Musc" = "#F8A627", "Prostate-Ep" = "#336699")


## ---- b: reference DMRs in the sorted samples --------------------------------

ref_betas <- fread("data/all_patients_ave_beta_kaplanDMR.tsv") %>%
  mutate(region = paste0(chr, ".", start, ".", end)) %>%
  distinct(region, .keep_all = TRUE) %>%
  select(region, target, ends_with("beta")) %>%
  na.omit() %>%
  rename(cell_type = target)

ref_mtx <- ref_betas %>%
  select(-cell_type) %>%
  column_to_rownames("region")
colnames(ref_mtx) <- gsub("_beta$", "", colnames(ref_mtx))

row_split <- factor(ref_betas$cell_type)

# tumour fraction and ploidy, in the same order as the heatmap columns
ichor <- fread("data/all_patients_ichor.txt")
colnames(ichor) <- c("sample", "TF", "ploidy")
ichor <- ichor %>%
  mutate(sample = gsub(".params.*", "", sample),
         sample = gsub(".dedup.sorted_ichor", "", sample),
         sample = gsub("_1_val_1_bismark_bt2_pe.deduplicated_sorted_ichor", "", sample)) %>%
  filter(sample %in% colnames(ref_mtx)) %>%
  mutate(TF = TF * 100) %>%
  arrange(match(sample, colnames(ref_mtx)))

row_ha <- rowAnnotation(cell_type = row_split,
                        col = list(cell_type = reference_colours),
                        width = unit(0.8, "cm"))

top_ha <- HeatmapAnnotation(
  TF = anno_barplot(ichor$TF, border = FALSE, gp = gpar(fill = "firebrick"), axis = TRUE),
  .spacer = anno_empty(border = FALSE, height = unit(1, "mm")),
  Ploidy = anno_points(ichor$ploidy, ylim = c(0, 4), pch = 16, size = unit(2, "mm"),
                       gp = gpar(col = "black"), baseline = 2,
                       baseline_gp = gpar(lty = 2, col = "black")),
  annotation_name_side = "left")

ht <- Heatmap(as.matrix(ref_mtx),
              name = "Beta",
              col = plasma(100),
              cluster_rows = TRUE, cluster_columns = TRUE,
              show_row_names = FALSE, show_column_names = TRUE,
              top_annotation = top_ha,
              left_annotation = row_ha,
              row_split = row_split,
              row_title = "Cell-specific hypomethylated DMRs")

pdf("EFig5b_sorted_samples_reference_DMRs.pdf", width = 9, height = 6)
draw(ht)
dev.off()


## ---- c: deconvolution of the sorted samples ---------------------------------

fractions <- fread("data/cell_sort_prostate_groups_deconv2.csv") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "Fraction") %>%
  mutate(Sample = gsub("\\.dedup\\.sorted.*", "", Sample),
         CellType = factor(CellType, levels = names(reference_colours))) %>%
  filter(Sample %in% ichor$sample)

fig5c <- ggplot(fractions, aes(x = Sample, y = Fraction, fill = CellType)) +
  geom_col() +
  scale_fill_manual(values = reference_colours) +
  labs(x = "Sample", y = "Fraction", fill = "Cell Type") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("EFig5c_sorted_samples_deconvolution.pdf", fig5c, width = 12, height = 6)
