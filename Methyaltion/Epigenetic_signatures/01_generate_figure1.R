#!/usr/bin/env Rscript

# Figure 1: paired tumour-benign methylation by MEvotype.
# A: 717 MEvotype-specific/divergent DMRs.
# B: prespecified neuronal DMRs.
# C: prespecified RhoGTPase DMRs.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplotify)
  library(patchwork)
  library(grid)
})

set.seed(20260904)
grDevices::pdf(file = NULL)
on.exit(grDevices::dev.off(), add = TRUE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}
source(file.path(script_dir, "shared_data.R"))
package_dir <- package_path(script_dir)
input_dir <- file.path(package_dir, "inputs")
figure_dir <- file.path(package_dir, "figures")
table_dir <- file.path(package_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

data_root <- Sys.getenv(
  "MEVOTYPE_DATA_ROOT", unset = "/Users/tobymilne-clark/Documents/Data"
)
dmr_file <- file.path(
  data_root, "DMR_calling", "ICGC", "DMRs_with_coords_350_samples.csv"
)
stopifnot(file.exists(dmr_file))

col_m1 <- MEVOTYPE_COLOURS[["MEvotype1"]]
col_m2 <- MEVOTYPE_COLOURS[["MEvotype2"]]
col_hyper <- "#F8A627"
col_hypo <- "#803E97"
col_present <- "#025196"
col_absent <- "#EDEDED"
col_na <- "#B2B2B2"
col_tumour <- "#F15A2B"

# Paired tumour-benign methylation matrix -----------------------------------
cohort <- read_csv(
  file.path(input_dir, "current_patient_analysis_any_relapse.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    mevotype = factor(mevotype, levels = MEVOTYPE_LEVELS),
    high_grade = as.integer(grade_group %in% c("GG4", "GG5"))
  )
stopifnot(nrow(cohort) == 163L, !anyDuplicated(cohort$patient))

purity <- read_csv(
  file.path(input_dir, "ICGC_sample_annotation.csv"), show_col_types = FALSE
) %>%
  transmute(patient = ICGC_sample_ID, purity) %>%
  filter(!is.na(patient)) %>%
  distinct(patient, .keep_all = TRUE)

dmr_raw <- fread(dmr_file)
dmr_ids <- paste(dmr_raw$chr, dmr_raw$start, dmr_raw$end, sep = ":")
stopifnot(nrow(dmr_raw) == 8915L, !anyDuplicated(dmr_ids))
benign_samples <- names(dmr_raw)[str_ends(names(dmr_raw), "_BEN")]
benign_map <- tibble(
  benign_sample = benign_samples,
  patient = str_extract(benign_samples, "CRUK_PC_[0-9]{4}")
) %>% distinct(patient, .keep_all = TRUE)

pairs <- cohort %>%
  select(patient, methylation_sample, mevotype, grade_group, high_grade,
         advanced_stage) %>%
  inner_join(benign_map, by = "patient") %>%
  left_join(purity, by = "patient")
stopifnot(
  nrow(pairs) == 163L,
  !anyDuplicated(pairs$patient),
  all(pairs$methylation_sample %in% names(dmr_raw)),
  all(pairs$benign_sample %in% names(dmr_raw))
)

tumour_columns <- pairs$methylation_sample
benign_columns <- pairs$benign_sample
tumour_matrix <- as.matrix(dmr_raw[, ..tumour_columns])
benign_matrix <- as.matrix(dmr_raw[, ..benign_columns])
storage.mode(tumour_matrix) <- "double"
storage.mode(benign_matrix) <- "double"
paired_delta <- tumour_matrix - benign_matrix
rownames(paired_delta) <- dmr_ids
colnames(paired_delta) <- pairs$patient

# Panel A -------------------------------------------------------------------
specific_changes <- read_csv(
  file.path(input_dir, "MEvotype_specific_change_DMRs_717.csv"),
  show_col_types = FALSE
)
stopifnot(
  nrow(specific_changes) == 717L,
  !anyDuplicated(specific_changes$dmr_id),
  all(specific_changes$dmr_id %in% rownames(paired_delta))
)
heat_a <- paired_delta[specific_changes$dmr_id, , drop = FALSE]
column_order <- order(pairs$mevotype, replace_na(pairs$purity, -Inf))
heat_a <- heat_a[, column_order, drop = FALSE]
meta_a <- pairs[column_order, ]

event_levels <- c(
  "MEvotype1-specific hypermethylation",
  "MEvotype1-specific hypomethylation",
  "MEvotype2-specific hypermethylation",
  "MEvotype2-specific hypomethylation",
  "Divergent/opposite change"
)
row_event <- factor(specific_changes$event_class, levels = event_levels)
top_a <- HeatmapAnnotation(
  MEvotype = meta_a$mevotype,
  Purity = meta_a$purity,
  `High grade` = factor(meta_a$high_grade, c(0, 1), c("No", "Yes")),
  col = list(
    MEvotype = MEVOTYPE_COLOURS,
    Purity = colorRamp2(c(0.2, 0.5, 0.8), c("#E7F4F5", "#6BC3A5", "#198F55")),
    `High grade` = c("No" = col_absent, "Yes" = col_tumour)
  ),
  annotation_name_gp = gpar(fontsize = 6.5, fontface = "bold"),
  annotation_legend_param = list(
    title_gp = gpar(fontface = "bold", fontsize = 7.5),
    labels_gp = gpar(fontsize = 6.5)
  )
)

ht_a <- Heatmap(
  heat_a,
  name = "Delta beta",
  col = colorRamp2(c(-0.45, 0, 0.45), c(col_hypo, "#FFFFFF", col_hyper)),
  top_annotation = top_a,
  row_split = row_event,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = meta_a$mevotype,
  column_gap = unit(1.0, "mm"),
  row_gap = unit(0.7, "mm"),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 6.0, fontface = "bold"),
  column_title_gp = gpar(fontsize = 7.5, fontface = "bold"),
  heatmap_legend_param = list(
    title = expression(Delta * beta), at = c(-0.4, 0, 0.4),
    title_gp = gpar(fontface = "bold", fontsize = 7.5),
    labels_gp = gpar(fontsize = 6.5)
  ),
  column_title = "MEvotype-specific and divergent tumour changes"
)
p_a <- as.ggplot(grid.grabExpr(
  draw(ht_a, heatmap_legend_side = "right", annotation_legend_side = "right"),
  wrap = TRUE
))

# Panels B-C ----------------------------------------------------------------
programme_statistics <- list()
make_programme_heatmap <- function(tab, programme_name, title, row_width_mm = 28) {
  tab <- tab %>%
    mutate(
      programme_order = factor(gene_order),
      row_label = if_else(
        dmr_count_gene > 1,
        paste0(gene_names, " [", dmr_within_gene, "]"),
        gene_names
      )
    ) %>%
    arrange(programme_order, dmr_within_gene)
  stopifnot(all(tab$dmr_id %in% rownames(paired_delta)))
  delta <- paired_delta[tab$dmr_id, , drop = FALSE]
  m2_idx <- pairs$mevotype == "MEvotype2"
  m1_idx <- pairs$mevotype == "MEvotype1"
  group_means <- cbind(
    MEvotype2 = rowMeans(delta[, m2_idx, drop = FALSE], na.rm = TRUE),
    MEvotype1 = rowMeans(delta[, m1_idx, drop = FALSE], na.rm = TRUE)
  )
  p_values <- vapply(seq_len(nrow(delta)), function(i) {
    suppressWarnings(wilcox.test(
      delta[i, m1_idx], delta[i, m2_idx], exact = FALSE
    )$p.value)
  }, numeric(1))
  fdr <- p.adjust(p_values, method = "BH")
  difference <- group_means[, "MEvotype1"] - group_means[, "MEvotype2"]
  direction <- ifelse(difference >= 0, "Higher in M1", "Higher in M2")
  programme_statistics[[programme_name]] <<- tibble(
    programme = programme_name,
    dmr_id = tab$dmr_id,
    gene_names = tab$gene_names,
    MEvotype2_mean_delta_beta = group_means[, "MEvotype2"],
    MEvotype1_mean_delta_beta = group_means[, "MEvotype1"],
    M1_minus_M2 = difference,
    wilcoxon_p = p_values,
    BH_FDR = fdr,
    direction = direction
  )

  rownames(group_means) <- make.unique(tab$row_label)
  fdr_score <- pmin(-log10(pmax(fdr, 1e-12)), 12)
  rna_limit <- max(abs(tab$RNA_logFC), na.rm = TRUE)
  if (!is.finite(rna_limit) || rna_limit == 0) rna_limit <- 1
  right_ann <- rowAnnotation(
    Direction = direction,
    `M1-M2 -log10 FDR` = fdr_score,
    `RNA logFC` = tab$RNA_logFC,
    col = list(
      Direction = c("Higher in M1" = col_hyper, "Higher in M2" = col_hypo),
      `M1-M2 -log10 FDR` = colorRamp2(
        c(0, 1.3, 6), c("#F2F2F2", "#A6CEE3", col_present)
      ),
      `RNA logFC` = colorRamp2(
        c(-rna_limit, 0, rna_limit), c(col_m2, "#FFFFFF", col_m1)
      )
    ),
    na_col = col_na,
    annotation_name_gp = gpar(fontsize = 6.2, fontface = "bold"),
    annotation_legend_param = list(
      title_gp = gpar(fontface = "bold", fontsize = 7),
      labels_gp = gpar(fontsize = 6.2)
    )
  )
  ht <- Heatmap(
    group_means,
    name = "Delta beta",
    col = colorRamp2(c(-0.45, 0, 0.45), c(col_hypo, "#FFFFFF", col_hyper)),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6.0),
    row_names_max_width = unit(row_width_mm, "mm"),
    column_names_gp = gpar(fontsize = 7.0, fontface = "bold"),
    column_names_rot = 25,
    right_annotation = right_ann,
    column_title = title,
    column_title_gp = gpar(fontsize = 8, fontface = "bold"),
    heatmap_legend_param = list(
      title = expression(Delta * beta), at = c(-0.4, 0, 0.4),
      title_gp = gpar(fontface = "bold", fontsize = 7),
      labels_gp = gpar(fontsize = 6.2)
    )
  )
  as.ggplot(grid.grabExpr(
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right"),
    wrap = TRUE
  ))
}

neuronal_tab <- read_csv(
  file.path(input_dir, "Neuronal_21_gene_DMR_source.csv"), show_col_types = FALSE
)
rho_tab <- read_csv(
  file.path(input_dir, "RhoGTPase_DMR_source.csv"), show_col_types = FALSE
)
stopifnot(nrow(neuronal_tab) == 30L, nrow(rho_tab) == 27L)
p_b <- make_programme_heatmap(
  neuronal_tab, "Neuronal", "Prespecified neuronal DMRs by MEvotype", 28
)
p_c <- make_programme_heatmap(
  rho_tab, "RhoGTPase", "Prespecified RhoGTPase DMRs by MEvotype", 30
)

write_csv(specific_changes, file.path(table_dir, "Figure_1_MEvotype_specific_DMRs.csv"))
write_csv(
  meta_a %>% select(patient, mevotype, purity, high_grade),
  file.path(table_dir, "Figure_1_heatmap_sample_annotation.csv")
)
write_csv(
  bind_rows(programme_statistics),
  file.path(table_dir, "Figure_1_programme_MEvotype_comparison.csv")
)
write_csv(
  tibble(
    check = c(
      "paired patients", "MEvotype1 patients", "MEvotype2 patients",
      "recurrence events", "source DMRs", "MEvotype-specific DMRs",
      "neuronal programme DMRs", "RhoGTPase programme DMRs",
      "missing purity values"
    ),
    value = c(
      nrow(pairs), sum(pairs$mevotype == "MEvotype1"),
      sum(pairs$mevotype == "MEvotype2"), sum(cohort$rfs_event),
      nrow(dmr_raw), nrow(specific_changes), nrow(neuronal_tab), nrow(rho_tab),
      sum(is.na(meta_a$purity))
    )
  ),
  file.path(table_dir, "Figure_1_input_audit.csv")
)

# Assemble A-C only ----------------------------------------------------------
layout_design <- c(
  area(t = 1, l = 1, b = 2, r = 1),
  area(t = 1, l = 2, b = 1, r = 2),
  area(t = 2, l = 2, b = 2, r = 2)
)
figure_1 <- p_a + p_b + p_c +
  plot_layout(
    design = layout_design,
    widths = c(1.35, 1),
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14, family = "Arial")
    )
  )

ggsave(
  file.path(figure_dir, "Figure_1_MEvotype_DMR_programmes.pdf"),
  figure_1, width = 16, height = 13, device = cairo_pdf,
  family = "Arial", bg = "white"
)
ggsave(
  file.path(figure_dir, "Figure_1_MEvotype_DMR_programmes.png"),
  figure_1, width = 16, height = 13, dpi = 400, bg = "white"
)
message("Generated Figure 1 (A-C, MEvotype1 versus MEvotype2)")
