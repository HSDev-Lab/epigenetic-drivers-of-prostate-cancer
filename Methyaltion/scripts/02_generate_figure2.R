#!/usr/bin/env Rscript

# Figure 2
# A: neuronal and RhoGTPase scores by MEvotype1 versus MEvotype2.
# B: fixed programme-high composition by MEvotype.
# C: recurrence-free survival by fixed programme burden in the full cohort.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(survival)
})

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

dat <- load_programme_data(file.path(
  input_dir, "current_patient_analysis_any_relapse.csv"
))
colours <- MEVOTYPE_COLOURS
col_high <- "#C44E52"
col_text <- "#222222"

theme_paper <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 11, colour = col_text),
    plot.subtitle = element_text(size = 8.7, colour = col_text),
    plot.caption = element_text(size = 7.6, colour = "#555555", hjust = 0),
    axis.title = element_text(face = "bold", colour = col_text),
    axis.text = element_text(colour = col_text),
    strip.background = element_rect(fill = "#F1F1F1", colour = NA),
    strip.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Panel A -------------------------------------------------------------------
score_long <- rbindlist(list(
  dat[, .(patient, mevotype, programme = "Neuronal", score = neuronal_score)],
  dat[, .(patient, mevotype, programme = "RhoGTPase", score = rho_score)]
))
score_long[, programme := factor(programme, levels = c("Neuronal", "RhoGTPase"))]
cutoff_lines <- data.table(
  programme = factor(c("Neuronal", "RhoGTPase"), levels = levels(score_long$programme)),
  cutoff = c(NEURONAL_CUTOFF, RHOGTPASE_CUTOFF)
)
comparison_labels <- score_long[, {
  test <- suppressWarnings(wilcox.test(
    score[mevotype == "MEvotype1"], score[mevotype == "MEvotype2"],
    exact = FALSE, alternative = "two.sided"
  ))
  .(
    W = unname(test$statistic),
    p_value = test$p.value,
    label = paste0("Wilcoxon p ", format.pval(test$p.value, digits = 2, eps = 1e-16))
  )
}, by = programme]
comparison_labels[, x := 1.5]

p_a <- ggplot(score_long, aes(mevotype, score, colour = mevotype)) +
  geom_hline(
    data = cutoff_lines, aes(yintercept = cutoff),
    inherit.aes = FALSE, linetype = 2, colour = "#333333", linewidth = 0.65
  ) +
  geom_boxplot(width = 0.56, outlier.shape = NA, alpha = 0.14, linewidth = 0.5) +
  geom_jitter(width = 0.12, height = 0, size = 1.15, alpha = 0.62) +
  geom_text(
    data = comparison_labels, aes(x = x, y = Inf, label = label),
    inherit.aes = FALSE, vjust = 1.35, size = 3.0, colour = col_text
  ) +
  facet_wrap(~programme, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = colours, guide = "none") +
  labs(
    title = "A  Methylation programme scores by MEvotype",
    subtitle = "Full paired cohort; dashed cut-offs: neuronal >=0.8 and RhoGTPase >=1.0",
    x = NULL, y = "Methylation programme score"
  ) +
  theme_paper

# Panel B -------------------------------------------------------------------
programme_counts <- dat[, .N, by = .(mevotype, high_programme_count)]
programme_counts[, fraction := N / sum(N), by = mevotype]
p_b <- ggplot(
  programme_counts,
  aes(mevotype, fraction, fill = factor(high_programme_count))
) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(N >= 4, N, "")),
    position = position_stack(vjust = 0.5), size = 3.0, colour = col_text
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c("0" = "#D9E8F5", "1" = "#9CC2E5", "2" = col_high),
    name = "High programmes"
  ) +
  labs(
    title = "B  Programme-high composition by MEvotype",
    subtitle = "Number of fixed high programmes per tumour",
    x = NULL, y = "Patients"
  ) +
  theme_paper

# Panel C -------------------------------------------------------------------
count_colours <- c(
  "0 high programmes" = "#4C78A8",
  "1 high programme" = "#D39C32",
  "2 high programmes" = col_high
)
dat[, programme_burden := factor(
  high_programme_count, levels = 0:2, labels = names(count_colours)
)]
fit_count_ordinal <- coxph(
  Surv(rfs_years, rfs_event) ~ high_programme_count, data = dat, x = TRUE
)
fit_count_factor <- coxph(
  Surv(rfs_years, rfs_event) ~ relevel(programme_burden, ref = "0 high programmes"),
  data = dat, x = TRUE
)
count_sm <- summary(fit_count_ordinal)
count_hr <- unname(count_sm$conf.int[1, "exp(coef)"])
count_low <- unname(count_sm$conf.int[1, "lower .95"])
count_high <- unname(count_sm$conf.int[1, "upper .95"])
count_p <- unname(count_sm$coefficients[1, "Pr(>|z|)"])

fit_count_km <- survfit(Surv(rfs_years, rfs_event) ~ programme_burden, data = dat)
count_km_sm <- summary(fit_count_km, censored = TRUE)
count_curve <- data.table(
  time = count_km_sm$time,
  survival = count_km_sm$surv,
  n_censor = count_km_sm$n.censor,
  programme_burden = sub("^programme_burden=", "", as.character(count_km_sm$strata))
)
count_starts <- data.table(
  time = 0, survival = 1, n_censor = 0,
  programme_burden = sub("^programme_burden=", "", names(fit_count_km$strata))
)
count_curve <- rbind(count_starts, count_curve, fill = TRUE)
count_curve[, programme_burden := factor(programme_burden, levels = names(count_colours))]

count_summary <- dat[, .(n = .N, events = sum(rfs_event)), by = programme_burden]
setorder(count_summary, programme_burden)
count_summary[, legend_label := sprintf(
  "%s (n=%d; events=%d)", programme_burden, n, events
)]
legend_labels <- setNames(count_summary$legend_label, as.character(count_summary$programme_burden))
risk_times <- seq(0, 8, by = 2)
count_risk_sm <- summary(fit_count_km, times = risk_times, extend = TRUE)
count_risk <- data.table(
  time = count_risk_sm$time,
  n_risk = count_risk_sm$n.risk,
  programme_burden = sub("^programme_burden=", "", as.character(count_risk_sm$strata))
)
count_risk[, programme_burden := factor(
  programme_burden, levels = rev(names(count_colours))
)]

p_c_curve <- ggplot(
  count_curve,
  aes(time, survival, colour = programme_burden, group = programme_burden)
) +
  geom_step(linewidth = 0.78) +
  geom_point(
    data = count_curve[n_censor > 0], shape = 3, size = 1.15, stroke = 0.5
  ) +
  scale_colour_manual(
    values = count_colours, labels = legend_labels,
    breaks = names(count_colours), drop = FALSE
  ) +
  scale_x_continuous(breaks = risk_times, expand = expansion(mult = c(0.02, 0.01))) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.25),
    labels = percent_format(accuracy = 1), expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(0, 8), ylim = c(0, 1)) +
  labs(
    title = "C  Recurrence-free survival by programme burden",
    subtitle = sprintf(
      "Ordinal HR %.2f per programme (95%% CI %.2f-%.2f); nominal p = %.3f",
      count_hr, count_low, count_high, count_p
    ),
    x = NULL, y = "Recurrence-free survival", colour = NULL
  ) +
  theme_paper +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7.0),
    legend.key.width = grid::unit(10, "pt"),
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  )

p_c_risk <- ggplot(
  count_risk,
  aes(time, programme_burden, label = n_risk, colour = programme_burden)
) +
  geom_text(size = 2.35, family = "Arial", show.legend = FALSE) +
  scale_colour_manual(values = count_colours, drop = FALSE) +
  scale_x_continuous(breaks = risk_times, expand = expansion(mult = c(0.04, 0.01))) +
  coord_cartesian(xlim = c(0, 8)) +
  labs(x = "Years", y = "Number at risk") +
  theme_classic(base_size = 7.6, base_family = "Arial") +
  theme(
    axis.text.y = element_text(colour = col_text, size = 6.7),
    axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    plot.margin = margin(0, 5.5, 5.5, 5.5)
  )
p_c <- p_c_curve / p_c_risk + plot_layout(heights = c(3.4, 1.1))

# Tables --------------------------------------------------------------------
fwrite(
  score_long,
  file.path(table_dir, "Figure_2_patient_programme_scores_by_MEvotype.csv")
)
fwrite(
  programme_counts,
  file.path(table_dir, "Figure_2_programme_count_composition_by_MEvotype.csv")
)
fwrite(
  comparison_labels[, .(
    programme, W, p_value,
    method = "Two-sided Wilcoxon rank-sum test; MEvotype1 versus MEvotype2"
  )],
  file.path(table_dir, "Figure_2_MEvotype_score_tests.csv")
)
count_factor_sm <- summary(fit_count_factor)
count_factor_terms <- rownames(count_factor_sm$coefficients)
count_survival_results <- rbindlist(list(
  data.table(
    comparison = "Per one-programme increase", n = nrow(dat),
    events = sum(dat$rfs_event), hazard_ratio = count_hr,
    conf_low = count_low, conf_high = count_high, p_value = count_p
  ),
  rbindlist(lapply(seq_along(count_factor_terms), function(i) {
    data.table(
      comparison = c("1 versus 0 programmes", "2 versus 0 programmes")[i],
      n = nrow(dat), events = sum(dat$rfs_event),
      hazard_ratio = unname(count_factor_sm$conf.int[i, "exp(coef)"]),
      conf_low = unname(count_factor_sm$conf.int[i, "lower .95"]),
      conf_high = unname(count_factor_sm$conf.int[i, "upper .95"]),
      p_value = unname(count_factor_sm$coefficients[i, "Pr(>|z|)"])
    )
  }))
))
fwrite(
  count_survival_results,
  file.path(table_dir, "Figure_2_programme_count_survival.csv")
)
fwrite(
  count_summary[, .(programme_burden, n, events)],
  file.path(table_dir, "Figure_2_programme_count_survival_groups.csv")
)
fwrite(
  data.table(
    programme = c("Neuronal", "RhoGTPase"),
    absolute_cutoff = c(NEURONAL_CUTOFF, RHOGTPASE_CUTOFF),
    outcome_used_to_select_cutoff = FALSE
  ),
  file.path(table_dir, "fixed_cutoff_definitions.csv")
)

figure_2 <- p_a / (p_b | p_c) + plot_layout(heights = c(1.02, 1.12))
ggsave(
  file.path(figure_dir, "Figure_2_programme_cutoffs_discrimination.pdf"),
  figure_2, width = 13.5, height = 9.2, device = cairo_pdf,
  family = "Arial", bg = "white"
)
ggsave(
  file.path(figure_dir, "Figure_2_programme_cutoffs_discrimination.png"),
  figure_2, width = 13.5, height = 9.2, dpi = 400, bg = "white"
)
message("Generated Figure 2 (MEvotype1 versus MEvotype2)")
