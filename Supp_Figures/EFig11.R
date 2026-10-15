## Relevant Extended Data Fig 11 panels


library(data.table)
library(ggpubr)
library(tidyverse)

metastatic <- fread("data/canary_epidish_w_buffycoat_TMC_deconvolve.tsv") %>%
  left_join(fread("data/twist_canary_stage_tumfrac.csv"), by = c("sample" = "CMDL_ID")) %>%
  filter(Stage == "met")

fraction_vs_tumour <- function(cell_type, label) {
  ggplot(metastatic, aes(x = .data[[cell_type]], y = ichorCNA_TF)) +
    geom_point(size = 3) +
    geom_smooth(method = lm, formula = y ~ x) +
    stat_cor(size = 8) +
    labs(title = paste("Tumour Fraction vs", label, "Fraction"),
         x = paste(label, "Fraction"), y = "Tumour Fraction") +
    theme_minimal() +
    theme(text = element_text(size = 18))
}

panels <- c(a = "luminal", b = "basal", c = "endo", d = "fibroblast", e = "immune", f = "buffycoat")
labels <- c(luminal = "Luminal", basal = "Basal", endo = "Endo", fibroblast = "Fibro",
            immune = "Immune", buffycoat = "Buffycoat")

for (p in names(panels)) {
  ggsave(paste0("EFig11", p, "_", panels[[p]], "_vs_tumour_fraction.pdf"),
         fraction_vs_tumour(panels[[p]], labels[[panels[[p]]]]), width = 6, height = 6)
}
