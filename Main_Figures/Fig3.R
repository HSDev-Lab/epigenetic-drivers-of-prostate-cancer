## Relevant Fig 3 panels

library(data.table)
library(annotatr)
library(ReactomePA)
library(org.Hs.eg.db)
library(tidyverse)

hyper <- fread("data/all_patients_ave_beta_hiQC.summary.hyper.celltype.tsv")
hypo  <- fread("data/all_patients_ave_beta_hiQC.summary.hypo.celltype.tsv")
all_dmrs <- rbind(hyper, hypo)


## ---- e: cell type assignment of the DMRs ------------------------------------

# a DMR counts towards every cell type it is assigned to; DMRs not assigned to any one are dropped
unassigned <- "Basal, Endo, Fibroblast, Immune, Luminal, SmoothMuscle"
cell_types <- c("Luminal", "Basal", "Immune", "Endo", "Fibroblast", "SmoothMuscle")

count_assignments <- function(dmrs) {
  assigned <- table(dmrs$cell_assigned) %>%
    as.data.frame() %>%
    dplyr::filter(Var1 != unassigned)
  sapply(cell_types, function(ct) sum(assigned$Freq[str_detect(assigned$Var1, ct)]))
}

assignment_counts <- data.frame(Celltype = cell_types,
                                Hyper = count_assignments(hyper),
                                Hypo = count_assignments(hypo)) %>%
  pivot_longer(c(Hyper, Hypo), names_to = "Type", values_to = "Count")

fig3e <- ggplot(assignment_counts, aes(x = Celltype, y = Count, fill = Type)) +
  geom_col() +
  scale_fill_manual(values = c("Hyper" = "#F8A627", "Hypo" = "#803E97")) +
  labs(x = "Cell type", y = "Count", fill = "Type") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Fig3e_celltype_assignments.pdf", fig3e, width = 6, height = 5)


## ---- f, g: Reactome pathways of the immune- and luminal-assigned genes ------

# the EPIC array genes are the background
annotations <- build_annotations(genome = "hg38", annotations = c("hg38_basicgenes", "hg38_cpgs"))

epic_panel <- fread("data/EPIC_panel_hg38.bed")
colnames(epic_panel) <- c("chr", "start", "end")
epic_annot <- annotate_regions(makeGRangesFromDataFrame(epic_panel), annotations = annotations,
                               ignore.strand = TRUE, quiet = TRUE) %>% as.data.frame()
universe <- clusterProfiler::bitr(unique(na.omit(epic_annot$annot.symbol)), fromType = "SYMBOL",
                                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

reactome_barplot <- function(cell_type, n_pathways, title) {
  genes <- all_dmrs %>% dplyr::filter(top1 == cell_type) %>% dplyr::pull(gene) %>% unique()
  entrez <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db)$ENTREZID
  enrichPathway(gene = entrez, organism = "human", pvalueCutoff = 0.05,
                readable = TRUE, universe = universe) %>%
    barplot(showCategory = n_pathways) +
    ggtitle(title) +
    theme(axis.text.y = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, size = 22))
}

ggsave("Fig3f_immune_Reactome.pdf", reactome_barplot("Immune", 5, "Immune assigned genes"), width = 7, height = 8)
ggsave("Fig3g_luminal_Reactome.pdf", reactome_barplot("Luminal", 6, "Luminal assigned genes"), width = 7, height = 8)
