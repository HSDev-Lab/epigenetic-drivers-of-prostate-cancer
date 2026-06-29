library(Matrix)
library(stringr)
library(data.table)
library(tidyverse)
library(ggpubr)
mat <- fread("/Users/tobymilne-clark/Documents/PhD_Thesis/Data/plasma.samples.beta.mtx")

colnames(mat) <- sub("([A-Za-z0-9]+)([0-9]{2})$", "\\1_\\2", colnames(mat))

plasma_samples <- mat

annotation_col = data.frame(
  group = str_extract(colnames(mat)[-1], "[^_]+")
) 
rownames(annotation_col) = colnames(mat)[-1]


Heatmap(as.matrix(mat[,-1]),
        cluster_columns = TRUE,
        col = viridis::plasma(100),
        raster_quality = 10,
        top_annotation = HeatmapAnnotation(df = annotation_col, 
                                           col = list(group = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B"))))



mat <- as.data.frame(mat[,-1])
labels = labels <- sub("_.*", "", colnames(mat))
labels_bin <- ifelse(labels == "ben", "Ben", "Tumour")
labels_bin  <- factor(labels_bin )
table(labels)
score <- colMeans(mat)

kruskal.test(score ~ labels)
pairwise.wilcox.test(score, labels)

install.packages("clinfun")
library(clinfun)

ord <- as.numeric(factor(labels,
                         levels = c("ben","cpg1","cpg5","met")))

jonckheere.test(score, ord, alternative = "increasing")


cor.test(score, ord, method = "spearman")

score_df <- data.frame(
  sample = colnames(mat),
  group  = labels,
  score = score
)

DMR_icgc$region <- paste0(DMR_icgc$seqnames.x,":", DMR_icgc$start.x,"-", DMR_icgc$end.x)

ggplot(score_df, aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(x = NULL, y = "Sum(DMR)", title = "Probability gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B"))



ggplot(score_df, aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(x = NULL, y = "P(tumour | Ben-vs-Met model)", title = "Probability gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B"))



tumour_idx <- labels %in% c("cpg1","cpg5","met")
ben_idx    <- labels == "ben"

mean_tumour <- rowMeans(mat[, tumour_idx])
mean_ben    <- rowMeans(mat[, ben_idx])

direction <- sign(mean_tumour - mean_ben)

score_dir <- colMeans(mat * direction)

pairwise.wilcox.test(score_dir, labels)


jonckheere.test(score_dir, ord, alternative = "increasing")

cor.test(score_dir, ord, method = "spearman")


score_df <- data.frame(
  sample = colnames(mat),
  group  = labels,
  score = score_dir
) %>% mutate(group_2 = case_when(
  group == "ben" ~ "ben",
  group == "cpg1" ~ "localised",
  group == "cpg5" ~ "localised",
  group == "met" ~ "met"
))

ggplot(score_df, aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  ylim(-10,40) +
  labs(x = NULL, y = "Methyaltion score", title = "Sum of methyaltion gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B")) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ben", label.y.npc = 0.95, 
                     label = "p.format")

ggplot(score_df, aes(group_2, score, fill = group_2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  ylim(-10,40) +
  labs(x = NULL, y = "Methyaltion score", title = "Sum of methyaltion gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#EAEAEA", "localised" = "#B9B9BA" ,"met" = "#717272")) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ben", label.y.npc = 0.95, 
                     label = "p.format") + 
  theme(legend.position = "none")

score_no_met <- score_df %>% filter(group != "met")
ggplot(score_no_met , aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(x = NULL, y = "Methyaltion score", title = "Sum of methyaltion gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B")) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ben", label.y.npc = 0.9, 
                     label = "p.format")

# Microenviromental score

mat_micro <- fread("/Users/tobymilne-clark/Documents/PhD_Thesis/Data/plasma.samples.beta.microenv.mtx")

colnames(mat_micro) <- sub("([A-Za-z0-9]+)([0-9]{2})$", "\\1_\\2", colnames(mat_micro))

mat_micro <- as.data.frame(mat_micro[,-1])

micro_mean_tumour <- rowMeans(mat_micro[, tumour_idx])
micro_mean_ben    <- rowMeans(mat_micro[, ben_idx])

micro_direction <- sign(micro_mean_tumour - micro_mean_ben)

micro_score_dir <- colMeans(mat_micro * micro_direction)

pairwise.wilcox.test(micro_score_dir, labels)


jonckheere.test(micro_score_dir, ord, alternative = "increasing")

cor.test(micro_score_dir, ord, method = "spearman")

micro_score_df <- data.frame(
  sample = colnames(mat_micro),
  group  = labels,
  score = micro_score_dir
) %>% mutate(group_2 = case_when(
  group == "ben" ~ "ben",
  group == "cpg1" ~ "localised",
  group == "cpg5" ~ "localised",
  group == "met" ~ "met"
))

ggplot(micro_score_df , aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(x = NULL, y = "Methyaltion score", title = "Sum of methyaltion gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B")) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ben", label.y.npc = 0.9, 
                     label = "p.format")

ggplot(micro_score_df, aes(group_2, score, fill = group_2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(x = NULL, y = "Methyaltion score", title = "Sum of methyaltion gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#EAEAEA", "localised" = "#B9B9BA" ,"met" = "#717272")) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ben", label.y.npc = 0.95, 
                     label = "p.format") + 
  theme(legend.position = "none")

# Luminal score

mat_lum <- fread("/Users/tobymilne-clark/Documents/PhD_Thesis/Data/plasma.samples.beta.luminal.mtx")

colnames(mat_lum) <- sub("([A-Za-z0-9]+)([0-9]{2})$", "\\1_\\2", colnames(mat_lum))

mat_lum<- as.data.frame(mat_lum[,-1])

lum_mean_tumour <- rowMeans(mat_lum[, tumour_idx])
lum_mean_ben    <- rowMeans(mat_lum[, ben_idx])

lum_direction <- sign(lum_mean_tumour - lum_mean_ben)

lum_score_dir <- colMeans(mat_lum * lum_direction)

pairwise.wilcox.test(lum_score_dir, labels)


jonckheere.test(lum_score_dir, ord, alternative = "increasing")

cor.test(lum_score_dir, ord, method = "spearman")

lum_score_df <- data.frame(
  sample = colnames(mat_lum),
  group  = labels,
  score = lum_score_dir
) %>% mutate(group_2 = case_when(
  group == "ben" ~ "ben",
  group == "cpg1" ~ "localised",
  group == "cpg5" ~ "localised",
  group == "met" ~ "met"
))

ggplot(lum_score_df , aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(x = NULL, y = "Methyaltion score", title = "Sum of methyaltion gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B")) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ben", label.y.npc = 0.9, 
                     label = "p.format")


ggplot(lum_score_df, aes(group_2, score, fill = group_2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  ylim(-10,40) +
  labs(x = NULL, y = "Methyaltion score", title = "Sum of methyaltion gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#EAEAEA", "localised" = "#B9B9BA" ,"met" = "#717272")) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ben", label.y.npc = 0.95, 
                     label = "p.format") + 
  theme(legend.position = "none")

Rho_GTPase_genes <- c(  "ARHGAP25", "ARHGAP26", "ARHGEF7", "BAIAP2", "BCR", "EPHA2",   "FAM13A", "FARP1", "ITGB1", "MCF2", "PAK6", "PLEKHG3",   "SRGAP3", "TAOK3", "VAV2", "VAV3")
neuronal_genes <- c("NRXN1", "NRXN3", "NTRK3", "SLITRK3", "SLITRK4", "GRIN2D", "GRIA3", "SHANK2", "CAMK2A", "CAMK2B", "CAMKK2", "GABRA5", "GABRB1", "GABRG3", "SLC32A1", "SYT7", "ADCY1", "RASGRF1", "CHAT", "KIF17", "PPFIBP2")



plasma_samples

Rho_DMR <- DMR_icgc %>%filter(gene_name %in% Rho_GTPase_genes) %>% distinct(gene_name, .keep_all = TRUE) %>% pull(region)



rho_mat <- plasma_samples %>% filter(V1 %in% Rho_DMR) %>% dplyr::select(-V1) %>% as.data.frame()


rho_mean_tumour <- rowMeans(rho_mat[, tumour_idx])
rho_mean_ben    <- rowMeans(rho_mat[, ben_idx])

rho_direction <- sign(rho_mean_tumour - rho_mean_ben)

rho_score_dir <- colMeans(rho_mat* rho_direction)

pairwise.wilcox.test(rho_score_dir, labels)


jonckheere.test(rho_score_dir, ord, alternative = "increasing")

cor.test(rho_score_dir, ord, method = "spearman")

rho_score_df <- data.frame(
  sample = colnames(rho_mat ),
  group  = labels,
  score = rho_score_dir
)

ggplot(rho_score_df , aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(x = NULL, y = "Methyaltion score", title = "Sum of methyaltion gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B")) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ben", label.y.npc = 0.9, 
                     label = "p.format")


neur_DMR <- DMR_icgc %>%filter(gene_name %in% neuronal_genes) %>% distinct(gene_name, .keep_all = TRUE) %>% pull(region)



neur_mat <- plasma_samples %>% filter(V1 %in% neur_DMR) %>% dplyr::select(-V1) %>% as.data.frame()


neur_mean_tumour <- rowMeans(neur_mat[, tumour_idx])
neur_mean_ben    <- rowMeans(neur_mat[, ben_idx])

neur_direction <- sign(neur_mean_tumour - neur_mean_ben)

neur_score_dir <- colMeans(neur_mat* neur_direction)

pairwise.wilcox.test(neur_score_dir, labels)


jonckheere.test(neur_score_dir, ord, alternative = "increasing")

cor.test(neur_score_dir, ord, method = "spearman")

neur_score_df <- data.frame(
  sample = colnames(neur_mat),
  group  = labels,
  score = neur_score_dir
)

ggplot(neur_score_df , aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(x = NULL, y = "Methyaltion score", title = "Sum of methyaltion gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B")) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ben", label.y.npc = 0.9, 
                     label = "p.format")






deconv <- read_tsv("/Users/tobymilne-clark/Documents/PhD_Thesis/Data/Epidish/canary_epidish_w_buffycoat_TMC_deconvolve.tsv")
tum_frac <- read_csv("/Users/tobymilne-clark/Documents/PhD_Thesis/Data/Epidish/twist_canary_stage_tumfrac.csv")
CMDL <- read_csv("/Users/tobymilne-clark/Documents/PhD_Thesis/Data/Epidish/CMDL_filenames.txt") 
CMDL$group_CMDL <- c(rep("ben", 20), rep("cpg1", 16), rep("cpg5", 16), rep("met", 12),"tmp")
CMDL$CMDL_ID <- sub("_.*", "", CMDL$x)


cbind(CMDL[-65,],micro_score_df) %>% left_join(deconv, by = c("CMDL_ID"="sample" )) %>% 
  left_join(tum_frac, by = c("CMDL_ID" )) %>% mutate(microenviroment = immune  + endo + fibroblast) %>%
  ggplot(aes(micro_score_dir, microenviroment, color = group_CMDL)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = "Microenviromental score", y = "B cell proportion", title = "Correlation between microenviromental score and B cell proportion") +
  scale_color_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B")) +
  stat_cor(method = "spearman", label.x.npc = 0.1, label.y.npc = 0.9)

   
cbind(CMDL[-65,],lum_score_df) %>% left_join(deconv, by = c("CMDL_ID"="sample" )) %>% 
  left_join(tum_frac, by = c("CMDL_ID" )) %>%
  ggplot(aes(lum_score_dir, luminal, color = group_CMDL)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = "Luminal Methylation score", y = "Methylation proportion", title = "Correlation between Luminal score and Luminal cell proportion") +
  scale_color_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B")) +
  stat_cor(method = "spearman", label.x.npc = 0.1, label.y.npc = 0.9)

cbind(CMDL[-65,],lum_score_df) %>% left_join(deconv, by = c("CMDL_ID"="sample" )) %>% 
  left_join(tum_frac, by = c("CMDL_ID" )) %>%
  ggplot(aes(lum_score_dir, ichorCNA_TF, color = group_CMDL)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = "Luminal Methylation score", y = "Methylation proportion", title = "Correlation between Luminal score and tumour fraction") +
  scale_color_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B")) +
  stat_cor(method = "spearman", label.x.npc = 0.1, label.y.npc = 0.9)




library(tidyverse)
CANCAP03 <- read_tsv("/Users/tobymilne-clark/Documents/PhD_Thesis/Data/Epidish/CANCAP03plasma_epidish_w_buffycoat_TMC_deconvolve.tsv")


library(dplyr)
library(stringr)

CANCAP03 <- CANCAP03 %>%
  mutate(
    stage = case_when(
      str_ends(sample, "Surgery")  ~ "Surgery",
      str_ends(sample, "FollowUp") ~ "FollowUp",
      TRUE                             ~ "Screen"
    )
  )  %>%
  mutate(
    stage = factor(stage,
                   levels = c("Screen", "Surgery", "FollowUp"),
                   ordered = TRUE)
  ) %>% 
  mutate(sampleID = str_extract(sample, "(?<=_).{4}"))

CANCAP03 %>% 
  pivot_longer(-c(sample,stage,sampleID), names_to = "cell", values_to = "proportion") %>%
  ggplot(aes(x=reorder(sample, proportion), y=proportion, fill=cell)) + geom_bar(stat="identity") +
  scale_colour_manual(
    values = c(
      "basal" = "#292BC1",
      "luminal" = "#9DCEE2",
      "fibroblast" = "#5CAE63",
      "endo" = "#986EAC",
      "immune" = "#CB1B16",
      "smoothmuscle" = "#F8A627",
      "buffycoat" = "#FEDFD4"
    )
  )+
  scale_fill_manual(
    values = c(
      "basal" = "#292BC1",
      "luminal" = "#9DCEE2",
      "fibroblast" = "#5CAE63",
      "endo" = "#986EAC",
      "immune" = "#CB1B16",
      "smoothmuscle" = "#F8A627",
      "buffycoat" = "#FEDFD4"
    )
  )+
  theme_bw() +  
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"), 
        text = element_text(size = 12),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size =6),
        axis.text.y = element_text()) +  
  facet_grid(.~stage, scales="free", space="free") +
  xlab("Samples") + ylab("Cell Proportion") +  
  labs(fill = "Cell Type")  



CANCAP03 %>% 
  pivot_longer(-c(sample,stage,sampleID), names_to = "cell", values_to = "proportion") %>%
  ggplot(aes(x=stage, y=proportion, fill=cell)) + geom_boxplot() +
  scale_colour_manual(
    values = c(
      "basal" = "#292BC1",
      "luminal" = "#9DCEE2",
      "fibroblast" = "#5CAE63",
      "endo" = "#986EAC",
      "immune" = "#CB1B16",
      "smoothmuscle" = "#F8A627",
      "buffycoat" = "#FEDFD4"
    )
  )+
  scale_fill_manual(
    values = c(
      "basal" = "#292BC1",
      "luminal" = "#9DCEE2",
      "fibroblast" = "#5CAE63",
      "endo" = "#986EAC",
      "immune" = "#CB1B16",
      "smoothmuscle" = "#F8A627",
      "buffycoat" = "#FEDFD4"
    )
  )+
  theme_bw() +  
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"), 
        text = element_text(size = 12),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size =6),
        axis.text.y = element_text()) +  
  facet_grid(.~stage, scales="free", space="free") +
  xlab("Samples") + ylab("Cell Proportion") +  
  labs(fill = "Cell Type")  




CANCAP03 %>% 
  pivot_longer(-c(sample,stage), names_to = "cell", values_to = "proportion") %>%
  filter(cell == "luminal") %>%
  ggplot(aes(x=stage, y=proportion, fill=cell)) + geom_boxplot() +
  scale_colour_manual(
    values = c(
      "basal" = "#292BC1",
      "luminal" = "#9DCEE2",
      "fibroblast" = "#5CAE63",
      "endo" = "#986EAC",
      "immune" = "#CB1B16",
      "smoothmuscle" = "#F8A627",
      "buffycoat" = "#FEDFD4"
    )
  )+
  scale_fill_manual(
    values = c(
      "basal" = "#292BC1",
      "luminal" = "#9DCEE2",
      "fibroblast" = "#5CAE63",
      "endo" = "#986EAC",
      "immune" = "#CB1B16",
      "smoothmuscle" = "#F8A627",
      "buffycoat" = "#FEDFD4"
    )
  )+
  theme_bw() +  
  theme(panel.background=element_blank(), 
        axis.text=element_text(color="black"), 
        text = element_text(size = 12),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size =6),
        axis.text.y = element_text()) +  
  xlab("Samples") + ylab("Cell Proportion") +  
  labs(fill = "Cell Type")  +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "Surgery", label.y.npc = 0.9, 
                     label = "p.format")


ggplot(CANCAP03,
       aes(stage, luminal, group = sampleID)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  stat_summary(aes(group = 1),
               fun = mean,
               geom = "line",
               linewidth = 1.5) +
  theme_classic() +
  labs(y = "Luminal proportion")

ggplot(CANCAP03,
       aes(stage, immune, group = sampleID)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  stat_summary(aes(group = 1),
               fun = mean,
               geom = "line",
               linewidth = 1.5) +
  theme_classic() +
  labs(y = "Immune proportion")


ggplot(CANCAP03,
       aes(stage, basal, group = sampleID)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  stat_summary(aes(group = 1),
               fun = mean,
               geom = "line",
               linewidth = 1.5) +
  theme_classic() +
  labs(y = "basal proportion")






ggplot(CANCAP03,
       aes(stage, luminal)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = sampleID), alpha = 0.2) +
  geom_point()

delta_df <- CANCAP03 %>%
  pivot_wider(names_from = stage,
              values_from = luminal) %>%
  mutate(
    delta_surg = Surgery - Screen,
    delta_follow = FollowUp - Surgery
  )

ggplot(delta_df, aes("", delta_surg)) +
  geom_boxplot() +
  geom_jitter(width = 0.1)


friedman.test(luminal ~ stage | sampleID,
              data = CANCAP03)

pairwise.wilcox.test(
  CANCAP03$luminal,
  CANCAP03$stage,
  paired = TRUE,
  p.adjust.method = "holm"
)

library(lme4)

model <- lmer(luminal ~ stage + (1 | sampleID),
              data = CANCAP03)

anova(model)


CANCAP_meth <- read_tsv("/Users/tobymilne-clark/Documents/PhD_Thesis/Data/CANCAP_plasma_ave_beta.tsv")


meta_cancap <- CANCAP_meth %>% mutate(region = paste0("chr",chr,":",start,"-",end)) %>% dplyr::select(-c(chr,start,end,difference,gene)) %>%
  pivot_longer(-region, names_to = "sample", values_to = "beta")  %>%
  mutate(
    stage = case_when(
      str_ends(sample, "Surgery_beta")  ~ "Surgery",
      str_ends(sample, "FollowUp_beta") ~ "FollowUp",
      TRUE                             ~ "Screen"
    )
  )  %>%
  mutate(
    stage = factor(stage,
                   levels = c("Screen", "Surgery", "FollowUp"),
                   ordered = TRUE)
  ) %>% 
  mutate(sampleID = str_extract(sample, "(?<=_).{4}")) %>%
  
  ggplot(aes(stage, beta, group = sampleID)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  stat_summary(aes(group = 1),
               fun = mean,
               geom = "line",
               linewidth = 1.5) 

cancap_meta <- CANCAP_meth %>% mutate(region = paste0("chr",chr,":",start,"-",end)) %>% dplyr::select(-c(chr,start,end,difference,gene,region)) %>% 
  t() %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(
    stage = case_when(
      str_ends(sample, "Surgery_beta")  ~ "Surgery",
      str_ends(sample, "FollowUp_beta") ~ "FollowUp",
      TRUE                             ~ "Screen"
    )
  )  %>%
  mutate(
    stage = factor(stage,
                   levels = c("Screen", "Surgery", "FollowUp"),
                   ordered = TRUE)
  ) %>% 
  mutate(sampleID = str_extract(sample, "(?<=_).{4}")) %>% dplyr::select(sample, stage, sampleID) %>%
  mutate(relapse = ifelse(sampleID %in% c(1009,1012,1015),"relpase","no relapse"))

cancap_region <- CANCAP_meth %>% mutate(region = paste0("chr",chr,":",start,"-",end)) %>% dplyr::select(-c(chr,start,end,difference,gene)) 

cancap_neur_mat <- cancap_region %>% filter(region %in% neur_DMR) %>% dplyr::select(-region) %>% as.data.frame()


cancap_neur_score_dir <- colMeans(cancap_neur_mat * neur_direction)

pairwise.wilcox.test(neur_score_dir, cancap_meta$stage)

labels_cancap <- sub("_.*", "", colnames(cancap_neur_mat))

jonckheere.test(neur_score_dir, ord, alternative = "increasing")

cor.test(neur_score_dir, ord, method = "spearman")

cancap_neur_score_df <- data.frame(
  sample = colnames(cancap_neur_mat),
  group  = cancap_meta$stage,
  sampleID = cancap_meta$sampleID,
  score = cancap_neur_score_dir,
  relapse = cancap_meta$relapse
)

ggplot(cancap_neur_score_df , aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(x = NULL, y = "Methyaltion score", title = "Sum of methyaltion gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B")) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ben", label.y.npc = 0.9, 
                     label = "p.format")


cancap_rho_mat <- cancap_region %>% filter(region %in% Rho_DMR) %>% dplyr::select(-region) %>% as.data.frame()


cancap_rho_score_dir <- colMeans(cancap_rho_mat* rho_direction)

cancap_rho_score_df <- data.frame(
  sample = colnames(cancap_rho_mat ),
  group  = cancap_meta$stage,
  sampleID = cancap_meta$sampleID,
  score = cancap_rho_score_dir,
  relapse = cancap_meta$relapse
  
)

ggplot(cancap_rho_score_df , aes(group, score, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(x = NULL, y = "Methyaltion score", title = "Sum of methyaltion gradient across disease states") +
  scale_fill_manual(values = c("ben" = "#27A9E1", "cpg1" = "#FBB041", "cpg5" = "#F8A592" ,"met" = "#F15A2B")) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "ben", label.y.npc = 0.9, 
                     label = "p.format")


ggplot(cancap_neur_score_df,
       aes(group, score, group = sampleID,color = relapse)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  stat_summary(aes(group = 1),
               fun = mean,
               geom = "line",
               linewidth = 1.5) +
  theme_classic() +
  labs(y = "Neuronal proportion")


ggplot(cancap_rho_score_df,
       aes(group, score, group = sampleID,colour = relapse)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 2) +
  stat_summary(aes(group = 1),
               fun = mean,
               geom = "line",
               linewidth = 1.5) +
  theme_classic() +
  labs(y = "RhoGTPase score")





deseq_genes_OVsOD_3 
library(tidyverse)
deseq_genes_PreVsOD_3 %>% as.data.frame() %>% rownames_to_column("ENSEMBL") %>% left_join(Table_Genes, by = "ENSEMBL") %>% write_csv("/Users/tobymilne-clark/Documents/Data/Atlas/PPI_Networks/CANCAP03/CANCAP03_DEG_Pre_vs_OD_buld_RNA.csv")



















