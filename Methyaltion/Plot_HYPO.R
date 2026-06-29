#############################
### RNA-seq prep 
#############################
library(preprocessCore)

load('/Users/tobymilne-clark/Library/Mobile Documents/com~apple~CloudDocs/PhD/Data/ICGC/RNA/salmon_gene_dataset.RData')
df <- as.data.frame(se$counts)
colnames(df) <- PrelimSampleInfo$File

matched_samples <- readxl::read_xlsx("/Users/tobymilne-clark/Documents/Documents/ICGC/scRNA/Samples_RNA_Meth_matched.xlsx")
keep <- matched_samples$Exclude == "Keep"
cts_matched = df[,matched_samples$File.y[keep]] 
cts_matched = cts_matched[rowSums(cts_matched) > 0,] 

expression.threshold <- 5
cts.filtered = cts_matched
cts.filtered[cts.filtered<expression.threshold]<-expression.threshold
cts.filtered = cts.filtered[rowSums(cts.filtered)>expression.threshold*ncol(cts.filtered),]

cts.qnorm.filtered=data.frame(normalize.quantiles(as.matrix(cts.filtered)),row.names=rownames(cts.filtered))
colnames(cts.qnorm.filtered)=colnames(cts.filtered)


suppressPackageStartupMessages(require(biomaRt))
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

Table_Genes = getBM(attributes=c("ensembl_gene_id",'entrezgene_id',"external_gene_name","description"), 
                    values=rownames(cts_2), 
                    mart=ensembl) # Creating gene reference table

colnames(Table_Genes) <- c("ENSEMBL","ENTREZ","SYMBOL","description")




neuronal_genes<- c(  "KCNAB2", "PTPRF", "SLC1A7", "CHAT", "KCNMA1", "SLC18A2", "PPFIBP2", "KCNA4",
                     "SYT7", "KCNK7", "SHANK2", "GRIK4", "ABCC9", "ALDH2", "NBEA", "ADCY4", "NRXN3",
                     "GABRA5", "RASGRF1", "NTRK3", "ADCY9", "KCNJ12", "CACNG4", "GNAL", "TUBB6",
                     "HCN2", "GNG7", "PTPRS", "CACNA1A", "LRFN1", "GRIN2D", "PPFIA3", "KCNS3",
                     "NRXN1", "EPB41L1", "SLC32A1", "KCNS1", "KCNB1", "KCNQ2", "KCNJ15", "SLC6A11",
                     "GRIP2", "CACNA2D2", "CACNA2D3", "ADCY5", "SLITRK3", "GABRB1", "ADCY2",
                     "KCNN2", "CAMK2A", "KCNK17", "KCNQ5", "RPS6KA2", "PRKAR1B", "CAMK2B", "ADCY1",
                     "KCNH2", "PRKAG2", "KCNV1", "GRIA3", "SLITRK4", "KIF17", "KCNK2", "KCNA1",
                     "PRKAG1", "CAMKK2", "GABRG3", "CACNG3", "KCNG2", "NRG1")

DMR_icgc <- read_csv("/Users/tobymilne-clark/Documents/Data/ICGC/DMRs_annotated.csv")


dmrs_withoutsamplesep <- DMR_icgc  %>% dplyr::rename(seqnames = seqnames.x, start = start.x, end = end.x)
dmrs_withoutsamplesep.dt <- dmrs_withoutsamplesep %>% data.table::as.data.table()
data.table::setkey(dmrs_withoutsamplesep.dt, seqnames, start, end)  

consensus_dt <- data.table::as.data.table(consensus_df)  
data.table::setkey(consensus_dt, seqnames, start, end)  

DMRs_CHIP_peaks_overlap  <- data.table::foverlaps(consensus_dt,dmrs_withoutsamplesep.dt, type="any",nomatch=0L)



chr11
67583459
67583520



chr19
48414873
48414948

chr5:112737838-112737916
icgc <- read_tsv("/Users/tobymilne-clark/Documents/Data/DMR_calling/ICGC/methylation_matrix_10x_90perc.tsv")
icgc_region <- icgc %>%
  dplyr::filter(chr == "chr19" & start >= 112735838 & start <= 112739916) 

region <- as.matrix(icgc_region %>%
                      dplyr::select(-chr, -start, -strand))

labels <- icgc_region$start
cpg_pos <- icgc_region$start
region_chr <- "chr19"

# -- 1) subset & covariance
S      <- cov(t(region))
n      <- nrow(S)

# -- 2) mask *upper* triangle + diagonal (keep strict lower),
#        then flip vertically just like np.flipud
S_masked    <- S
S_masked[upper.tri(S_masked, diag = TRUE)] <- NA
S_flip      <- S_masked[n:1, , drop = FALSE]

# -- 3) rebuild the Python mesh exactly:
rows   <- n:0
cols   <- 0:n
# Python: for r in rows, for c in cols: coords.append((c, r))
pts <- do.call(rbind, lapply(rows, function(r) {
  data.frame(col = cols, row = r)
}))
t_mat <- matrix(c(
  1,   0.5,   # [ [1,  0.5],
  -1,   0.5    #   [-1, 0.5] ]
), nrow = 2, byrow = TRUE)

xy <- as.matrix(pts) %*% t_mat
# According to Python pcolormesh:
#   X_mesh = coordinate_matrix[:,1]  = xy[,2]
#   Y_mesh = coordinate_matrix[:,0]  = xy[,1]
X_raw <- xy[, 2]
Y_raw <- xy[, 1]

# reshape into (n+1)x(n+1) so [i, j] matches Python’s grid
Xgrid <- matrix(X_raw, nrow = n+1, byrow = TRUE)
Ygrid <- matrix(Y_raw, nrow = n+1, byrow = TRUE)

# -- 4) rescale Xgrid from [0, n] → [1, n] so it lines up with CpG ticks
Xgrid <- (Xgrid - min(Xgrid)) / (max(Xgrid) - min(Xgrid)) * (n - 1) + 1

# -- 5) extract each valid cell from the LOWER‐triangle of S_flip
#    (that’s what Python’s flip + mask shows) as a 4-corner diamond:
poly_list <- vector("list", n * (n-1) / 2)
k <- 1
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    val <- S_flip[i, j]
    if (!is.na(val)) {
      xs <- c(
        Xgrid[i,   j],
        Xgrid[i,   j+1],
        Xgrid[i+1, j+1],
        Xgrid[i+1, j]
      )
      ys <- c(
        Ygrid[i,   j],
        Ygrid[i,   j+1],
        Ygrid[i+1, j+1],
        Ygrid[i+1, j]
      )
      poly_list[[k]] <- data.frame(
        x    = xs,
        y    = ys,
        cov  = val,
        cell = k
      )
      k <- k + 1
    }
  }
}
heat_df <- bind_rows(poly_list)

# -- 6) top panel: time-series of each CpG
line_df <- as.data.frame(region) %>%
  mutate(pos = seq_along(labels)) %>%
  pivot_longer(-pos, names_to = "sample", values_to = "value") %>% 
  mutate(group = sub(".*\\_", "", sample))

p1 <- ggplot(line_df, aes(x = pos, y = value, group = sample, colour = group)) +
  geom_line(alpha = 0.3) +
  scale_x_continuous(
    breaks = 1:n,
    labels = labels,
    limits = c(1, n),
    expand = c(0, 0)
  ) +
  geom_smooth(
    aes(group = group),
    #    method = "loess",
    se     = FALSE,
    alpha  = 0.5
  ) +
  ylab("Methylation") +
  theme_minimal(base_size = 9) +
  theme(
    axis.title.x    = element_blank(),
#    axis.text.x     = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.x = element_blank(),
    axis.ticks.x    = element_blank(),
    legend.position = "none"
  ) 


# -- 7) bottom panel: the rotated, tip-down triangle
p2 <- ggplot(heat_df, aes(x = x, y = y, group = cell, fill = cov)) +
  geom_polygon(color = NA) +
  scale_fill_viridis_c(na.value = NA) +
  scale_x_continuous(
    breaks = 1:n,
    labels = labels,
    limits = c(1, n),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = range(heat_df$y),
    expand = c(0, 0)
  ) +
  coord_fixed() +
  theme_void()

# 2) Turn them into gtables:
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)

# 3) Locate the panel columns in each gtable:
panel_col1 <- g1$layout[g1$layout$name == "panel", "l"]
panel_col2 <- g2$layout[g2$layout$name == "panel", "l"]

# 4) Pull out the width spec for p1’s panel, and assign it to p2:
g2$widths[panel_col2] <- g1$widths[panel_col1]


# ─── 0) your CpG info ───────────────────────────────────────────────────────────
#   cpg_pos: numeric vector of genomic positions (bp) for each CpG, length = n
#   labels : character vector of same length for the x-axis tick labels
#   region_chr: the chromosome name (e.g. "chr1")
n   <- length(cpg_pos)

# ─── 1) find all bigwigs & sample names ─────────────────────────────────────────
bw_dir     <- "/Users/tobymilne-clark/Documents/Data/Atlas/Chromatin_and_hic/ENCODE/bigWigs/"
bw_files   <- list.files(bw_dir, pattern = "\\.bigWig$", full.names = TRUE)
sample_ids <- sub("\\_.*", "", tools::file_path_sans_ext(basename(bw_files)))


# ─── 2) define the genomic window covering your CpGs ────────────────────────────
region_gr <- GRanges(
  seqnames = region_chr,
  ranges   = IRanges(start = min(cpg_pos), end = max(cpg_pos))
)

# ─── 3) import each bigWig & extract at CpG positions ───────────────────────────
bw_df_list <- map2(
  bw_files, sample_ids,
  ~ {
    # import only the region as an RleList
    rles <- rtracklayer::import.bw(.x, which = region_gr, as = "RleList")
    sig  <- as.numeric(rles[[region_chr]][cpg_pos])  # extract at each CpG
    
    data.frame(
      pos    = seq_len(n),
      signal = sig,
      sample = .y
    )
  }
) %>% set_names(sample_ids)

# combine into one long data.frame
df_bw <- bind_rows(df_list)

# ─── 4) build the overlaid line plot with end‐labels ────────────────────────────
# compute where to put each sample’s label (at pos = n)
label_df <- df_bw %>%
  group_by(sample) %>%
  filter(pos == max(pos), !is.na(signal))

p_bw <- ggplot(df_bw, aes(x = pos, y = signal, color = sample)) +
  geom_line(size = 0.8) +
  # hang a text label off the right end of each line
  geom_text_repel(
    data            = label_df,
    aes(label = sample),
    nudge_x         = 0.5,
    direction       = "y",
    hjust           = 0,
    segment.size    = 0.2,
    segment.linetype= "dashed",
    size            = 3
  ) +
  scale_x_continuous(
    breaks = seq_len(n),
    labels = labels,
    limits = c(1, n + 1),    # extend a bit for the labels
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(angle = 90, vjust = 0.5, size = 6)
  ) +
  labs(y = "BigWig signal")

common_scale <- scale_x_continuous(
  breaks = seq_len(n),
  labels = labels,
  limits = c(1, n),
  expand = c(0, 0)
)


# ─── build one ggplot per sample ─────────────────────────────────────────────────
p_bw_list <- map(bw_df_list, ~{
  ggplot(.x, aes(x = pos, y = signal)) +
    geom_line() +
    common_scale +
    labs(title = unique(.x$sample), y = "Signal") +
    theme_minimal(base_size = 9) +
    theme(
      axis.title.x    = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      plot.title      = element_text(hjust = 0.5, size = 8),
      plot.margin     = margin(0, 0, 0, 0)
    )
})

peak_dir     <- "/Users/tobymilne-clark/Documents/Data/Atlas/Chromatin_and_hic/ENCODE/peak_plots"
peak_files   <- list.files(peak_dir , pattern = "\\.bed$", full.names = TRUE)
sample_ids <- sub("\\_.*", "", tools::file_path_sans_ext(basename(bw_files)))

load_narrowpeaks <- function(bed_dir,
                             pattern = "\\.(bed|narrowPeak)$") {
  files <- list.files(bed_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    warning("No peak files found in ", bed_dir)
    return(GRangesList())
  }
  
  sample_ids <- tools::file_path_sans_ext(basename(files))
  out_list   <- vector("list", length(files))
  names(out_list) <- sample_ids
  
  for (i in seq_along(files)) {
    dt <- data.table::fread(files[i], header = FALSE, data.table = TRUE)
    nc <- ncol(dt)
    if (nc < 3) stop("File ", files[i], " has fewer than 3 columns")
    
    # ensure at least 6 columns
    if (nc < 4) dt[, V4 := "."]
    if (nc < 5) dt[, V5 := NA_real_]
    if (nc < 6) dt[, V6 := "*"]
    # cap at 10 columns if more
    if (nc > 10) dt <- dt[, .SD, .SDcols = 1:10]
    
    # assign names up to however many columns we have
    cn <- c("chr","start","end","name","score",
            "strand","signalValue","pValue","qValue","peak")
    data.table::setnames(dt, cn[1:ncol(dt)])
    
    # coerce types
    dt[, `:=`(
      start = as.integer(start),
      end   = as.integer(end),
      score = as.numeric(score)
    )]
    # optional extra numeric columns
    if ("signalValue" %in% names(dt)) dt[, signalValue := as.numeric(signalValue)]
    if ("pValue"      %in% names(dt)) dt[, pValue      := as.numeric(pValue)]
    if ("qValue"      %in% names(dt)) dt[, qValue      := as.numeric(qValue)]
    if ("peak"        %in% names(dt)) dt[, peak        := as.integer(peak)]
    
    # force strand to be one of "+", "-", "*"
    dt[, strand := as.character(strand)]
    dt[ ! strand %in% c("+","-","*"), strand := "*" ]
    
    # build GRanges (BED is 0-based, half-open)
    gr <- GRanges(
      seqnames = dt$chr,
      ranges   = IRanges(start = dt$start + 1, end = dt$end),
      strand   = dt$strand
    )
    # attach other columns as metadata
    mcols(gr) <- dt[, setdiff(names(dt), c("chr","start","end","strand")), with = FALSE]
    
    out_list[[i]] <- gr
  }
  
  GRangesList(out_list)
}

peak_list  <- load_narrowpeaks(peak_dir, pattern = ".bed")
names(peak_list) <- sample_ids

peak_idx_list <- lapply(sample_ids, function(sid) {
  gr <- peak_list[[sid]]
  
  # guard: NULL or not GRanges or zero-length → empty tibble
  if (is.null(gr) || !inherits(gr, "GRanges") || length(gr) == 0) {
    return(tibble(start_idx = integer(), end_idx = integer()))
  }
  
  # subset to region_chr and overlapping CpG window
  gr_sub <- gr[ seqnames(gr) == region_chr &
                  end(gr)   >= min(cpg_pos) &
                  start(gr) <= max(cpg_pos) ]
  if (length(gr_sub) == 0) {
    return(tibble(start_idx = integer(), end_idx = integer()))
  }
  
  # map genomic → CpG‐index (1…n)
  tibble(
    start_idx = findInterval(start(gr_sub), cpg_pos),
    end_idx   = findInterval(end(gr_sub),   cpg_pos)
  ) %>%
    filter(start_idx > 0, end_idx > 0)
}) 
names(peak_idx_list) <- sample_ids

p_bw_list_highlight <- map(sample_ids, function(sid) {
  df_bw    <- bw_df_list[[sid]]
  peaks_df <- peak_idx_list[[sid]]
  
  p <- ggplot(df_bw, aes(x = pos, y = signal))
  
  # add grey bands only if there are peaks
  if (nrow(peaks_df) > 0) {
    p <- p +
      geom_rect(data = peaks_df, inherit.aes = FALSE,
                aes(xmin = start_idx - 0.5, xmax = end_idx + 0.5),
                ymin = -Inf, ymax = Inf,
                fill = "grey80", alpha = 0.4)
  }
  
  p + geom_line(color = "steelblue") +
    common_scale +
    labs(title = sid, y = "BigWig") +
    theme_minimal(base_size = 9) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title   = element_text(hjust = 0.5, size = 8),
      plot.margin  = margin(0, 0, 0, 0)
    )
})


all_plots <- c(list(p1), p_bw_list_highlight, list(p2))

# 1) Convert each ggplot to a gtable
gt_list <- lapply(all_plots, ggplotGrob)

# 2) Figure out which column in the gtable layout is the "panel"
#    (they should all be the same, so we just look at the first)
panel_col <- gt_list[[1]]$layout[gt_list[[1]]$layout$name == "panel", "l"]

# 3) Force every other plot to use the same panel width as the first
for (i in seq_along(gt_list)[-1]) {
  gt_list[[i]]$widths[panel_col] <- gt_list[[1]]$widths[panel_col]
}

height_factors <- c(1, rep(0.25, length(p_bw_list_highlight)), 1)

for(i in seq_along(gt_list)) {
  # find which row is the panel
  panel_row <- gt_list[[i]]$layout[
    gt_list[[i]]$layout$name == "panel", "t"
  ]
  # scale its height unit
  gt_list[[i]]$heights[panel_row] <-
    gt_list[[i]]$heights[panel_row] * height_factors[i]
}

# 5) Stack them top‐to‐bottom with rbind (rbind.gtable from gridExtra)
stacked <- gt_list[[1]]
for(i in seq_along(gt_list)[-1]) {
  stacked <- rbind(stacked, gt_list[[i]])
}

# 6) Draw the final aligned stack
grid.newpage()
grid.draw(stacked)

