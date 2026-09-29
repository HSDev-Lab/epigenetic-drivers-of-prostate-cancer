#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(survival)
  library(ggplot2)
  library(patchwork)
  library(GenomicRanges)
  library(IRanges)
  library(SummarizedExperiment)
})

# =============================================================================
# Submission figure rebuild from raw data
# =============================================================================
# This script deliberately does not source prior analysis scripts and does not
# read prior burden-score or survival-result exports. It reconstructs the ICGC
# 142-patient cohort and TCGA-PRAD scores from the primary methylation,
# annotation and clinical inputs listed below.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
output_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
figure_dir <- file.path(output_root, "figures")
panel_dir <- file.path(figure_dir, "panels")
table_dir <- file.path(output_root, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260911)

inputs <- list(
  icgc_all_dmr = paste0(
    "/Users/tobymilne-clark/Documents/Data/DMR_calling/ICGC/",
    "DMRs_all_samples.csv"
  ),
  icgc_curated_dmr = paste0(
    "/Users/tobymilne-clark/Documents/Data/DMR_calling/ICGC/",
    "DMRs_with_coords_350_samples.csv"
  ),
  icgc_analysis_cohort = paste0(
    "/Users/tobymilne-clark/Documents/Codex/2026-08-04/i/work/",
    "mevotype_points_1_2/analysis_cohort.csv"
  ),
  icgc_evotype = paste0(
    "/Users/tobymilne-clark/Documents/Data/ICGC/Mevotypes/",
    "evotype_final_data_v3test.csv"
  ),
  icgc_mevotype = paste0(
    "/Users/tobymilne-clark/Documents/Data/ICGC/Mevotypes/",
    "mevotype_list.csv"
  ),
  icgc_genomic_map = paste0(
    "/Users/tobymilne-clark/Documents/Data/fCpGs/Deconvolution/",
    "ICGC_with_ID.csv"
  ),
  icgc_sample_metadata = paste0(
    "/Users/tobymilne-clark/Documents/Documents/ICGC/",
    "Downsampled_ICGC_EPIC/Metadata_ICGC_EPIC_Downsampling.csv"
  ),
  icgc_clinical = paste0(
    "/Users/tobymilne-clark/Documents/Data/ICGC/Mevotypes/",
    "20230920_ICGC_UK_ClinData.xlsx"
  ),
  neuronal_definition = paste0(
    "/Users/tobymilne-clark/Documents/Paper/FINAL_PAPER_FIGURES/tables/",
    "Neuronal_score_definitions.csv"
  ),
  mevotype_dmr_definition = paste0(
    "/Users/tobymilne-clark/Documents/Paper/FINAL_PAPER_FIGURES/tables/",
    "MEvotype_specific_change_DMRs_717.csv"
  ),
  tcga_se = paste0(
    "/Users/tobymilne-clark/Documents/Codex/2026-08-04/i/outputs/",
    "mevotype_reviewer_points_1_2/TCGA_PRAD_GDC_450K/",
    "TCGA_PRAD_450K_ICGC_DMR_overlap_SE.rds"
  ),
  tcga_survival = paste0(
    "/Users/tobymilne-clark/Documents/Codex/2026-08-04/i/outputs/",
    "mevotype_reviewer_points_1_2/TCGA_PRAD_GDC_450K/",
    "TCGA_PanCancer_CDR_survival.tsv"
  ),
  tcga_clinical = paste0(
    "/Users/tobymilne-clark/Documents/Codex/2026-08-04/i/outputs/",
    "mevotype_reviewer_points_1_2/TCGA_PRAD_GDC_450K/",
    "TCGA_PRAD_GDC_patient_clinical.tsv"
  )
)

missing_inputs <- unlist(inputs)[!file.exists(unlist(inputs))]
if (length(missing_inputs)) {
  stop("Missing input files:\n", paste(missing_inputs, collapse = "\n"))
}

# =============================================================================
# Utility functions
# =============================================================================

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

extract_icgc_patient <- function(x) {
  x <- as.character(x)
  out <- sub(".*(CRUK_PC_[0-9]{4}).*", "\\1", x)
  out[!grepl("CRUK_PC_[0-9]{4}", x)] <- NA_character_
  out
}

normalise_chr <- function(x) {
  x <- as.character(x)
  ifelse(
    grepl("^chr", x, ignore.case = TRUE),
    paste0("chr", sub("^chr", "", x, ignore.case = TRUE)),
    paste0("chr", x)
  )
}

as_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  z <- as.character(x)
  out <- suppressWarnings(as.Date(z, format = "%d/%m/%Y"))
  missing <- is.na(out) & !is.na(z) & nzchar(z)
  out[missing] <- suppressWarnings(as.Date(z[missing]))
  out
}

canonical_genomic_sample <- function(x) {
  x <- as.character(x)
  paste0("PD", sub("[^0-9].*$", "", sub("^PD", "", x)))
}

methylation_core_key <- function(x) {
  sub("_T[^_]+_", "_", as.character(x))
}

collapse_candidates <- function(x) {
  x <- sort(unique(na.omit(as.character(x))))
  if (!length(x)) NA_character_ else paste(x, collapse = "|")
}

first_candidate <- function(x) {
  x <- sort(unique(na.omit(as.character(x))))
  if (!length(x)) NA_character_ else x[1]
}

make_sample_map <- function(sample_names) {
  data.table(
    methylation_sample = sample_names,
    patient = extract_icgc_patient(sample_names),
    tissue = fcase(
      grepl("_TUM$", sample_names), "Tumour",
      grepl("_BEN$", sample_names), "Benign",
      grepl("_BPH$", sample_names), "BPH",
      default = "Other"
    )
  )
}

score_directional_burden <- function(delta, dmr_ids, directions, threshold) {
  keep <- dmr_ids %in% colnames(delta)
  ids <- dmr_ids[keep]
  directions <- directions[keep]
  if (!length(ids)) stop("No scoreable DMRs for directional burden.")
  oriented_delta <- sweep(delta[, ids, drop = FALSE], 2L, directions, "*")
  hits <- oriented_delta >= threshold
  hits[is.na(hits)] <- FALSE
  setNames(rowSums(hits) / length(ids), rownames(delta))
}

fmt_p <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 0.001) return(formatC(x, format = "e", digits = 1))
  formatC(x, format = "f", digits = 3)
}

compact_contrast <- function(x) {
  x <- gsub("ISUP low / neuronal high", "ISUP-L/N-high", x, fixed = TRUE)
  x <- gsub("ISUP low / neuronal low", "ISUP-L/N-low", x, fixed = TRUE)
  x <- gsub(" vs ISUP-L/N-low", "", x, fixed = TRUE)
  x <- gsub("MEvotype 2 / MMB low", "ME2/MMB-low", x, fixed = TRUE)
  x <- gsub("MEvotype 2 / MMB high", "ME2/MMB-high", x, fixed = TRUE)
  x <- gsub("MEvotype 1 / MMB low", "ME1/MMB-low", x, fixed = TRUE)
  x <- gsub("MEvotype 1 / MMB high", "ME1/MMB-high", x, fixed = TRUE)
  x <- gsub("Canonical / MMB low", "Can/MMB-low", x, fixed = TRUE)
  x <- gsub("Canonical / MMB high", "Can/MMB-high", x, fixed = TRUE)
  x <- gsub("Alternative / MMB low", "Alt/MMB-low", x, fixed = TRUE)
  x <- gsub("Alternative / MMB high", "Alt/MMB-high", x, fixed = TRUE)
  x <- gsub(" vs ME2/MMB-low", "", x, fixed = TRUE)
  x <- gsub(" vs Can/MMB-low", "", x, fixed = TRUE)
  x
}

# =============================================================================
# DMR definitions
# =============================================================================

neuronal_raw <- fread(inputs$neuronal_definition)[
  score == "All neuronal 30-DMR"
]
neuronal_raw[, direction := fcase(
  orientation == "Lower methylation = higher score", -1L,
  orientation == "Higher methylation = higher score", 1L,
  default = NA_integer_
)]
if (nrow(neuronal_raw) != 30L || anyNA(neuronal_raw$direction) ||
    any(sign(neuronal_raw$endpoint_change) != neuronal_raw$direction)) {
  stop("Neuronal DMR definition audit failed.")
}

mevotype_raw <- fread(inputs$mevotype_dmr_definition)
mevotype_raw[, module := fcase(
  grepl("^MEvotype1-specific", event_class), "M1",
  grepl("^MEvotype2-specific", event_class), "M2",
  default = NA_character_
)]
mevotype_raw[, direction := fcase(
  grepl("hypermethylation", event_class, fixed = TRUE), 1L,
  grepl("hypomethylation", event_class, fixed = TRUE), -1L,
  default = NA_integer_
)]
mevotype_sets <- mevotype_raw[!is.na(module), .(
  module,
  dmr_id = as.character(dmr_id),
  direction
)]
if (mevotype_sets[module == "M1", .N] != 390L ||
    mevotype_sets[module == "M2", .N] != 297L ||
    anyNA(mevotype_sets$direction) ||
    anyDuplicated(mevotype_sets[, paste(module, dmr_id)])) {
  stop("M1/M2 DMR definition audit failed.")
}

neuronal_set <- neuronal_raw[, .(
  dmr_id = as.character(dmr_id), direction
)]
m1_set <- mevotype_sets[module == "M1", .(dmr_id, direction)]
m2_set <- mevotype_sets[module == "M2", .(dmr_id, direction)]

score_modules <- function(delta) {
  out <- data.table(patient = rownames(delta))
  out[, M1_burden := as.numeric(score_directional_burden(
    delta, m1_set$dmr_id, m1_set$direction, 0.05
  ))]
  out[, M2_burden := as.numeric(score_directional_burden(
    delta, m2_set$dmr_id, m2_set$direction, 0.20
  ))]
  out[, neuronal_burden := as.numeric(score_directional_burden(
    delta, neuronal_set$dmr_id, neuronal_set$direction, 0.15
  ))]
  out[, MMB := M1_burden + M2_burden]
  out[]
}

# =============================================================================
# Historical 126-patient ICGC reference: rebuild fixed numerical cutoffs
# =============================================================================

message("Rebuilding the historical 126-patient ICGC threshold reference...")
icgc_curated <- fread(inputs$icgc_curated_dmr)
icgc_curated[, dmr_id := paste(
  normalise_chr(chr), as.integer(start), as.integer(end), sep = ":"
)]
if (nrow(icgc_curated) != 8915L || anyDuplicated(icgc_curated$dmr_id)) {
  stop("Curated ICGC DMR matrix audit failed.")
}

historical_cohort <- fread(
  inputs$icgc_analysis_cohort,
  colClasses = list(character = "patient")
)
if (nrow(historical_cohort) != 175L ||
    anyDuplicated(historical_cohort$patient)) {
  stop("Historical analysis-cohort audit failed.")
}

curated_benign_samples <- grep("_BEN$", names(icgc_curated), value = TRUE)
curated_benign_map <- unique(data.table(
  patient = extract_icgc_patient(curated_benign_samples),
  benign_sample = curated_benign_samples
), by = "patient")

historical_pairs <- merge(
  historical_cohort[, .(
    patient, methylation_sample, mevotype, evotype
  )],
  curated_benign_map,
  by = "patient",
  all = FALSE,
  sort = FALSE
)
historical_pairs <- historical_pairs[
  !is.na(evotype) & nzchar(trimws(as.character(evotype)))
]
historical_pairs <- historical_pairs[
  match(historical_cohort$patient, patient), nomatch = NULL
]
if (nrow(historical_pairs) != 126L ||
    anyDuplicated(historical_pairs$patient) ||
    !all(historical_pairs$methylation_sample %in% names(icgc_curated)) ||
    !all(historical_pairs$benign_sample %in% names(icgc_curated))) {
  stop("Historical 126-patient pair reconstruction failed.")
}

historical_tumour_columns <- historical_pairs$methylation_sample
historical_benign_columns <- historical_pairs$benign_sample
historical_tumour <- t(as.matrix(
  icgc_curated[, ..historical_tumour_columns]
))
historical_benign <- t(as.matrix(
  icgc_curated[, ..historical_benign_columns]
))
storage.mode(historical_tumour) <- "double"
storage.mode(historical_benign) <- "double"
rownames(historical_tumour) <- historical_pairs$patient
rownames(historical_benign) <- historical_pairs$patient
colnames(historical_tumour) <- icgc_curated$dmr_id
colnames(historical_benign) <- icgc_curated$dmr_id
historical_delta <- historical_tumour - historical_benign
historical_scores <- score_modules(historical_delta)

fixed_cuts <- c(
  M1 = median(historical_scores$M1_burden),
  M2 = median(historical_scores$M2_burden),
  MMB = median(historical_scores$MMB),
  neuronal = median(historical_scores$neuronal_burden)
)
expected_cuts <- c(
  M1 = 0.535897435897436,
  M2 = 0.452861952861953,
  MMB = 1.04260554260554,
  neuronal = 0.433333333333333
)
if (any(abs(fixed_cuts - expected_cuts) > 1e-12)) {
  stop("Historical cutoff reconstruction did not reproduce expected values.")
}

# =============================================================================
# Expanded 142-patient ICGC cohort: rebuild matching, scores and endpoints
# =============================================================================

message("Reconstructing the expanded 142-patient ICGC cohort...")
icgc_all <- fread(inputs$icgc_all_dmr)
icgc_all[, dmr_id := paste(
  normalise_chr(chr), as.integer(start), as.integer(end), sep = ":"
)]
if (nrow(icgc_all) != 8915L || anyDuplicated(icgc_all$dmr_id)) {
  stop("Complete ICGC DMR matrix audit failed.")
}

all_samples <- grep("^CRUK_PC_", names(icgc_all), value = TRUE)
all_sample_map <- make_sample_map(all_samples)
curated_samples <- grep("^CRUK_PC_", names(icgc_curated), value = TRUE)
curated_sample_map <- make_sample_map(curated_samples)

curated_tumour <- curated_sample_map[tissue == "Tumour", .(
  selected_tumour_sample = first_candidate(methylation_sample)
), by = patient]
curated_benign <- curated_sample_map[tissue == "Benign", .(
  selected_benign_sample = first_candidate(methylation_sample)
), by = patient]
curated_pairs <- merge(curated_tumour, curated_benign, by = "patient")
if (nrow(curated_pairs) != 163L || anyDuplicated(curated_pairs$patient)) {
  stop("Expected 163 curated ICGC tumour-benign pairs.")
}

all_tumour <- all_sample_map[tissue == "Tumour", .(
  tumour_candidates = collapse_candidates(methylation_sample),
  n_tumour_candidates = uniqueN(methylation_sample)
), by = patient]
all_benign <- all_sample_map[tissue == "Benign", .(
  benign_candidates = collapse_candidates(methylation_sample),
  n_benign_candidates = uniqueN(methylation_sample)
), by = patient]

evotype_source <- fread(inputs$icgc_evotype)
evotype_source[, genomic_core := canonical_genomic_sample(Sample)]
evotype_source[, evotype := factor(
  ifelse(Evotype == 1L, "Alternative", "Canonical"),
  levels = c("Canonical", "Alternative")
)]

genomic_map <- fread(inputs$icgc_genomic_map)
genomic_map[, genomic_core := canonical_genomic_sample(Sample_ID)]
genomic_map[, has_follow_up := as.integer(!is.na(Age_at_last_followup))]
setorder(genomic_map, genomic_core, -has_follow_up)
genomic_one <- unique(genomic_map[, .(
  genomic_core,
  patient = ICGC_ID_short,
  genomic_sample_map = Sample_ID,
  genomic_core_reported = CorePD
)], by = "genomic_core")

target159 <- merge(
  evotype_source[, .(
    genomic_core,
    genomic_sample_evotype = Sample,
    evotype
  )],
  genomic_one,
  by = "genomic_core",
  all.x = TRUE
)
if (nrow(target159) != 159L || anyNA(target159$patient) ||
    anyDuplicated(target159$patient)) {
  stop("Evotype-to-patient mapping did not produce 159 unique patients.")
}

sample_metadata <- fread(inputs$icgc_sample_metadata)
sample_metadata[, `:=`(
  patient = extract_icgc_patient(sample_names),
  mapped_genomic_core = canonical_genomic_sample(WGS_Local_Sample_ID),
  methylation_key = methylation_core_key(sample_names)
)]
exact_tumour_metadata <- sample_metadata[
  Condition == "TUM" & !is.na(patient) & !is.na(mapped_genomic_core),
  .(
    exact_tumour_metadata_samples = collapse_candidates(sample_names),
    exact_tumour_keys = collapse_candidates(methylation_key)
  ),
  by = .(patient, genomic_core = mapped_genomic_core)
]

target159 <- merge(target159, all_tumour, by = "patient", all.x = TRUE)
target159 <- merge(target159, all_benign, by = "patient", all.x = TRUE)
target159 <- merge(target159, curated_tumour, by = "patient", all.x = TRUE)
target159 <- merge(target159, curated_benign, by = "patient", all.x = TRUE)
target159 <- merge(
  target159,
  exact_tumour_metadata,
  by = c("patient", "genomic_core"),
  all.x = TRUE
)

candidate_rows <- all_sample_map[tissue == "Tumour"]
candidate_rows[, methylation_key := methylation_core_key(methylation_sample)]
exact_keys_long <- sample_metadata[
  Condition == "TUM" & !is.na(patient) & !is.na(mapped_genomic_core),
  .(
    patient,
    genomic_core = mapped_genomic_core,
    methylation_key
  )
]
candidate_rows <- merge(
  candidate_rows,
  target159[, .(patient, genomic_core)],
  by = "patient"
)
candidate_rows <- merge(
  candidate_rows,
  exact_keys_long[, .(
    patient,
    genomic_core,
    methylation_key,
    exact_genomic_core_match = TRUE
  )],
  by = c("patient", "genomic_core", "methylation_key"),
  all.x = TRUE
)
candidate_rows[
  is.na(exact_genomic_core_match), exact_genomic_core_match := FALSE
]
candidate_rows[, in_curated_350 := methylation_sample %chin% curated_samples]
setorder(
  candidate_rows,
  patient,
  -in_curated_350,
  -exact_genomic_core_match,
  methylation_sample
)
chosen_tumour <- candidate_rows[, .SD[1], by = patient][, .(
  patient,
  tumour_sample = methylation_sample,
  tumour_exact_genomic_core = exact_genomic_core_match,
  tumour_from_curated_350 = in_curated_350
)]

benign_rows <- all_sample_map[tissue == "Benign"]
benign_rows[, in_curated_350 := methylation_sample %chin% curated_samples]
setorder(benign_rows, patient, -in_curated_350, methylation_sample)
chosen_benign <- benign_rows[, .SD[1], by = patient][, .(
  patient,
  benign_sample = methylation_sample,
  benign_from_curated_350 = in_curated_350
)]

target159 <- merge(target159, chosen_tumour, by = "patient", all.x = TRUE)
target159 <- merge(target159, chosen_benign, by = "patient", all.x = TRUE)
target159[, scoreable_pair := !is.na(tumour_sample) & !is.na(benign_sample)]
target159[, historical_126 := patient %chin% historical_pairs$patient]
target159[, matching_status := fcase(
  historical_126,
  "Historical 126 curated same-patient pair",
  scoreable_pair & tumour_exact_genomic_core,
  "Additional exact genomic-core tumour with same-patient benign",
  scoreable_pair,
  "Additional same-patient alternate tumour core with same-patient benign",
  !is.na(tumour_sample) & is.na(benign_sample),
  "Same-patient tumour available; no same-patient benign",
  default = "No same-patient methylation tumour available"
)]
if (sum(target159$scoreable_pair) != 142L ||
    sum(!is.na(target159$tumour_sample)) != 147L ||
    sum(target159$historical_126) != 126L) {
  stop("Expanded same-patient matching audit failed.")
}

icgc_pairs142 <- target159[scoreable_pair == TRUE]
icgc_tumour_columns <- icgc_pairs142$tumour_sample
icgc_benign_columns <- icgc_pairs142$benign_sample
icgc_tumour142 <- t(as.matrix(icgc_all[, ..icgc_tumour_columns]))
icgc_benign142 <- t(as.matrix(icgc_all[, ..icgc_benign_columns]))
storage.mode(icgc_tumour142) <- "double"
storage.mode(icgc_benign142) <- "double"
rownames(icgc_tumour142) <- icgc_pairs142$patient
rownames(icgc_benign142) <- icgc_pairs142$patient
colnames(icgc_tumour142) <- icgc_all$dmr_id
colnames(icgc_benign142) <- icgc_all$dmr_id
icgc_delta142 <- icgc_tumour142 - icgc_benign142
icgc_scores142 <- score_modules(icgc_delta142)

clinical <- as.data.table(read_excel(inputs$icgc_clinical, sheet = "d"))
clinical[, patient := extract_icgc_patient(ICGC_ID)]
clinical <- unique(clinical[!is.na(patient)], by = "patient")
clinical[, `:=`(
  surgery_date = as_date_safe(date_of_sample_collection),
  relapse_date = as_date_safe(date_of_first_relapse),
  censor_date = as_date_safe(last_follow_up_date)
)]
clinical[, recognised_event := as.integer(
  !is.na(relapse_date) & !is.na(relapse_type) & nzchar(relapse_type) &
    tolower(relapse_type) != "na/unknown"
)]
clinical[, recognised_years := as.numeric(
  fifelse(recognised_event == 1L, relapse_date, censor_date) - surgery_date
) / 365.25]
clinical[, `:=`(
  gleason_1 = as_num(gleason_pattern_1),
  gleason_2 = as_num(gleason_pattern_2),
  age = as_num(age_at_diagnosis)
)]
clinical[, grade_group := fcase(
  gleason_1 + gleason_2 <= 6, "GG1",
  gleason_1 == 3 & gleason_2 == 4, "GG2",
  gleason_1 == 4 & gleason_2 == 3, "GG3",
  gleason_1 + gleason_2 == 8, "GG4",
  gleason_1 + gleason_2 >= 9, "GG5",
  default = NA_character_
)]
clinical[, grade_high := as.integer(
  as.integer(sub("GG", "", grade_group, fixed = TRUE)) >= 3L
)]

mevotype_source <- fread(inputs$icgc_mevotype)
mevotype_source[, patient := extract_icgc_patient(sample)]
mevotype1_patients <- unique(mevotype_source$patient)
historical_pairs[, recovered_mevotype := ifelse(
  patient %chin% mevotype1_patients, "MEvotype1", "MEvotype2"
)]
if (any(historical_pairs$recovered_mevotype != historical_pairs$mevotype)) {
  stop("MEvotype source did not reproduce the historical assignments.")
}

icgc142 <- merge(
  icgc_pairs142[, .(
    patient,
    evotype,
    tumour_sample,
    benign_sample,
    tumour_exact_genomic_core,
    tumour_from_curated_350,
    benign_from_curated_350,
    historical_126,
    matching_status
  )],
  icgc_scores142,
  by = "patient"
)
icgc142 <- merge(
  icgc142,
  clinical[, .(
    patient,
    time = recognised_years,
    event = recognised_event,
    grade_group,
    grade_high,
    age,
    relapse_type
  )],
  by = "patient",
  all.x = TRUE
)
icgc142[, mevotype := factor(
  ifelse(patient %chin% mevotype1_patients, "MEvotype 1", "MEvotype 2"),
  levels = c("MEvotype 2", "MEvotype 1")
)]
icgc142[, evotype := factor(
  as.character(evotype), levels = c("Canonical", "Alternative")
)]
if (nrow(icgc142) != 142L || anyDuplicated(icgc142$patient) ||
    sum(icgc142$event) != 37L ||
    anyNA(icgc142[, .(
      time, event, grade_high, age, mevotype, evotype,
      M1_burden, M2_burden, MMB, neuronal_burden
    )])) {
  stop("Expanded 142-patient analysis audit failed.")
}

# =============================================================================
# TCGA-PRAD: aggregate probes to the same DMRs and rebuild scores
# =============================================================================

message("Aggregating TCGA HM450 probes and rebuilding TCGA scores...")
tcga_se <- readRDS(inputs$tcga_se)
tcga_beta <- assay(tcga_se, 1)
storage.mode(tcga_beta) <- "double"

tcga_coords <- icgc_curated[, .(
  chr = normalise_chr(chr),
  start = as.integer(start),
  end = as.integer(end),
  dmr_id
)]
dmr_ranges <- GRanges(tcga_coords$chr, IRanges(tcga_coords$start, tcga_coords$end))
probe_hits <- findOverlaps(rowRanges(tcga_se), dmr_ranges, ignore.strand = TRUE)
probe_map <- unique(data.table(
  probe = queryHits(probe_hits),
  dmr = subjectHits(probe_hits)
))
split_probe <- split(probe_map$probe, probe_map$dmr)
tcga_dmr_beta <- t(vapply(split_probe, function(rows) {
  z <- colMeans(tcga_beta[rows, , drop = FALSE], na.rm = TRUE)
  z[!is.finite(z)] <- NA_real_
  z
}, numeric(ncol(tcga_beta))))
rownames(tcga_dmr_beta) <- tcga_coords$dmr_id[as.integer(names(split_probe))]
colnames(tcga_dmr_beta) <- colnames(tcga_se)

tcga_sample_meta <- data.table(
  sample = colnames(tcga_dmr_beta),
  patient = substr(colnames(tcga_dmr_beta), 1, 12),
  type_code = suppressWarnings(as.integer(substr(colnames(tcga_dmr_beta), 14, 15))),
  completeness = colMeans(is.finite(tcga_dmr_beta))
)
pick_one <- function(x) {
  setorder(x, patient, -completeness, sample)
  x[, .SD[1], by = patient]
}
tcga_tumour_meta <- pick_one(tcga_sample_meta[type_code == 1L])
tcga_normal_meta <- pick_one(tcga_sample_meta[type_code == 11L])
if (nrow(tcga_tumour_meta) < 490L || nrow(tcga_normal_meta) < 45L) {
  stop("TCGA sample selection audit failed.")
}

tcga_normal_reference <- apply(
  tcga_dmr_beta[, tcga_normal_meta$sample, drop = FALSE],
  1L,
  median,
  na.rm = TRUE
)
eligible <- is.finite(tcga_normal_reference) &
  rowMeans(is.finite(
    tcga_dmr_beta[, tcga_tumour_meta$sample, drop = FALSE]
  )) >= 0.90
eligible_ids <- rownames(tcga_dmr_beta)[eligible]
tcga_tumour_beta <- t(
  tcga_dmr_beta[eligible_ids, tcga_tumour_meta$sample, drop = FALSE]
)
rownames(tcga_tumour_beta) <- tcga_tumour_meta$patient
tcga_delta <- sweep(
  tcga_tumour_beta,
  2L,
  tcga_normal_reference[eligible_ids],
  "-"
)
tcga_scores <- score_modules(tcga_delta)

tcga_survival <- fread(inputs$tcga_survival)
tcga_survival[, patient := substr(as.character(sample), 1, 12)]
tcga_survival <- tcga_survival[
  `cancer type abbreviation` == "PRAD" & PFI %in% 0:1 &
    is.finite(PFI.time) & PFI.time > 0,
  .(
    patient,
    event = as.integer(PFI),
    PFI_days = as.numeric(PFI.time)
  )
]
setorder(tcga_survival, patient, -PFI_days)
tcga_survival <- tcga_survival[, .SD[1], by = patient]
tcga_survival[, time := PFI_days / 365.25]

tcga_clinical <- fread(inputs$tcga_clinical)
tcga_clinical[, `:=`(
  patient = as.character(submitter_id),
  primary_pattern = as_num(sub("Pattern ", "", primary_gleason_grade)),
  age = as_num(age_at_diagnosis) / 365.25
)]
tcga_clinical[, grade_high := fcase(
  as_num(gleason_score) >= 8, 1L,
  as_num(gleason_score) == 7 & primary_pattern == 4, 1L,
  as_num(gleason_score) <= 6, 0L,
  as_num(gleason_score) == 7 & primary_pattern == 3, 0L,
  default = NA_integer_
)]
tcga_clinical <- unique(
  tcga_clinical[, .(patient, age, grade_high)],
  by = "patient"
)

tcga <- merge(
  tcga_scores,
  tcga_survival[, .(patient, time, event)],
  by = "patient"
)
tcga <- merge(tcga, tcga_clinical, by = "patient", all.x = TRUE)
if (nrow(tcga) != 498L || anyDuplicated(tcga$patient) ||
    sum(tcga$event) != 93L || anyNA(tcga$grade_high) ||
    anyNA(tcga[, .(MMB, neuronal_burden)])) {
  stop("TCGA-PRAD analysis audit failed.")
}

# =============================================================================
# Fixed groups transferred from the historical ICGC 126-patient cohort
# =============================================================================

icgc142[, `:=`(
  M1_high = M1_burden >= fixed_cuts["M1"],
  M2_high = M2_burden >= fixed_cuts["M2"],
  MMB_high = MMB >= fixed_cuts["MMB"],
  neuronal_high = neuronal_burden >= fixed_cuts["neuronal"]
)]
tcga[, `:=`(
  MMB_high = MMB >= fixed_cuts["MMB"],
  neuronal_high = neuronal_burden >= fixed_cuts["neuronal"]
)]

make_isup_neuronal_group <- function(grade_high, neuronal_high) {
  factor(
    fcase(
      grade_high == 1L, "ISUP high",
      grade_high == 0L & neuronal_high, "ISUP low / neuronal high",
      grade_high == 0L & !neuronal_high, "ISUP low / neuronal low",
      default = NA_character_
    ),
    levels = c(
      "ISUP low / neuronal low",
      "ISUP low / neuronal high",
      "ISUP high"
    )
  )
}
icgc142[, isup_neuronal_group := make_isup_neuronal_group(
  grade_high, neuronal_high
)]
tcga[, isup_neuronal_group := make_isup_neuronal_group(
  grade_high, neuronal_high
)]

icgc142[, mevo_mmb_group := factor(
  paste(
    as.character(mevotype),
    ifelse(MMB_high, "MMB high", "MMB low"),
    sep = " / "
  ),
  levels = c(
    "MEvotype 2 / MMB low",
    "MEvotype 2 / MMB high",
    "MEvotype 1 / MMB low",
    "MEvotype 1 / MMB high"
  )
)]
icgc142[, evo_mmb_group := factor(
  paste(
    as.character(evotype),
    ifelse(MMB_high, "MMB high", "MMB low"),
    sep = " / "
  ),
  levels = c(
    "Canonical / MMB low",
    "Canonical / MMB high",
    "Alternative / MMB low",
    "Alternative / MMB high"
  )
)]

# =============================================================================
# Survival plotting and statistics
# =============================================================================

blue <- "#4C78A8"
red <- "#D62728"
gold <- "#E69F00"
green <- "#008837"
purple <- "#7B3294"

make_km_panel <- function(d, group, reference, palette, panel,
                          x_label, y_label, adjust_vars,
                          annotation_size = 1.9,
                          base_size = 8.5,
                          legend_alias = NULL,
                          legend_counts = TRUE,
                          legend_rows = 1L) {
  dd <- copy(d)
  dd[, group__ := relevel(factor(group), ref = reference)]
  dd <- dd[
    is.finite(time) & time >= 0 & event %in% 0:1 & !is.na(group__)
  ]
  if (nlevels(droplevels(dd$group__)) < 2L) {
    stop("Fewer than two groups in panel ", panel)
  }

  fit <- survfit(Surv(time, event) ~ group__, data = dd)
  strata <- rep(sub("^group__=", "", names(fit$strata)), fit$strata)
  curve <- data.table(
    time = fit$time,
    survival = fit$surv,
    n_censor = fit$n.censor,
    group = strata
  )
  group_levels <- levels(dd$group__)
  curve <- rbind(
    data.table(time = 0, survival = 1, n_censor = 0,
               group = group_levels),
    curve
  )
  curve[, group := factor(group, levels = group_levels)]

  counts <- dd[, .(n = .N, events = sum(event)), by = group__]
  displayed_group <- as.character(counts$group__)
  if (!is.null(legend_alias)) {
    displayed_group <- unname(legend_alias[displayed_group])
  }
  legend_labels <- setNames(
    if (legend_counts) paste0(displayed_group, " (n=", counts$n, ")")
    else displayed_group,
    counts$group__
  )

  unadjusted <- coxph(Surv(time, event) ~ group__, data = dd, x = TRUE)
  unadjusted_summary <- summary(unadjusted)
  terms <- rownames(unadjusted_summary$coefficients)
  adjusted_formula <- as.formula(paste0(
    "Surv(time, event) ~ group__ + ", paste(adjust_vars, collapse = " + ")
  ))
  adjusted <- coxph(adjusted_formula, data = dd, x = TRUE)
  adjusted_summary <- summary(adjusted)

  stats <- rbindlist(lapply(terms, function(term) {
    b <- unname(coef(unadjusted)[term])
    se <- unname(sqrt(vcov(unadjusted)[term, term]))
    contrast <- sub("^group__", "", term)
    data.table(
      panel,
      contrast = paste0(contrast, " vs ", reference),
      n = nrow(dd),
      events = sum(dd$event),
      HR = exp(b),
      CI_low = exp(b - 1.96 * se),
      CI_high = exp(b + 1.96 * se),
      Cox_p = unname(
        unadjusted_summary$coefficients[term, "Pr(>|z|)"]
      ),
      overall_Cox_p = unname(unadjusted_summary$logtest["pvalue"]),
      adjusted_for = paste(adjust_vars, collapse = "+"),
      adjusted_p = unname(
        adjusted_summary$coefficients[term, "Pr(>|z|)"]
      ),
      adjusted_model_n = adjusted$n,
      Harrell_C_index = unname(unadjusted_summary$concordance[1])
    )
  }))

  hr_lines <- vapply(seq_len(nrow(stats)), function(i) {
    z <- stats[i]
    paste0(
      compact_contrast(z$contrast), ": HR=", sprintf("%.2f", z$HR),
      " (", sprintf("%.2f", z$CI_low), "–", sprintf("%.2f", z$CI_high),
      "); p=", fmt_p(z$Cox_p), "; padj=", fmt_p(z$adjusted_p)
    )
  }, character(1))
  extra <- c(
    if (nrow(stats) > 1L) {
      paste0("Overall Cox p=", fmt_p(stats$overall_Cox_p[1]))
    } else {
      character()
    },
    paste0("C-index=", sprintf("%.3f", stats$Harrell_C_index[1]))
  )
  annotation <- paste(c(hr_lines, extra), collapse = "\n")

  plot <- ggplot(curve, aes(time, survival, colour = group)) +
    geom_step(linewidth = 0.78) +
    geom_point(
      data = curve[n_censor > 0],
      shape = 3,
      size = 0.85,
      stroke = 0.45
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_colour_manual(
      values = palette,
      labels = legend_labels,
      drop = FALSE
    ) +
    guides(colour = guide_legend(nrow = legend_rows, byrow = TRUE)) +
    labs(
      tag = panel,
      x = x_label,
      y = y_label,
      colour = NULL
    ) +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(face = "bold", size = base_size + 2),
      plot.tag.position = "topleft",
      legend.position = "bottom",
      legend.text = element_text(size = base_size - 1.3),
      legend.key.width = unit(0.9, "lines"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )

  list(plot = plot, stats = stats, counts = counts)
}

binary_group <- function(high) {
  factor(ifelse(high, "High", "Low"), levels = c("Low", "High"))
}
binary_palette <- c("Low" = blue, "High" = red)
signal_legend_alias <- function(signal) {
  c("Low" = paste(signal, "low"), "High" = paste(signal, "high"))
}
hierarchy_palette <- c(
  "ISUP low / neuronal low" = blue,
  "ISUP low / neuronal high" = gold,
  "ISUP high" = red
)
hierarchy_legend_alias <- c(
  "ISUP low / neuronal low" = "ISUP low / N-low",
  "ISUP low / neuronal high" = "ISUP low / N-high",
  "ISUP high" = "ISUP high"
)

# Source survival analyses. The requested submission layout is assembled below.
figure1_panels <- list(
  A = make_km_panel(
    icgc142, binary_group(icgc142$M1_high), "Low", binary_palette, "A",
    "Years from surgery", "Recurrence-free survival",
    c("grade_high", "age"),
    legend_alias = signal_legend_alias("M1")
  ),
  B = make_km_panel(
    icgc142, binary_group(icgc142$M2_high), "Low", binary_palette, "B",
    "Years from surgery", "Recurrence-free survival",
    c("grade_high", "age"),
    legend_alias = signal_legend_alias("M2")
  ),
  C = make_km_panel(
    icgc142, binary_group(icgc142$MMB_high), "Low", binary_palette, "C",
    "Years from surgery", "Recurrence-free survival",
    c("grade_high", "age"),
    legend_alias = signal_legend_alias("MMB")
  ),
  D = make_km_panel(
    tcga, binary_group(tcga$MMB_high), "Low", binary_palette, "D",
    "Years from diagnosis", "Progression-free interval",
    c("grade_high", "age"),
    legend_alias = signal_legend_alias("MMB")
  ),
  E = make_km_panel(
    icgc142, binary_group(icgc142$neuronal_high), "Low", binary_palette, "E",
    "Years from surgery", "Recurrence-free survival",
    c("grade_high", "age"),
    legend_alias = signal_legend_alias("Neuronal")
  ),
  F = make_km_panel(
    icgc142, icgc142$isup_neuronal_group,
    "ISUP low / neuronal low", hierarchy_palette, "F",
    "Years from surgery", "Recurrence-free survival", "age",
    annotation_size = 1.55,
    legend_alias = hierarchy_legend_alias,
    legend_rows = 2L
  ),
  G = make_km_panel(
    tcga, binary_group(tcga$neuronal_high), "Low", binary_palette, "G",
    "Years from diagnosis", "Progression-free interval",
    c("grade_high", "age"),
    legend_alias = signal_legend_alias("Neuronal")
  ),
  H = make_km_panel(
    tcga, tcga$isup_neuronal_group,
    "ISUP low / neuronal low", hierarchy_palette, "H",
    "Years from diagnosis", "Progression-free interval", "age",
    annotation_size = 1.55,
    legend_alias = hierarchy_legend_alias,
    legend_rows = 2L
  )
)



# =============================================================================
# Additional categorical survival and neuronal hyper/hypo analyses
# =============================================================================

figure2_survival <- list(
  A = make_km_panel(
    icgc142,
    icgc142$mevotype,
    "MEvotype 2",
    c("MEvotype 2" = green, "MEvotype 1" = purple),
    "A",
    "Years from surgery",
    "Recurrence-free survival",
    c("grade_high", "age"),
    annotation_size = 2.5,
    base_size = 10
  ),
  B = make_km_panel(
    icgc142,
    icgc142$evotype,
    "Canonical",
    c("Canonical" = blue, "Alternative" = red),
    "B",
    "Years from surgery",
    "Recurrence-free survival",
    c("grade_high", "age"),
    annotation_size = 2.5,
    base_size = 10
  )
)

# Hyper/hypo status is based on the sign of the mean tumour-minus-benign delta
# beta across all 142 patients, independent of MEvotype or Evotype labels.
neuronal_ids <- intersect(neuronal_set$dmr_id, colnames(icgc_delta142))
neuronal_cohort_mean <- colMeans(
  icgc_delta142[, neuronal_ids, drop = FALSE], na.rm = TRUE
)
if (any(!is.finite(neuronal_cohort_mean)) ||
    any(neuronal_cohort_mean == 0)) {
  stop("Unable to classify neuronal DMRs as hyper or hypo.")
}
neuronal_hyper_ids <- names(neuronal_cohort_mean)[neuronal_cohort_mean > 0]
neuronal_hypo_ids <- names(neuronal_cohort_mean)[neuronal_cohort_mean < 0]
if (!length(neuronal_hyper_ids) || !length(neuronal_hypo_ids)) {
  stop("Neuronal hyper/hypo split is empty.")
}

neuronal_patient_values <- data.table(
  patient = rownames(icgc_delta142),
  neuronal_hyper_mean_delta_beta = rowMeans(
    icgc_delta142[, neuronal_hyper_ids, drop = FALSE], na.rm = TRUE
  ),
  neuronal_hypo_mean_delta_beta = rowMeans(
    icgc_delta142[, neuronal_hypo_ids, drop = FALSE], na.rm = TRUE
  )
)
neuronal_patient_values <- merge(
  neuronal_patient_values,
  icgc142[, .(
    patient, mevotype, evotype, MMB_high,
    mevo_mmb_group, evo_mmb_group
  )],
  by = "patient"
)
if (nrow(neuronal_patient_values) != 142L ||
    anyNA(neuronal_patient_values)) {
  stop("Neuronal boxplot patient-data audit failed.")
}

make_box_panel <- function(d, group_column, value_column, palette, panel,
                           y_label) {
  dd <- copy(d)
  dd[, group__ := factor(get(group_column))]
  dd[, value__ := get(value_column)]
  counts <- dd[, .N, by = group__]
  labels <- setNames(
    paste0(counts$group__, "\n(n=", counts$N, ")"), counts$group__
  )
  n_groups <- nlevels(droplevels(dd$group__))
  if (n_groups == 2L) {
    test_result <- wilcox.test(value__ ~ group__, data = dd, exact = FALSE)
    test_label <- "Wilcoxon"
  } else {
    test_result <- kruskal.test(value__ ~ group__, data = dd)
    test_label <- "Kruskal–Wallis"
  }
  summary_stats <- dd[, .(
    n = .N,
    median = median(value__),
    Q1 = quantile(value__, 0.25),
    Q3 = quantile(value__, 0.75),
    mean = mean(value__),
    SD = sd(value__)
  ), by = group__]
  summary_stats[, `:=`(
    panel = panel,
    classification = group_column,
    neuronal_direction = value_column,
    test = test_label,
    omnibus_p = test_result$p.value,
    Wilcoxon_p = if (n_groups == 2L) test_result$p.value else NA_real_,
    Kruskal_Wallis_p = if (n_groups > 2L) test_result$p.value else NA_real_
  )]

  plot <- ggplot(dd, aes(group__, value__, fill = group__)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
    geom_boxplot(
      width = 0.58,
      outlier.shape = NA,
      alpha = 0.78,
      linewidth = 0.55
    ) +
    geom_jitter(
      width = 0.13,
      height = 0,
      size = 1.05,
      alpha = 0.5,
      shape = 21,
      stroke = 0.15
    ) +
    annotate(
      "text",
      x = (n_groups + 1) / 2,
      y = Inf,
      label = paste0(test_label, " p=", fmt_p(test_result$p.value)),
      vjust = 1.3,
      size = 2.8
    ) +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
    labs(tag = panel, x = NULL, y = y_label, fill = NULL) +
    theme_bw(base_size = 8.5) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(face = "bold", size = 11),
      plot.tag.position = "topleft",
      axis.text.x = element_text(size = 7.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )

  list(plot = plot, stats = summary_stats)
}

mevotype_palette <- c("MEvotype 2" = green, "MEvotype 1" = purple)
evotype_palette <- c("Canonical" = blue, "Alternative" = red)
figure2_boxplots <- list(
  C = make_box_panel(
    neuronal_patient_values,
    "mevotype",
    "neuronal_hyper_mean_delta_beta",
    mevotype_palette,
    "C",
    "Neuronal hyper mean Δβ"
  ),
  D = make_box_panel(
    neuronal_patient_values,
    "mevotype",
    "neuronal_hypo_mean_delta_beta",
    mevotype_palette,
    "D",
    "Neuronal hypo mean Δβ"
  ),
  E = make_box_panel(
    neuronal_patient_values,
    "evotype",
    "neuronal_hyper_mean_delta_beta",
    evotype_palette,
    "E",
    "Neuronal hyper mean Δβ"
  ),
  F = make_box_panel(
    neuronal_patient_values,
    "evotype",
    "neuronal_hypo_mean_delta_beta",
    evotype_palette,
    "F",
    "Neuronal hypo mean Δβ"
  )
)



# =============================================================================
# Joint MEvotype/Evotype and fixed MMB analyses
# =============================================================================

mevo_mmb_palette <- c(
  "MEvotype 2 / MMB low" = "#80C89D",
  "MEvotype 2 / MMB high" = green,
  "MEvotype 1 / MMB low" = "#B68BC1",
  "MEvotype 1 / MMB high" = purple
)
evo_mmb_palette <- c(
  "Canonical / MMB low" = "#8FB3D1",
  "Canonical / MMB high" = blue,
  "Alternative / MMB low" = "#F08A86",
  "Alternative / MMB high" = red
)
mevo_mmb_legend_alias <- c(
  "MEvotype 2 / MMB low" = "ME2 low",
  "MEvotype 2 / MMB high" = "ME2 high",
  "MEvotype 1 / MMB low" = "ME1 low",
  "MEvotype 1 / MMB high" = "ME1 high"
)
evo_mmb_legend_alias <- c(
  "Canonical / MMB low" = "Canonical / MMB low",
  "Canonical / MMB high" = "Canonical / MMB high",
  "Alternative / MMB low" = "Alternative / MMB low",
  "Alternative / MMB high" = "Alternative / MMB high"
)

figure3_survival <- list(
  A = make_km_panel(
    icgc142,
    icgc142$mevo_mmb_group,
    "MEvotype 2 / MMB low",
    mevo_mmb_palette,
    "A",
    "Years from surgery",
    "Recurrence-free survival",
    c("grade_high", "age"),
    annotation_size = 1.8,
    base_size = 10,
    legend_alias = mevo_mmb_legend_alias,
    legend_counts = FALSE,
    legend_rows = 2L
  ),
  B = make_km_panel(
    icgc142,
    icgc142$evo_mmb_group,
    "Canonical / MMB low",
    evo_mmb_palette,
    "B",
    "Years from surgery",
    "Recurrence-free survival",
    c("grade_high", "age"),
    annotation_size = 1.8,
    base_size = 10,
    legend_alias = evo_mmb_legend_alias,
    legend_rows = 2L
  )
)

neuronal_patient_values[, mevo_mmb_box_group := factor(
  as.character(mevo_mmb_group),
  levels = names(mevo_mmb_palette),
  labels = c(
    "ME2 / MMB low", "ME2 / MMB high",
    "ME1 / MMB low", "ME1 / MMB high"
  )
)]
neuronal_patient_values[, evo_mmb_box_group := factor(
  as.character(evo_mmb_group),
  levels = names(evo_mmb_palette),
  labels = c(
    "Can. / MMB low", "Can. / MMB high",
    "Alt. / MMB low", "Alt. / MMB high"
  )
)]
mevo_mmb_box_palette <- setNames(
  unname(mevo_mmb_palette), levels(neuronal_patient_values$mevo_mmb_box_group)
)
evo_mmb_box_palette <- setNames(
  unname(evo_mmb_palette), levels(neuronal_patient_values$evo_mmb_box_group)
)

figure3_boxplots <- list(
  C = make_box_panel(
    neuronal_patient_values,
    "mevo_mmb_box_group",
    "neuronal_hyper_mean_delta_beta",
    mevo_mmb_box_palette,
    "C",
    "Neuronal hyper mean Δβ"
  ),
  D = make_box_panel(
    neuronal_patient_values,
    "mevo_mmb_box_group",
    "neuronal_hypo_mean_delta_beta",
    mevo_mmb_box_palette,
    "D",
    "Neuronal hypo mean Δβ"
  ),
  E = make_box_panel(
    neuronal_patient_values,
    "evo_mmb_box_group",
    "neuronal_hyper_mean_delta_beta",
    evo_mmb_box_palette,
    "E",
    "Neuronal hyper mean Δβ"
  ),
  F = make_box_panel(
    neuronal_patient_values,
    "evo_mmb_box_group",
    "neuronal_hypo_mean_delta_beta",
    evo_mmb_box_palette,
    "F",
    "Neuronal hypo mean Δβ"
  )
)



# =============================================================================
# Clean, auditable result exports
# =============================================================================

figure2_boxplot_stats <- rbindlist(
  lapply(figure2_boxplots, `[[`, "stats"), fill = TRUE
)
figure2_boxplot_stats[, analysis_set := "MEvotype/Evotype"]
figure3_boxplot_stats <- rbindlist(
  lapply(figure3_boxplots, `[[`, "stats"), fill = TRUE
)
figure3_boxplot_stats[, analysis_set := "Joint classification/MMB"]
boxplot_stats <- rbindlist(
  list(figure2_boxplot_stats, figure3_boxplot_stats),
  fill = TRUE
)
setnames(boxplot_stats, "group__", "group")
setcolorder(boxplot_stats, c(
  "analysis_set", "panel", "classification", "neuronal_direction", "group",
  setdiff(names(boxplot_stats), c(
    "analysis_set", "panel", "classification", "neuronal_direction", "group"
  ))
))
fwrite(boxplot_stats, file.path(table_dir, "additional_neuronal_boxplot_statistics.csv"))

# Diagnostics for whether an unweighted M1 + M2 score is statistically
# defensible. Standardised coefficients are reported per one cohort SD.
icgc142[, `:=`(
  z_M1 = as.numeric(scale(M1_burden)),
  z_M2 = as.numeric(scale(M2_burden)),
  z_MMB = as.numeric(scale(MMB))
)]
mmb_diagnostic_models <- list(
  `M1 only, unadjusted` = coxph(
    Surv(time, event) ~ z_M1, data = icgc142, x = TRUE
  ),
  `M2 only, unadjusted` = coxph(
    Surv(time, event) ~ z_M2, data = icgc142, x = TRUE
  ),
  `MMB only, unadjusted` = coxph(
    Surv(time, event) ~ z_MMB, data = icgc142, x = TRUE
  ),
  `M1 and M2, unadjusted` = coxph(
    Surv(time, event) ~ z_M1 + z_M2, data = icgc142, x = TRUE
  ),
  `MMB, grade and age` = coxph(
    Surv(time, event) ~ z_MMB + grade_high + age,
    data = icgc142,
    x = TRUE
  ),
  `M1, M2, grade and age` = coxph(
    Surv(time, event) ~ z_M1 + z_M2 + grade_high + age,
    data = icgc142,
    x = TRUE
  )
)
mmb_model_diagnostics <- rbindlist(lapply(
  names(mmb_diagnostic_models),
  function(model_name) {
    fit <- mmb_diagnostic_models[[model_name]]
    s <- summary(fit)
    coefficients <- as.data.table(s$coefficients, keep.rownames = "term")
    confidence <- as.data.table(s$conf.int, keep.rownames = "term")
    merge(
      coefficients[, .(
        term,
        coefficient = coef,
        HR = `exp(coef)`,
        Wald_p = `Pr(>|z|)`
      )],
      confidence[, .(
        term,
        CI_low = `lower .95`,
        CI_high = `upper .95`
      )],
      by = "term"
    )[, `:=`(
      model = model_name,
      n = fit$n,
      events = fit$nevent,
      AIC = AIC(fit),
      Harrell_C_index = unname(s$concordance[1])
    )]
  }
), fill = TRUE)
setcolorder(mmb_model_diagnostics, c(
  "model", "term", "n", "events", "coefficient", "HR", "CI_low",
  "CI_high", "Wald_p", "AIC", "Harrell_C_index"
))
fwrite(
  mmb_model_diagnostics,
  file.path(table_dir, "MMB_continuous_model_diagnostics.csv")
)

raw_separate_fit <- coxph(
  Surv(time, event) ~ M1_burden + M2_burden,
  data = icgc142,
  x = TRUE
)
raw_equal_weight_fit <- coxph(
  Surv(time, event) ~ I(M1_burden + M2_burden),
  data = icgc142,
  x = TRUE
)
equal_weight_lrt <- 2 * (
  as.numeric(logLik(raw_separate_fit)) -
    as.numeric(logLik(raw_equal_weight_fit))
)
mmb_summary_diagnostics <- data.table(
  n = nrow(icgc142),
  events = sum(icgc142$event),
  Pearson_M1_M2 = cor(icgc142$M1_burden, icgc142$M2_burden),
  Spearman_M1_M2 = cor(
    icgc142$M1_burden, icgc142$M2_burden, method = "spearman"
  ),
  equal_weight_constraint_LR_chisq = equal_weight_lrt,
  equal_weight_constraint_df = 1L,
  equal_weight_constraint_p = pchisq(
    equal_weight_lrt, df = 1L, lower.tail = FALSE
  ),
  MMB_proportional_hazards_global_p = cox.zph(
    mmb_diagnostic_models[["MMB only, unadjusted"]]
  )$table["GLOBAL", "p"]
)
fwrite(
  mmb_summary_diagnostics,
  file.path(table_dir, "MMB_score_summary_diagnostics.csv")
)

mmb_component_counts <- icgc142[, .(
  n = .N,
  events = sum(event),
  event_fraction = mean(event),
  median_MMB = median(MMB)
), by = .(
  M1_group = ifelse(M1_high, "M1 high", "M1 low"),
  M2_group = ifelse(M2_high, "M2 high", "M2 low"),
  MMB_group = ifelse(MMB_high, "MMB high", "MMB low")
)]
setorder(mmb_component_counts, MMB_group, M1_group, M2_group)
fwrite(
  mmb_component_counts,
  file.path(table_dir, "MMB_component_group_counts.csv")
)

within_lineage_mmb <- rbindlist(lapply(
  c("mevotype", "evotype"),
  function(classification) {
    rbindlist(lapply(levels(icgc142[[classification]]), function(class_value) {
      dd <- icgc142[get(classification) == class_value]
      unadjusted <- coxph(
        Surv(time, event) ~ MMB_high,
        data = dd,
        x = TRUE
      )
      adjusted <- coxph(
        Surv(time, event) ~ MMB_high + grade_high + age,
        data = dd,
        x = TRUE
      )
      b <- unname(coef(unadjusted)["MMB_highTRUE"])
      se <- sqrt(vcov(unadjusted)["MMB_highTRUE", "MMB_highTRUE"])
      data.table(
        classification,
        class = class_value,
        n = nrow(dd),
        events = sum(dd$event),
        MMB_low_n = sum(!dd$MMB_high),
        MMB_high_n = sum(dd$MMB_high),
        MMB_high_vs_low_HR = exp(b),
        CI_low = exp(b - 1.96 * se),
        CI_high = exp(b + 1.96 * se),
        Cox_p = summary(unadjusted)$coefficients[
          "MMB_highTRUE", "Pr(>|z|)"
        ],
        adjusted_for = "grade_high+age",
        adjusted_p = summary(adjusted)$coefficients[
          "MMB_highTRUE", "Pr(>|z|)"
        ]
      )
    }))
  }
))
fwrite(
  within_lineage_mmb,
  file.path(table_dir, "MMB_within_lineage_survival.csv")
)

fwrite(
  data.table(
    signal = names(fixed_cuts),
    threshold = c(0.05, 0.20, NA_real_, 0.15),
    historical_126_median_cutoff = as.numeric(fixed_cuts),
    rule = c(
      "Fraction of directional M1 DMRs meeting oriented delta beta >= 0.05",
      "Fraction of directional M2 DMRs meeting oriented delta beta >= 0.20",
      "M1 burden + M2 burden",
      "Fraction of directional neuronal DMRs meeting oriented delta beta >= 0.15"
    )
  ),
  file.path(table_dir, "fixed_thresholds.csv")
)

fwrite(
  icgc142[, .(
    patient,
    tumour_sample,
    benign_sample,
    matching_status,
    historical_126,
    mevotype = as.character(mevotype),
    evotype = as.character(evotype),
    M1_burden,
    M2_burden,
    MMB,
    neuronal_burden,
    M1_high,
    M2_high,
    MMB_high,
    neuronal_high,
    grade_group,
    grade_high,
    age,
    recognised_recurrence_years = time,
    recognised_recurrence_event = event,
    relapse_type,
    isup_neuronal_group,
    mevo_mmb_group,
    evo_mmb_group
  )],
  file.path(table_dir, "ICGC_142_patient_analysis.csv")
)

fwrite(
  tcga[, .(
    patient,
    M1_burden,
    M2_burden,
    MMB,
    neuronal_burden,
    MMB_high,
    neuronal_high,
    grade_high,
    age,
    PFI_years = time,
    PFI_event = event,
    isup_neuronal_group
  )],
  file.path(table_dir, "TCGA_patient_analysis.csv")
)

fwrite(
  neuronal_patient_values,
  file.path(table_dir, "ICGC_142_neuronal_hyper_hypo_patient_values.csv")
)
fwrite(
  data.table(
    dmr_id = names(neuronal_cohort_mean),
    mean_tumour_minus_benign_delta_beta = as.numeric(neuronal_cohort_mean),
    direction = ifelse(neuronal_cohort_mean > 0, "Hyper", "Hypo")
  ),
  file.path(table_dir, "ICGC_142_neuronal_DMR_hyper_hypo_classification.csv")
)
fwrite(
  target159[, .(
    patient,
    evotype = as.character(evotype),
    genomic_sample_evotype,
    genomic_sample_map,
    tumour_sample,
    benign_sample,
    tumour_exact_genomic_core,
    tumour_from_curated_350,
    benign_from_curated_350,
    historical_126,
    scoreable_pair,
    matching_status
  )],
  file.path(table_dir, "ICGC_159_matching_audit.csv")
)

cohort_audit <- data.table(
  cohort = c("Historical ICGC reference", "Expanded ICGC", "TCGA-PRAD"),
  n = c(126L, nrow(icgc142), nrow(tcga)),
  events = c(NA_integer_, sum(icgc142$event), sum(tcga$event)),
  endpoint = c(
    "Threshold calibration only",
    "Recognised postoperative recurrence",
    "PanCancer progression-free interval"
  ),
  missing_grade = c(NA_integer_, sum(is.na(icgc142$grade_high)),
                    sum(is.na(tcga$grade_high))),
  missing_age = c(NA_integer_, sum(is.na(icgc142$age)),
                  sum(is.na(tcga$age))),
  tumour_samples_selected = c(126L, 142L, nrow(tcga_tumour_meta)),
  benign_or_normal_reference_samples = c(
    126L, 142L, nrow(tcga_normal_meta)
  ),
  eligible_DMRs = c(8915L, 8915L, length(eligible_ids))
)
fwrite(cohort_audit, file.path(table_dir, "cohort_audit.csv"))

coverage_audit <- data.table(
  cohort = rep(c("ICGC 142", "TCGA-PRAD"), each = 3),
  module = rep(c("M1", "M2", "Neuronal"), 2),
  source_DMRs = rep(c(nrow(m1_set), nrow(m2_set), nrow(neuronal_set)), 2),
  scoreable_DMRs = c(
    sum(m1_set$dmr_id %in% colnames(icgc_delta142)),
    sum(m2_set$dmr_id %in% colnames(icgc_delta142)),
    sum(neuronal_set$dmr_id %in% colnames(icgc_delta142)),
    sum(m1_set$dmr_id %in% eligible_ids),
    sum(m2_set$dmr_id %in% eligible_ids),
    sum(neuronal_set$dmr_id %in% eligible_ids)
  )
)
coverage_audit[, fraction_scoreable := scoreable_DMRs / source_DMRs]
fwrite(coverage_audit, file.path(table_dir, "DMR_coverage_audit.csv"))

source_manifest <- data.table(
  source = names(inputs),
  path = unlist(inputs, use.names = FALSE),
  bytes = as.numeric(file.info(unlist(inputs))$size)
)
fwrite(source_manifest, file.path(table_dir, "source_manifest.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(table_dir, "R_session_info.txt")
)

message(
  "Final figures complete. ICGC n=", nrow(icgc142),
  " (", sum(icgc142$event), " events); TCGA n=", nrow(tcga),
  " (", sum(tcga$event), " events); neuronal hyper/hypo DMRs=",
  length(neuronal_hyper_ids), "/", length(neuronal_hypo_ids), "."
)


# =============================================================================
# Submission figures: A B C / D E F and A B
# All plots below use the same raw-data-derived cohorts, cutoffs and models.
# Statistical annotations and all other figure text are 8 pt at export size.
# =============================================================================

panel_statistics_label <- function(stats) {
  if (nrow(stats) == 1L) {
    contrast <- stats$contrast[1]
    contrast <- sub("MEvotype 1 vs MEvotype 2", "MEvotype 1 vs 2",
                    contrast, fixed = TRUE)
    return(paste(
      contrast,
      sprintf("HR=%.2f (%.2f–%.2f)", stats$HR[1],
              stats$CI_low[1], stats$CI_high[1]),
      sprintf("p=%s; padj=%s", fmt_p(stats$Cox_p[1]),
              fmt_p(stats$adjusted_p[1])),
      sprintf("C-index=%.3f", stats$Harrell_C_index[1]),
      sep = "\n"
    ))
  }
  if (nrow(stats) != 3L) stop("Unexpected number of Cox contrasts.")
  group_names <- c(
    "MEvotype 2 / MMB high" = "ME2/H",
    "MEvotype 1 / MMB low" = "ME1/L",
    "MEvotype 1 / MMB high" = "ME1/H"
  )
  comparison_groups <- sub(" vs .*", "", stats$contrast)
  short_names <- unname(group_names[comparison_groups])
  if (anyNA(short_names)) stop("Unexpected MEvotype/MMB contrast.")
  contrast_lines <- unlist(lapply(seq_len(nrow(stats)), function(i) {
    c(
      sprintf("%s HR=%.2f (%.2f–%.2f)", short_names[i], stats$HR[i],
              stats$CI_low[i], stats$CI_high[i]),
      sprintf("p=%s; padj=%s", fmt_p(stats$Cox_p[i]),
              fmt_p(stats$adjusted_p[i]))
    )
  }))
  paste(c(
    contrast_lines,
    sprintf("Overall p=%s; C=%.3f", fmt_p(stats$overall_Cox_p[1]),
            stats$Harrell_C_index[1])
  ), collapse = "\n")
}

clean_panel <- function(x, letter, heading) {
  p <- x$plot +
    labs(title = heading, tag = tolower(letter)) +
    annotate("text", x = Inf, y = 0.018,
             label = panel_statistics_label(x$stats),
             hjust = 1.02, vjust = 0, size = 8 / 2.845276,
             family = "Helvetica", lineheight = 0.85) +
    theme(
      text = element_text(family = "Helvetica", size = 8),
      plot.title = element_text(size = 8, face = "bold", hjust = 0,
                                margin = margin(b = 1.5, unit = "mm")),
      plot.tag = element_text(family = "Helvetica", size = 14,
                              face = "bold"),
      plot.tag.position = "topleft",
      plot.tag.location = "margin",
      plot.title.position = "panel",
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.key.height = unit(3.1, "mm"),
      legend.key.width = unit(5, "mm"),
      legend.spacing.x = unit(1, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 0.25, colour = "black"),
      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )
  if (nrow(x$counts) == 2L) {
    p <- p + guides(colour = guide_legend(ncol = 1))
  }
  p
}

requested <- list(
  Figure_1 = list(
    A = list(source = figure2_survival$A, heading = "MEvotype · ICGC",
             analysis = "MEvotype", cohort = "ICGC 142",
             endpoint = "Recognised postoperative recurrence",
             stem = "Figure_1_A_MEvotype_ICGC"),
    B = list(source = figure1_panels$A, heading = "M1 burden · ICGC",
             analysis = "M1 burden", cohort = "ICGC 142",
             endpoint = "Recognised postoperative recurrence",
             stem = "Figure_1_B_M1_burden_ICGC"),
    C = list(source = figure1_panels$B, heading = "M2 burden · ICGC",
             analysis = "M2 burden", cohort = "ICGC 142",
             endpoint = "Recognised postoperative recurrence",
             stem = "Figure_1_C_M2_burden_ICGC"),
    D = list(source = figure1_panels$C, heading = "MMB · ICGC",
             analysis = "MMB", cohort = "ICGC 142",
             endpoint = "Recognised postoperative recurrence",
             stem = "Figure_1_D_MMB_burden_ICGC"),
    E = list(source = figure3_survival$A, heading = "MEvotype / MMB · ICGC",
             analysis = "MEvotype/MMB joint groups", cohort = "ICGC 142",
             endpoint = "Recognised postoperative recurrence",
             stem = "Figure_1_E_MEvotype_MMB_ICGC"),
    F = list(source = figure1_panels$D, heading = "MMB · TCGA",
             analysis = "MMB", cohort = "TCGA-PRAD",
             endpoint = "PanCancer progression-free interval",
             stem = "Figure_1_F_MMB_TCGA")
  ),
  Figure_2 = list(
    A = list(source = figure1_panels$E, heading = "Neuronal burden · ICGC",
             analysis = "Neuronal burden", cohort = "ICGC 142",
             endpoint = "Recognised postoperative recurrence",
             stem = "Figure_2_A_neuronal_burden_ICGC"),
    B = list(source = figure1_panels$G, heading = "Neuronal burden · TCGA",
             analysis = "Neuronal burden", cohort = "TCGA-PRAD",
             endpoint = "PanCancer progression-free interval",
             stem = "Figure_2_B_neuronal_burden_TCGA")
  )
)

save_pair <- function(plot, dir, stem, width_mm, height_mm) {
  ggsave(file.path(dir, paste0(stem, ".pdf")), plot,
         device = cairo_pdf, width = width_mm, height = height_mm,
         units = "mm", bg = "white", limitsize = FALSE)
  ggsave(file.path(dir, paste0(stem, ".png")), plot,
         width = width_mm, height = height_mm, units = "mm",
         dpi = 300, bg = "white", limitsize = FALSE)
}

figure_dimensions <- data.table(
  figure = character(), panel = character(), file_stem = character(),
  width_mm = numeric(), height_mm = numeric(), dpi_png = integer()
)
requested_stats <- list()
requested_counts <- list()
for (fig in names(requested)) {
  panel_width <- if (fig == "Figure_1") 60 else 90
  panel_height <- if (fig == "Figure_1") 85 else 75
  panel_plots <- list()
  for (letter in names(requested[[fig]])) {
    spec <- requested[[fig]][[letter]]
    p <- clean_panel(spec$source, letter, spec$heading)
    panel_plots[[letter]] <- p
    save_pair(p, panel_dir, spec$stem, panel_width, panel_height)
    figure_dimensions <- rbind(
      figure_dimensions,
      data.table(figure = fig, panel = letter, file_stem = spec$stem,
                 width_mm = panel_width, height_mm = panel_height,
                 dpi_png = 300L)
    )
    st <- copy(spec$source$stats)
    st[, `:=`(figure = fig, panel = letter, analysis = spec$analysis,
              cohort = spec$cohort, endpoint = spec$endpoint)]
    requested_stats[[paste(fig, letter)]] <- st
    ct <- copy(spec$source$counts)
    setnames(ct, "group__", "group")
    setorder(ct, group)
    ct[, `:=`(figure = fig, panel = letter, analysis = spec$analysis,
              cohort = spec$cohort, endpoint = spec$endpoint)]
    requested_counts[[paste(fig, letter)]] <- ct
  }
  combined <- wrap_plots(panel_plots, ncol = if (fig == "Figure_1") 3 else 2)
  width_mm <- 180
  height_mm <- if (fig == "Figure_1") 170 else 75
  save_pair(combined, figure_dir, fig, width_mm, height_mm)
  figure_dimensions <- rbind(
    figure_dimensions,
    data.table(figure = fig, panel = "combined", file_stem = fig,
               width_mm = width_mm, height_mm = height_mm, dpi_png = 300L)
  )
}

requested_stats <- rbindlist(requested_stats, fill = TRUE)
setcolorder(requested_stats, c(
  "figure", "panel", "analysis", "cohort", "endpoint",
  setdiff(names(requested_stats), c(
    "figure", "panel", "analysis", "cohort", "endpoint"
  ))
))
fwrite(requested_stats, file.path(table_dir, "survival_statistics.csv"))
requested_counts <- rbindlist(requested_counts, fill = TRUE)
setcolorder(requested_counts,
            c("figure", "panel", "analysis", "cohort", "endpoint",
              "group", "n", "events"))
fwrite(requested_counts, file.path(table_dir, "survival_group_counts.csv"))
fwrite(figure_dimensions, file.path(table_dir, "figure_dimensions.csv"))

count_text <- function(fig, letter) {
  d <- requested_counts[figure == fig & panel == letter]
  paste(sprintf("%s, n=%d", d$group, d$n), collapse = "; ")
}
legend_lines <- c(
  "# Figure legends",
  "",
  "## Figure 1. MEvotype and methylation burden in ICGC and TCGA",
  "",
  paste0(
    "Kaplan–Meier estimates are arranged a–c above d–f. ",
    "Panels a–e show recurrence-free survival in the 142-patient ICGC ",
    "same-patient tumour–benign cohort; panel f shows progression-free ",
    "interval in 498 TCGA-PRAD tumours. ",
    "(a) MEvotype (", count_text("Figure_1", "A"), "). ",
    "(b) M1 burden (", count_text("Figure_1", "B"), "). ",
    "(c) M2 burden (", count_text("Figure_1", "C"), "). ",
    "(d) MMB, the sum of M1 and M2 burdens (",
    count_text("Figure_1", "D"), "). ",
    "(e) Joint MEvotype/MMB groups (",
    count_text("Figure_1", "E"), "). ",
    "(f) TCGA MMB (", count_text("Figure_1", "F"), "). ",
    "In panel e, ME1/ME2 denote MEvotype 1/2; L/H in the ",
    "statistical annotation denotes MMB-low/high, and ",
    "MEvotype 2/MMB-low is the hazard-ratio reference. ",
    "High/low status uses the fixed historical 126-patient ICGC medians: ",
    "M1 0.5358974, M2 0.4528620 and MMB 1.0426055; the same ",
    "numerical MMB cutoff is transferred to TCGA. Crosses mark censoring. ",
    "Panels show Cox hazard ratios, 95% confidence intervals, two-sided ",
    "p-values and Harrell’s C-index. The value labelled padj is the Cox ",
    "p-value after adjustment for ISUP grade group and age, not a ",
    "multiple-testing correction. Full estimates are in ",
    "tables/survival_statistics.csv."
  ),
  "",
  "## Figure 2. Neuronal methylation burden in ICGC and TCGA",
  "",
  paste0(
    "Kaplan–Meier estimates for (a) recurrence-free survival in the ",
    "142-patient ICGC cohort (", count_text("Figure_2", "A"), ") ",
    "and (b) progression-free interval in 498 TCGA-PRAD tumours (",
    count_text("Figure_2", "B"), "). ",
    "Neuronal burden is dichotomised at the fixed historical 126-patient ",
    "ICGC median of 0.4333333 in both cohorts. Crosses mark censoring. ",
    "Panels show Cox hazard ratios, 95% confidence intervals, two-sided ",
    "p-values and Harrell’s C-index. The value labelled padj is the Cox ",
    "p-value after adjustment for ISUP grade group and age, not a ",
    "multiple-testing correction. Full estimates are in ",
    "tables/survival_statistics.csv. The two cohorts use different ",
    "endpoints and methylation reference designs, as described in METHODS.md."
  ),
  ""
)
writeLines(legend_lines, file.path(output_root, "FIGURE_LEGENDS.md"),
           useBytes = TRUE)
message("Submission figures and standalone panels exported to ", figure_dir)
