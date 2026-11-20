# analyze reprocessed results
# detect co-occurring aas

library(tidyverse)
library(Biostrings)
library(stringr)
library(data.table)
library(janitor)

`%notin%` <- function(a, b) !(a %in% b)

setwd("/home/lam.wayn/aas_stabilization/reprocessing")
results_dir <- "/scratch/lam.wayn/MS_SEARCH/results"

# ============================================================
# 1. Helper functions
# ============================================================

# classify an entry name as wt / single / multi -- single canonical
# implementation, used both when building feature_df and when reading
# abundance files, so the two can never silently disagree
get_type <- function(entry_names) {
  muts <- sub("^.*?_HUMAN_", "", entry_names)
  dplyr::case_when(
    is.na(entry_names)                 ~ NA_character_,
    grepl("_REFERENCE$", entry_names)  ~ "wt",
    stringr::str_count(muts, "_") >= 1 ~ "multi",
    TRUE                               ~ "single"
  )
}

# extract a patient/case ID from a TMT sample column name (first token before "_")
# NOTE: computed later and stored on long_abundance, but not yet consumed
# by any downstream step in this script -- kept here in case a future
# per-patient / repeated-measures analysis is added
get_patient <- function(sample_names) {
  vapply(strsplit(as.character(sample_names), "_", fixed = TRUE),
         function(x) x[1], character(1))
}

# check if a given aa is at the given site within the peptide
check_mutation <- function(peptide, aa, site, start, end) {
  in_range <- !is.na(site) & site >= start & site <= end
  rel_pos  <- site - start + 1
  res <- in_range & substring(peptide, rel_pos, rel_pos) == aa
  res[is.na(res)] <- FALSE
  res
}

# check all mutations of one protein, then rebuild its name with only the ones that pass
rebuild_entry <- function(entry_name, peptide, start, end) {
  prefix    <- str_match(entry_name, "^(.*?_HUMAN)_")[, 2]
  mutations <- str_split(str_match(entry_name, "_HUMAN_(.*)$")[, 2], "_")[[1]]
  parsed    <- str_match(mutations, "^([A-Za-z*])(\\d+)([A-Za-z*])$")
  # peptides come from mutated sequences, so check for the mutant residue
  keep <- check_mutation(peptide, parsed[, 4], as.integer(parsed[, 3]), start, end)
  if (any(keep)) paste(c(prefix, mutations[keep]), collapse = "_") else NA_character_
}

# if no mutation is in the peptide, fall back to the same gene's REFERENCE entry
# when it is listed in Mapped Proteins; otherwise NA
reference_entry <- function(entry_name, mapped_proteins) {
  ref    <- paste0(str_match(entry_name, "^(.*?_HUMAN)_")[, 2], "_REFERENCE")
  mapped <- sub("^.*\\|", "", str_trim(str_split(mapped_proteins, ",")[[1]]))
  if (ref %in% mapped) ref else NA_character_
}

# ============================================================
# 2. Build feature_df from the FASTA (canonical gene/entry/type catalog)
# ============================================================
# NOTE: this must run BEFORE section 3 (verification), which reads
# feature_names.csv and expects a gene/accession/entry_name/type schema.
# Previously section 3 ran first and either errored on a fresh run or
# silently used a stale file.
fasta_path <- "/scratch/lam.wayn/FASTA/multi_singleton_ref_supp7.fasta"
fasta_seqs <- readAAStringSet(fasta_path)
features <- names(fasta_seqs)
m <- str_match(features, "^[^|]+\\|([^|]+)\\|(\\S+)")
g <- str_match(features, "GN=(\\S+)")[, 2]

feature_df <- data.frame(
  gene       = g,
  accession  = m[, 2],   # Q02790
  entry_name = m[, 3]    # FKBP4_HUMAN_E21N_E301C
)
feature_df$type <- get_type(feature_df$entry_name)

write.csv(feature_df, "feature_names.csv", row.names = FALSE)
feature_df <- read.csv("feature_names.csv", header = TRUE)

# ============================================================
# 3. Verify every multi has its corresponding singles in the FASTA
# ============================================================
multi_substitutions <- feature_df %>%
  dplyr::filter(type == "multi") %>%
  # split "GENE_HUMAN_mut1_mut2..." into gene and the mutation string
  dplyr::mutate(
    gene     = stringr::str_match(entry_name, "^(.*?)_HUMAN_")[, 2],
    mut_list = stringr::str_match(entry_name, "_HUMAN_(.*)$")[, 2]
  ) %>%
  # one row per mutation
  tidyr::separate_rows(mut_list, sep = "_") %>%
  dplyr::filter(!is.na(mut_list), mut_list != "") %>%
  # parse e.g. N83R -> from = N, pos = 83, to = R
  tidyr::extract(mut_list, into = c("from", "pos", "to"),
                 regex = "^([A-Za-z*])(\\d+)([A-Za-z*])$",
                 remove = FALSE, convert = TRUE) %>%
  dplyr::rename(mutation = mut_list) %>%
  dplyr::select(gene, mutation, from, pos, to, dplyr::everything(), -entry_name)

single_substitutions <- feature_df %>%
  dplyr::filter(type == "single") %>%
  dplyr::mutate(
    gene     = stringr::str_match(entry_name, "^(.*?)_HUMAN_")[, 2],
    mut_list = stringr::str_match(entry_name, "_HUMAN_(.*)$")[, 2]
  ) %>%
  tidyr::separate_rows(mut_list, sep = "_") %>%
  dplyr::filter(!is.na(mut_list), mut_list != "") %>%
  tidyr::extract(mut_list, into = c("from", "pos", "to"),
                 regex = "^([A-Za-z*])(\\d+)([A-Za-z*])$",
                 remove = FALSE, convert = TRUE) %>%
  dplyr::rename(mutation = mut_list) %>%
  dplyr::select(gene, mutation, from, pos, to, dplyr::everything(), -entry_name)

missing <- anti_join(multi_substitutions, single_substitutions,
                      by = join_by("gene", "mutation", "accession"))
if (nrow(missing) > 0) {
  message("fasta warning: some multi constituent singles are missing from the FASTA")
}

# ============================================================
# 4. Plex lists and output check
# ============================================================
plex_list <- list.files(results_dir)
ucec  <- plex_list[grepl("UCEC", plex_list)]
ccrcc <- plex_list[grepl("CCRCC", plex_list)]
brca  <- plex_list[grepl("BCprospective", plex_list)]
luad  <- plex_list[grepl("LUAD", plex_list)]
lscc  <- plex_list[grepl("LSCC", plex_list)]
pdac  <- plex_list[grepl("PDA", plex_list)]
length(c(ucec, ccrcc, brca, luad, lscc, pdac))

# NOTE: defined once here; previously redefined identically later in the
# "statistical test" section
experiment_list <- list(ucec = ucec, ccrcc = ccrcc, brca = brca,
                        luad = luad, lscc = lscc, pdac = pdac)
all_plexes <- unlist(experiment_list)

failed <- list()
for (plex in plex_list) {
  if (!file.exists(file.path(results_dir, plex, "tmt-report", "abundance_protein_MD.tsv"))) {
    failed <- append(failed, plex)
  }
}
if (length(failed) > 0) message("Do not proceed")

# did we recover all singles from experiment?

# ============================================================
# 5. Shared constants used by every abundance-reading step below
# ============================================================
meta_cols <- c("Entry Name", "Index", "Gene", "ProteinID", "Peptide", "SequenceWindow",
               "Start", "End", "MaxPepProb", "Spectrum Number", "Protein",
               "Protein Description", "Mapped Genes", "Mapped Proteins",
               "ReferenceIntensity", "type")
keys <- c("updated_entry_name", "SequenceWindow")

# only genes with at least one declared multi are needed for the statistical test
relevant_genes <- feature_df %>% dplyr::filter(type == "multi") %>% dplyr::pull(gene) %>% unique()

# ============================================================
# 6. Aggregate peptide-level abundance across all plexes (summary stats)
# ============================================================
plex_summaries <- list()

for (i in seq_along(experiment_list)) {
  experiment_name <- names(experiment_list)[i]
  experiment      <- experiment_list[[i]]

  for (j in seq_along(experiment)) {
    plex        <- experiment[[j]]
    plex_number <- substr(plex, 1, 2)

    abundance <- fread(file.path(results_dir, plex, "tmt-report", "abundance_peptide_MD.tsv"),
                       sep = "\t")

    # type computed directly with the canonical get_type() -- no merge
    # against a separate type_labels table needed
    abundance[, type := get_type(`Entry Name`)]
    abundance[, updated_entry_name := `Entry Name`]

    # wt keeps its name; singles and multis keep only mutations found in the peptide,
    # or become the gene's REFERENCE (if mapped) when none are found, else NA
    checks <- unique(abundance[type != "wt",
                               .(`Entry Name`, Peptide, Start, End, `Mapped Proteins`)])
    if (nrow(checks) > 0) {
      checks[, new_name := unlist(Map(rebuild_entry, `Entry Name`, Peptide, Start, End))]
      checks[is.na(new_name),
             new_name := unlist(Map(reference_entry, `Entry Name`, `Mapped Proteins`))]
      abundance[checks, on = c("Entry Name", "Peptide", "Start", "End", "Mapped Proteins"),
                updated_entry_name := i.new_name]
    }

    sample_cols <- setdiff(names(abundance), c(meta_cols, keys))

    # peptide info: sums and counts so plexes can be combined later
    pep_info <- abundance[, .(
      Index          = first(Index),
      Gene           = first(Gene),
      ProteinID      = first(ProteinID),
      Peptide        = first(Peptide),
      Start          = first(Start),
      End            = first(End),
      pep_sum        = sum(MaxPepProb,        na.rm = TRUE),
      pep_n          = sum(!is.na(MaxPepProb)),
      spectrum_sum   = sum(`Spectrum Number`, na.rm = TRUE),
      spectrum_n     = sum(!is.na(`Spectrum Number`))
    ), by = keys]

    # abundances: one value per sample per sequence window (averages rows that collapsed to
    # the same window), then sum and count across samples, ignoring missing values
    per_sample <- abundance[, lapply(.SD, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)),
                            by = keys, .SDcols = sample_cols]
    vals <- as.matrix(per_sample[, ..sample_cols])
    per_sample <- per_sample[, ..keys][
      , `:=`(abundance_sum = rowSums(vals, na.rm = TRUE),
             n_samples     = as.integer(rowSums(!is.na(vals))))]

    plex_summaries[[plex]] <-
      merge(pep_info, per_sample, by = keys)[, `:=`(experiment = experiment_name, plex = plex_number)]
  }
}

# combine plex summaries into one row per unique entry name + sequence window
all_updated <- rbindlist(plex_summaries, use.names = TRUE)
write.csv(all_updated, "all_updated.csv", row.names = FALSE)

# ============================================================
# 7. Summary-level analyses (re-entry point if resuming from disk)
# ============================================================
all_updated <- read.csv("all_updated.csv", header = TRUE)
feature_df  <- read.csv("feature_names.csv", header = TRUE)

all_updated <- data.table(all_updated)
all_updated[, type := get_type(updated_entry_name)]

# function to aggregate at various levels
average_at <- function(dt, by) {
  data.table(dt)[, .(
    avg_pep       = sum(pep_sum)       / sum(pep_n),
    avg_spectrum  = sum(spectrum_sum)  / sum(spectrum_n),
    avg_abundance = sum(abundance_sum) / sum(n_samples),
    n_samples     = sum(n_samples),
    n_experiments = uniqueN(experiment),
    n_plexes      = .N
  ), by = c("updated_entry_name")]
}

overall <- average_at(all_updated, by = NULL)
overall[, type := get_type(updated_entry_name)]
by_peptide <- overall[type == "multi"]

by_experiment <- average_at(all_updated, by = "experiment")
by_experiment[, type := get_type(updated_entry_name)]
by_experiment[type == "multi", uniqueN(updated_entry_name)]
by_experiment[type == "multi", table(n_experiments)]

# overall abundance of multis by experiment
selected_proteins <- c("1433B_HUMAN_V31A_H36V", "VDAC2_HUMAN_E95M_D89M")

# one point per protein, pooled over all its sequence windows and plexes
multis <- all_updated[type == "multi"]
protein_level <- multis[, .(
  avg_abundance = sum(abundance_sum) / sum(n_samples),
  total_samples = sum(n_samples)
), by = .(experiment, updated_entry_name)]
overall_point <- multis[, .(
  updated_entry_name = "overall",
  avg_abundance = sum(abundance_sum) / sum(n_samples),
  total_samples = sum(n_samples)
), by = .(experiment)]
protein_level <- protein_level[updated_entry_name %in% selected_proteins]
plotting <- rbind(protein_level, overall_point)

ggplot(plotting, aes(x = updated_entry_name, y = experiment)) +
  geom_point(aes(colour = avg_abundance, size = log10(total_samples))) +
  scale_colour_viridis_c(name = expression("Avg abundance (log"[2]*")")) +
  scale_size_continuous(name = expression("log"[10]*" total samples"), range = c(1.5, 6)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "grey92"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# 8. Compare multi/single abundance to wt, by experiment
# ============================================================
entry_level <- all_updated[, .(
  abundance_sum = sum(abundance_sum), # in log2
  n_samples     = sum(n_samples)
), by = .(Gene, updated_entry_name)]
entry_level[, type := get_type(updated_entry_name)]
entry_level[, avg_abundance := abundance_sum / n_samples]

# wt abundance per gene
wt <- entry_level[type == "wt", .(
  wt_abundance = avg_abundance, # in log2
  wt_samples   = sum(n_samples)
), by = .(Gene)]

# ratio of each single/multi to its gene's wt (abundances are log2, so subtract)
variant_ratios <- merge(entry_level[type %in% c("single", "multi")], wt, by = "Gene")
variant_ratios[, log2_ratio := avg_abundance - wt_abundance] # in log2

# mean ratio per gene, per type
gene_ratios <- variant_ratios[, .(
  mean_log2_ratio = mean(log2_ratio),   # >0 -> SAAP>BP
  mean_ratio      = 2^mean(log2_ratio), # >1 -> SAPP>BP
  n_variants      = .N
), by = .(Gene, type)]

# one row per gene, with single and multi side by side
gene_ratios_wide <- dcast(gene_ratios, Gene ~ type,
                          value.var = c("mean_log2_ratio", "mean_ratio", "n_variants"))

mean(gene_ratios_wide$mean_ratio_multi, na.rm = TRUE)
mean(gene_ratios_wide$mean_ratio_single, na.rm = TRUE)

all_updated[type == "multi", uniqueN(Peptide)]

# ============================================================
# 9. Build per-sample long-format abundances (needed for the paired test)
# ============================================================
# takes a plex name, returns a data.table with Gene, updated_entry_name,
# sample_col, abundance, plex -- one row per entry observed in one sample
build_long_abundance <- function(plex) {
  abundance <- fread(file.path(results_dir, plex, "tmt-report", "abundance_peptide_MD.tsv"), sep = "\t")
  abundance <- abundance[Gene %in% relevant_genes]
  if (nrow(abundance) == 0) return(NULL)

  abundance[, type := get_type(`Entry Name`)]
  abundance[, updated_entry_name := `Entry Name`]

  checks <- unique(abundance[type != "wt", .(`Entry Name`, Peptide, Start, End, `Mapped Proteins`)])
  if (nrow(checks) > 0) {
    checks[, new_name := unlist(Map(rebuild_entry, `Entry Name`, Peptide, Start, End))]
    checks[is.na(new_name), new_name := unlist(Map(reference_entry, `Entry Name`, `Mapped Proteins`))]
    abundance[checks, on = c("Entry Name", "Peptide", "Start", "End", "Mapped Proteins"),
              updated_entry_name := i.new_name]
  }

  sample_cols <- setdiff(names(abundance), c(meta_cols, keys))
  per_sample <- abundance[, lapply(.SD, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)),
                          by = c(keys, "Gene"), .SDcols = sample_cols]

  # pivots sample cols to sample rows
  long <- melt(per_sample, id.vars = c(keys, "Gene"), measure.vars = sample_cols,
               variable.name = "sample_col", value.name = "abundance")
  # drop missing abundances
  long <- long[!is.na(abundance)]

  # SequenceWindow = site-centered view, not a distinct peptide -> average
  # across windows sharing the same updated_entry_name (site-level views
  # of the same molecule), rather than summing them
  long <- long[, .(abundance = mean(abundance, na.rm = TRUE)),
               by = .(Gene, updated_entry_name, sample_col)]
  long[, plex := plex]
  long
}

long_list <- lapply(all_plexes, build_long_abundance)
long_abundance <- rbindlist(long_list, use.names = TRUE)

# filter out references/QC samples
long_abundance <- long_abundance[!grepl("Not_Reported|QC|Korean|WU", sample_col)]

write.csv(long_abundance, "long_abundance.csv", row.names = FALSE)
long_abundance <- as.data.table(read.csv("long_abundance.csv", header = TRUE))
long_abundance[, type := get_type(updated_entry_name)]
setcolorder(long_abundance, c("Gene", "updated_entry_name", "type"))

long_abundance[, patient := get_patient(sample_col)]

# ============================================================
# 10. Observed feature list
# ============================================================
# NOTE: this must come BEFORE section 11, since multi_parts there is
# built from feature_cols. Previously this whole block lived at the very
# end of the script, ~90 lines after multi_parts already tried to use it.
feature_cols <- long_abundance[type != "wt" & Gene %in% relevant_genes,
                               .(Gene, updated_entry_name, type)]
feature_cols <- unique(feature_cols)
setorder(feature_cols, Gene, type)

feature_cols <- merge(feature_cols, unique(as.data.table(feature_df)[, .(gene, accession)]),
                      by.x = "Gene", by.y = "gene", all.x = TRUE)
feature_cols <- feature_cols[, .(Gene, accession, updated_entry_name, type)]
write.csv(feature_cols, "observed_features.csv", row.names = FALSE)
feature_cols <- as.data.table(read.csv("observed_features.csv", header = TRUE))

# ============================================================
# 11. Component-presence categories per multi, across ALL samples
# ============================================================
# Built from presence matrix 
# input: long_abundance with type, feature_cols

# --- 1. inputs and the full sample list -----------------------------------
long_abundance <- as.data.table(read.csv("long_abundance.csv", header = TRUE))
long_abundance[, type := get_type(updated_entry_name)]

feature_df <- read.csv("feature_names.csv", header = TRUE)
relevant_genes <- feature_df %>% dplyr::filter(type == "multi") %>% dplyr::pull(gene) %>% unique()
feature_cols <- long_abundance[type != "wt" & Gene %in% relevant_genes,
                               .(Gene, updated_entry_name, type)]

samples <- sort(unique(as.character(long_abundance$sample_col)))
# long abundance has per peptide (repeats of entry names in same sample)
# unique entry name, sample col pairs (1257025)
obs <- unique(long_abundance[!is.na(updated_entry_name),
                             .(entry = updated_entry_name, sample = as.character(sample_col))])

# --- 2. map singles to multis ---------------------------------------------
gene_of <- unique(feature_cols[type == "multi", .(Gene, multi_name = updated_entry_name)])
multis  <- sort(gene_of$multi_name)
# all possible constituents
cons <- rbindlist(lapply(multis, function(m) {
  prefix <- sub("^(.*?_HUMAN_).*$", "\\1", m, perl = TRUE)
  muts   <- strsplit(sub("^.*?_HUMAN_", "", m, perl = TRUE), "_", fixed = TRUE)[[1]]
  data.table(multi = m, single = unique(paste0(prefix, muts)))
}))
k <- cons[, .N, by = multi][match(multis, multi), N]   # singles per multi, in multis order
stopifnot(all(k >= 2))

# presence binary matrix: entries x samples
P <- unclass(table(obs$entry, factor(obs$sample, levels = samples))) > 0

# get P for given entry names
rows_of <- function(nms) {
  out <- P[match(nms, rownames(P)), , drop = FALSE]
  out[is.na(out)] <- FALSE
  rownames(out) <- nms 
  out
}
never_seen <- setdiff(cons$single, rownames(P))
if (length(never_seen) > 0){
  message(length(never_seen), " constituent single(s) never observed; counted as absent in every sample")
}
# presence of multi
multi_present <- rows_of(multis)
# number of singles present in sample for multi
# rows_of(cons$single) * 1 returns P of singles as 1/0
# rowsum grouped by cons$multi adds rows with same multi (based on df order)
# [multis] orders by multis list
n_present     <- rowsum(rows_of(cons$single) * 1L, group = cons$multi)[multis, , drop = FALSE]
n_missing     <- k - n_present          # k recycles down the rows: row i uses k[i]
K             <- matrix(k, nrow = length(multis), ncol = length(samples))

# categorize every (multi, sample)
# codes follow the order of cat_levels; code 1 (all_components) = testable in section 12
cat_levels <- c("all_components", "missing_single", "missing_singles", "missing_all_singles",
                "missing_multi", "missing_multi_and_single", "missing_all")
cat_code <- ifelse(multi_present,
                   ifelse(n_missing == 0, 1L, ifelse(n_missing == K, 4L, ifelse(n_missing == 1, 2L, 3L))),
                   ifelse(n_missing == 0, 5L, ifelse(n_missing == K, 7L, 6L)))

# --- 6. count samples per category per multi -------------------------------
counts <- t(apply(cat_code, 1, tabulate, nbins = 7L))
colnames(counts) <- cat_levels

component_category_tab <- merge(gene_of, data.table(multi_name = multis, counts), by = "multi_name")
component_category_tab[, n_samples := rowSums(.SD), .SDcols = cat_levels]
setcolorder(component_category_tab, c("Gene", "multi_name", "n_samples", cat_levels))
setorder(component_category_tab, Gene, multi_name)
stopifnot(all(component_category_tab$n_samples == length(samples)))

write.csv(component_category_tab, "stattest_viability.csv", row.names = FALSE)
component_category_tab <- read.csv("stattest_viability.csv", header = T)
component_category_tab <- component_category_tab %>% mutate(proportion = all_components/n_samples)
sum(component_category_tab$all_components == 0)/length(component_category_tab$all_components)
sum(component_category_tab$all_components >577)/length(component_category_tab$all_components) # >50% of samples
# ============================================================
# 12. Statistical test: observed vs. expected multi probability
# ============================================================
# Does each multi show up more (or less) often than we'd expect if its
#   single substitutions occurred independently of each other?
#
#   1. Score every substitution as a probability against its gene's wt:
#        prob = abundance_sub / (abundance_wt + abundance_sub)
#   2. Observed = prob of the multi itself
#      Expected = product of the probs of its singles (independence)
#   3. diff = observed - expected      (> 0: multi more common than expected)
#   4. Samples missing the multi, any single, or the wt are NA and dropped.
#   5. Paired t-test of observed vs expected across the remaining samples,
#      then BH correction across multis.
#
# INPUTS FROM SECTION 11
#   samples, multis, cons (multi -> single map), k (singles per multi),
#   P (presence matrix), cat_code / cat_levels (sample categories), gene_of
#
# OUTPUTS
#   observed_vs_expected.csv  one row per usable (multi, sample)
#   stattest_results.csv      one row per multi

# --- settings ---------------------------------------------------------------
min_n <- 2      # minimum usable samples needed to run the test on a multi
alpha <- 0.05   # significance cutoff on the BH-adjusted p-value


# --- step 1: raw abundance matrix (entries x samples) -----------------------
# Same rows and columns as P. If an entry has more than one log2 value in a
# sample (e.g. the same sample column in two plexes), average them first.
abund_long <- long_abundance[
  !is.na(updated_entry_name) & !is.na(abundance),
  .(log2_abundance = mean(abundance)),
  by = .(entry = updated_entry_name, sample = as.character(sample_col))
]

abund_raw <- matrix(NA_real_, nrow = nrow(P), ncol = ncol(P), dimnames = dimnames(P))
cell <- cbind(match(abund_long$entry, rownames(P)), match(abund_long$sample, colnames(P)))
abund_raw[cell] <- 2^abund_long$log2_abundance     # log2 -> raw

# sanity check: abundances exist exactly where P says the entry was detected
stopifnot(all(!is.na(abund_raw) == P))

# look up raw abundances for a vector of entry names
# (one row per name; entries never observed come back as all-NA rows)
get_abundance <- function(entry_names) {
  out <- abund_raw[match(entry_names, rownames(abund_raw)), , drop = FALSE]
  rownames(out) <- entry_names
  out
}

# probability of a substitution relative to wt
sub_prob <- function(sub_abund, wt_abund) sub_abund / (wt_abund + sub_abund)


# --- step 2: wt abundance for each multi ------------------------------------
# The wt is the gene's REFERENCE entry, e.g. FKBP4_HUMAN_E21N_E301C ->
# FKBP4_HUMAN_REFERENCE (same rule as reference_entry() in section 1).
wt_names <- sub("^(.*?_HUMAN)_.*$", "\\1_REFERENCE", multis, perl = TRUE)
wt_abund <- get_abundance(wt_names)                  # multis x samples

never_seen_wt <- unique(wt_names[rowSums(!is.na(wt_abund)) == 0])
if (length(never_seen_wt) > 0) {
  message(length(never_seen_wt), " REFERENCE entr(ies) never observed; ",
          "their multis cannot be tested")
}


# --- step 3: observed probability (the multi) -------------------------------
multi_abund   <- get_abundance(multis)               # multis x samples
observed_prob <- sub_prob(multi_abund, wt_abund)

# --- step 4: expected probability (product of the singles) ------------------
# One row per (multi, single) pair in cons; each single is scored against
# its multi's wt.
single_abund <- get_abundance(cons$single)
single_wt    <- wt_abund[match(cons$multi, multis), , drop = FALSE]
single_prob  <- sub_prob(single_abund, single_wt)    # nrow(cons) x samples

# multiply each multi's singles together, sample by sample.
# prod() gives NA if any single is missing, so expected is NA too.
expected_prob <- matrix(NA_real_, nrow = length(multis), ncol = length(samples),
                        dimnames = list(multis, samples))
for (m in multis) {
  rows_for_m <- which(cons$multi == m)
  expected_prob[m, ] <- apply(single_prob[rows_for_m, , drop = FALSE], 2, prod)
}


# --- step 5: difference, and blank out unusable samples ---------------------
diff_prob <- observed_prob - expected_prob

# a sample is usable only if the multi AND all of its singles were detected
# ("all_components" in section 11). A missing wt already gives NA above.
not_usable <- cat_code != match("all_components", cat_levels)
observed_prob[not_usable] <- NA_real_
expected_prob[not_usable] <- NA_real_
diff_prob[not_usable]     <- NA_real_


# --- step 6: one table over ALL samples, then keep the usable ones ----------
# as.vector() reads a matrix column by column, so multis vary fastest
obs_exp_all <- data.table(
  multi_name    = rep(multis,  times = length(samples)),
  sample_col    = rep(samples, each  = length(multis)),
  category      = cat_levels[as.vector(cat_code)],
  observed_prob = as.vector(observed_prob),
  expected_prob = as.vector(expected_prob),
  diff          = as.vector(diff_prob)
)
obs_exp_all <- merge(gene_of, obs_exp_all, by = "multi_name")
setcolorder(obs_exp_all, c("Gene", "multi_name", "sample_col"))
stopifnot(nrow(obs_exp_all) == length(multis) * length(samples))

paired <- obs_exp_all[!is.na(diff)]
write.csv(paired, "observed_vs_expected.csv", row.names = FALSE)

# --- step 7: paired t-test per multi ----------------------------------------
run_test <- function(observed, expected) {
  n <- length(observed)
  p <- NA_real_
  if (n >= min_n) {
    # t.test() errors if every difference is identical; keep p = NA then
    p <- tryCatch(t.test(observed, expected, paired = TRUE)$p.value,
                  error = function(e) NA_real_)
  }
  list(mean_observed_prob = mean(observed),
       mean_expected_prob = mean(expected),
       mean_diff          = mean(observed - expected),
       sd_diff            = if (n > 1) sd(observed - expected) else NA_real_,
       p                  = p)
}

test_results <- paired[, run_test(observed_prob, expected_prob), by = .(Gene, multi_name)]


# --- step 8: sample counts per multi ----------------------------------------
sample_summary <- data.table(
  multi_name       = multis,
  n_singles        = k,
  n_total          = length(samples),
  n_all_components = rowSums(!not_usable),          # multi + all singles detected
  n_viable         = rowSums(!is.na(diff_prob))     # ... and wt detected
)
sample_summary[, n_missing_wt := n_all_components - n_viable]
sample_summary[, prop_viable  := n_viable / n_total]


# --- step 9: combine, adjust p-values, and call each multi ------------------
results <- merge(gene_of, sample_summary, by = "multi_name")
results <- merge(results, test_results, by = c("Gene", "multi_name"), all.x = TRUE)

# BH correction; multis that weren't tested (p = NA) are not counted
results[, p_adj := p.adjust(p, method = "BH")]

results[, call := fcase(
  n_viable == 0,                  "no viable samples",
  n_viable < min_n,               "insufficient samples",
  is.na(p),                       "not testable",
  p_adj < alpha & mean_diff > 0,  "enriched",
  p_adj < alpha & mean_diff < 0,  "depleted",
  default =                       "not significant"
)]
results <- results %>% select(-c("n_total", "n_viable", "n_missing_wt", "prop_viable"))

setcolorder(results, c("Gene", "multi_name"))
setorder(results, p_adj, na.last = TRUE)
write.csv(results, "stattest_results.csv", row.names = FALSE)

# visualize
stattest_results <- read.csv("stattest_results.csv", header = T)
ggplot(stattest_results, aes(x=mean_observed_prob, y=mean_expected_prob, col=factor(call)) +
         geom_point()

# ============================================================================
# 13. Are singles from original paper really multis? need to do an open search
# ============================================================================
supp7 <- read.table("/home/lam.wayn/aas_stabilization/Supplemental_Data_7.SAAP_coordinates.tsv", sep = "\t", header = T)
supp7 <- supp7 %>% filter(name %in% relevant_genes)




