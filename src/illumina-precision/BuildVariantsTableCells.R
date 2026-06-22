#!/usr/bin/env Rscript
# BuildVariantsTableCells.R
#
# Reads the per-cell / per-branch classified variant files produced by
# ClassifyVariants.py and saves variants_table_cells.rds — the precision
# input for Figures.BALL.PTA.Rmd (Figure 4D / Figure S6C).
#
# Run after all CountandClassify array jobs have finished:
#   Rscript src/illumina-precision/BuildVariantsTableCells.R
#
# ClassifyVariants.py output columns (tab-separated, one line per variant):
#   1  chr:variant
#   2  hetSNP
#   3  branch
#   4  category        (TP / FP / NA)
#   5  REF_ALT
#   6  REF_REF
#   7  ALT_REF
#   8  ALT_ALT
#   9  matching_reads
#   10 nonmatching_reads
#   11 total_reads
#
# Output columns in the RDS:
#   branch   — variant-branch label from the somatic BED (e.g. "Invitro_A1")
#   category — "TP", "FP", or "NA"
#   reads    — total read depth at the variant position (col 11)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
})

MYDIR  <- "/gpfs/commons/groups/landau_lab/skao/pta/invitro_comparison/Illumina"
OUTFILE <- file.path(MYDIR, "variants_table_cells.rds")

# ---------------------------------------------------------------------------
# Helper: parse a .counts.classified.txt file produced by ClassifyVariants.py
#
# Keeps only data lines (starting with "chr") and skips the header line,
# path echo, and summary lines at the end.
# ---------------------------------------------------------------------------
parse_classified <- function(path) {
  lines     <- readLines(path, warn = FALSE)
  chr_lines <- grep("^chr", lines, value = TRUE)
  if (length(chr_lines) == 0L) return(NULL)

  mat <- str_split_fixed(chr_lines, "\t", 11)
  tibble(
    locus       = mat[, 1],
    branch      = mat[, 3],
    category    = mat[, 4],
    reads       = as.integer(mat[, 11])
  )
}

# ---------------------------------------------------------------------------
# 1. External cells  (A1 – H2)
# ---------------------------------------------------------------------------
external_cells <- c("A1","A2","B1","B2","C1","C2","D1","D2",
                    "E1","E2","F1","F2","G1","G2","H1","H2")

ext_classified <- map_dfr(external_cells, function(cell) {
  f <- file.path(MYDIR, "results", "classified",
                 paste0(cell, "_Illumina.counts.classified.txt"))
  if (!file.exists(f)) { message("Missing: ", f); return(NULL) }
  parse_classified(f)
})

# ---------------------------------------------------------------------------
# 2. Internal branches
# ---------------------------------------------------------------------------
branch_map_file <- file.path(MYDIR, "internal", "Internal_branches_bams.txt")
branch_map      <- read_tsv(branch_map_file, col_names = c("branch", "bam"),
                             show_col_types = FALSE)

int_classified <- map_dfr(branch_map$branch, function(b) {
  f <- file.path(MYDIR, "internal", "results", "classified",
                 paste0(b, ".counts.classified.txt"))
  if (!file.exists(f)) { message("Missing: ", f); return(NULL) }
  parse_classified(f)
})

# ---------------------------------------------------------------------------
# 3. Combine and save RDS
# ---------------------------------------------------------------------------
variants_table_cells <- bind_rows(ext_classified, int_classified) %>%
  select(branch, category, reads)

message("Rows: ",   nrow(variants_table_cells))
message("Cells/branches: ", n_distinct(variants_table_cells$branch))
message("TP: ", sum(variants_table_cells$category == "TP", na.rm = TRUE))
message("FP: ", sum(variants_table_cells$category == "FP", na.rm = TRUE))
message("NA: ", sum(variants_table_cells$category == "NA", na.rm = TRUE))
message("Missing reads (NA): ", sum(is.na(variants_table_cells$reads)))

saveRDS(variants_table_cells, file = OUTFILE)
message("Saved: ", OUTFILE)
