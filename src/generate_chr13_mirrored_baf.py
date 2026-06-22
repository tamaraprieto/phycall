#!/usr/bin/env python
"""
Regenerate region_1_chr13_5k_imbalance.csv for Figure 4F (Illumina/Invitro_hg38).

Reads dict_chr13.pkl from the BAF_5000 folder and extracts mirrored BAF
(|BAF - 0.5|) for each cell across chr13:40,000,000-43,000,000 in 5kb windows.

Output CSV: 16 rows (cells, no Bulk) x 600 columns (5kb windows).
After t(fread(...)) in R, becomes 600 rows x 16 columns.
Figures.BALL.PTA.Rmd assigns rownames seq(from=40000000, by=5000, ...) and
plots only the sub-region 41,400,000–42,900,000 (containing the chr13q14.11 deletion
at 41,875,000–42,400,000).

Usage:
    python generate_chr13_mirrored_baf.py [--out /path/to/output.csv]
"""

import pickle
import numpy as np
import pandas as pd
import argparse
import os

# ── Paths ──────────────────────────────────────────────────────────────────────
BAF_DIR  = "/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/Invitro/RESULTS/BAF/BAF_5000/"
PKL_FILE = os.path.join(BAF_DIR, "dict_chr13.pkl")
OUT_CSV  = "/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/Invitro/RESULTS/BAF/region_1_chr13_5k_imbalance.csv"

# ── Parameters ────────────────────────────────────────────────────────────────
# CHR_START / CHR_END define the genomic window extracted into the CSV.
# The Rmd assigns rownames as seq(from=CHR_START, by=BIN_SIZE, ...) and then
# filters the plot to interval > 41,400,000 and < 42,900,000, so the CSV must
# start at CHR_START and extend past 42,900,000.
# Deletion of interest: chr13:41,875,000–42,400,000 (chr13q14.11).
PATIENT      = "Invitro"
BIN_SIZE     = 5000
CHR_START    = 40_000_000   # must match seq(from=...) in Figures.BALL.PTA.Rmd
CHR_END      = 44_000_000   # covers the plot filter (<42,900,000) with margin
FILTER_SNPS  = 2    # minimum heterozygous SNPs per bin
FILTER_COUNTS= 10   # minimum total read counts per bin

# Cell order: no Bulk (matches colnames in Figures.BALL.PTA.Rmd)
CELLS = ["A1", "A2", "B1", "B2", "C1", "C2",
         "D1", "D2", "E1", "E2", "F1", "F2",
         "G1", "G2", "H1", "H2"]


def main(out_path):
    print(f"Loading {PKL_FILE} ...")
    with open(PKL_FILE, "rb") as f:
        d = pickle.load(f, encoding="latin1")

    # ── Determine bin indices for the region ──────────────────────────────────
    # Bins are 0-indexed from the start of the chromosome
    # A bin covers [i*BIN_SIZE, (i+1)*BIN_SIZE)
    idx_start = CHR_START // BIN_SIZE
    idx_end   = CHR_END   // BIN_SIZE
    n_bins    = idx_end - idx_start
    print(f"Extracting bins {idx_start}–{idx_end} ({n_bins} bins, "
          f"chr13:{CHR_START:,}–{CHR_END:,})")

    mirrored = {}

    for cell in CELLS:
        key = f"{PATIENT}_{cell}"
        if key not in d:
            print(f"  WARNING: key '{key}' not found in pkl — filling with NaN")
            mirrored[cell] = [np.nan] * n_bins
            continue

        entry = d[key]
        baf_raw    = np.array(entry["BAF"]).flatten()
        snps_raw   = np.array(entry["SNP"]).flatten()
        counts_raw = np.array(entry["COUNTS"]).flatten()

        baf_region    = baf_raw[idx_start:idx_end]
        snps_region   = snps_raw[idx_start:idx_end]
        counts_region = counts_raw[idx_start:idx_end]

        # Apply filters: bins with too few SNPs or reads become NaN
        baf_filtered = np.where(
            (snps_region >= FILTER_SNPS) & (counts_region >= FILTER_COUNTS),
            baf_raw[idx_start:idx_end],
            np.nan
        )

        # Mirror: |BAF - 0.5|
        mirrored[cell] = np.abs(baf_filtered - 0.5)

    # ── Build DataFrame: rows = bins, columns = cells ─────────────────────────
    df = pd.DataFrame(mirrored)  # shape: n_bins x 16
    # Transpose so cells are rows (matching Jean's original CSV format)
    df_out = df.T  # shape: 16 x n_bins

    print(f"Output shape: {df_out.shape}  (rows=cells, cols=windows)")
    print(f"Writing to {out_path} ...")
    df_out.to_csv(out_path, index=False, header=False)
    print("Done.")
    print()
    print("In Figures.BALL.PTA.Rmd, update the Illumina path to:")
    print(f"  {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=OUT_CSV,
                        help="Output CSV path (default: %(default)s)")
    args = parser.parse_args()
    main(args.out)
