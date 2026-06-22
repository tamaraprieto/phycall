#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=GetClassifiedVariants
#SBATCH --output=ClassifiedVariantsTable.out
#SBATCH --mail-user skao@nygenome.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 5:00:00
#SBATCH --mem 1G

# Concatenate all per-cell and per-branch classified variant files into a
# single table written to ClassifiedVariantsTable.out (via SLURM --output).
#
# Output format (one line per variant):
#   chr<N>:<pos>  <cell_or_branch>  TP|FP|NA
#
# This file is then read by BuildVariantsTableCells.R to create
# variants_table_cells.rds.

MYDIR="/gpfs/commons/groups/landau_lab/skao/pta/invitro_comparison/Illumina"

# External cells (A1-H2)
for f in ${MYDIR}/results/classified/*.txt; do
    grep "^chr" "$f"
done

# Internal branches
for f in ${MYDIR}/internal/results/classified/*.txt; do
    grep "^chr" "$f"
done
