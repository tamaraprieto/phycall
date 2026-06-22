#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=CountandClassifyIllumina
#SBATCH --output=CountandClassifyIllumina.out
#SBATCH --mail-user skao@nygenome.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 150:00:00
#SBATCH --mem 20G

# Per-cell read counting and TP/FP classification for external cells (A1-H2).
# Run as an array job:
#   sbatch --array=1-16 CountandClassifyIlluminaReads.sh
#
# Requires:
#   CountReads.py    — counts haplotype-phased reads at each variant position
#   ClassifyVariants.py — classifies each site as TP / FP / NA
# Both scripts are at: /gpfs/commons/groups/landau_lab/skao/tools/MyScripts/
# Copy them to this directory before running.
#
# Input (per cell):
#   Invitro_${SAMPLE}.recal.bam            — Illumina WGS BAM
#   bed/expanded_hetSNPs/${SAMPLE}.expanded.hetSNPs.bed  — from GetCalls.sh
#
# Output (per cell):
#   results/counts/${SAMPLE}_Illumina.counts.txt
#   results/classified/${SAMPLE}_Illumina.counts.classified.txt

MYDIR="/gpfs/commons/groups/landau_lab/skao/pta/invitro_comparison/Illumina"
BAMDIR="/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/Invitro/RESULTS"
SCRIPTS="$(dirname $0)"   # expects CountReads.py and ClassifyVariants.py here

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /gpfs/commons/groups/landau_lab/skao/pta/invitro_comparison/cells.txt)

echo "Processing cell: ${SAMPLE}"

module load python

# Step 1: count phased reads at each variant site
python ${SCRIPTS}/CountReads.py \
    --bamfile ${BAMDIR}/Invitro_${SAMPLE}.recal.bam \
    --bedfile ${MYDIR}/bed/expanded_hetSNPs/${SAMPLE}.expanded.hetSNPs.bed \
> ${MYDIR}/results/counts/${SAMPLE}_Illumina.counts.txt

module unload python
module load python/3.11.0

# Step 2: classify each site as TP / FP / NA
python ${SCRIPTS}/ClassifyVariants.py \
    --countsfile ${MYDIR}/results/counts/${SAMPLE}_Illumina.counts.txt \
    --bedfile    ${MYDIR}/bed/expanded_hetSNPs/${SAMPLE}.expanded.hetSNPs.bed \
> ${MYDIR}/results/classified/${SAMPLE}_Illumina.counts.classified.txt
