#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=CountandClassifyInternal
#SBATCH --output=CountandClassifyInternal.out
#SBATCH --mail-user skao@nygenome.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 150:00:00
#SBATCH --mem 25G

# Same logic as ../CountandClassifyIlluminaReads.sh but for internal branches.
# Each internal branch is a merged BAM covering the cells that subtend it.
# The sample-to-BAM mapping is in Internal_branches_bams.txt (tab-separated:
# branch_name \t bam_suffix).
#
# Run as array job:
#   sbatch --array=1-N CountandClassifyReads.sh
# where N = number of lines in Internal_branches_bams.txt
#
# Requires CountReads.py and ClassifyVariants.py (see ../CountandClassifyIlluminaReads.sh).

MYDIR="/gpfs/commons/groups/landau_lab/skao/pta/invitro_comparison/Illumina/internal"
BAMDIR="/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/Invitro/RESULTS"
SCRIPTS="$(dirname $0)/.."   # CountReads.py and ClassifyVariants.py live one level up

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${MYDIR}/Internal_branches_bams.txt | cut -f1)
BAM=$(sed    "${SLURM_ARRAY_TASK_ID}q;d" ${MYDIR}/Internal_branches_bams.txt | cut -f2)

echo "Branch: ${SAMPLE}  BAM: ${BAM}"

module load python

python ${SCRIPTS}/CountReads.py \
    --bamfile ${BAMDIR}/Invitro_${BAM}.recal.bam \
    --bedfile ${MYDIR}/bed/expanded_hetSNPs/${SAMPLE}.expanded.hetSNPs.bed \
> ${MYDIR}/results/counts/${SAMPLE}.counts.txt

module unload python
module load python/3.11.0

python ${SCRIPTS}/ClassifyVariants.py \
    --countsfile ${MYDIR}/results/counts/${SAMPLE}.counts.txt \
    --bedfile    ${MYDIR}/bed/expanded_hetSNPs/${SAMPLE}.expanded.hetSNPs.bed \
> ${MYDIR}/results/classified/${SAMPLE}.counts.classified.txt
