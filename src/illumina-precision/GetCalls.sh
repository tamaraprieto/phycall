#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=GetCalls
#SBATCH --output=GetCalls.out
#SBATCH --mail-user skao@nygenome.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem 25G

# For each external cell (A1-H2):
#   1. Expand the cell's somatic-variant BED by 100 bp
#   2. Intersect with the genome-wide het-SNP VCF to get a phasing-context BED
#
# Output: bed/expanded_hetSNPs/${cell}.expanded.hetSNPs.bed

module load bedtools

MYDIR="/gpfs/commons/groups/landau_lab/skao/pta/invitro_comparison/Illumina"
HETSNP_VCF="/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/Invitro/RESULTS/HaplotypeCaller.Invitro.vqsr.heterozygous.vcf.gz"
CHROMSIZES="/gpfs/commons/groups/landau_lab/skao/pta/hg38.chrom.sizes"

for cell in A1 A2 B1 B2 C1 C2 D1 D2 E1 E2 F1 F2 G1 G2 H1 H2
do
    # expand bed intervals by 100 bp
    bedtools slop \
        -i ${MYDIR}/bed/${cell}.bed \
        -g ${CHROMSIZES} \
        -header \
        -b 100 \
        > ${MYDIR}/bed/expanded/${cell}.expanded.bed

    # intersect expanded intervals with het-SNP VCF
    bedtools intersect \
        -wa \
        -wb \
        -a ${MYDIR}/bed/expanded/${cell}.expanded.bed \
        -b ${HETSNP_VCF} \
        -header \
        > ${MYDIR}/bed/expanded_hetSNPs/${cell}.expanded.hetSNPs.bed

done
