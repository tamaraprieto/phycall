#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=GetCallsInternal
#SBATCH --output=GetCallsInternal.out
#SBATCH --mail-user skao@nygenome.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem 25G

# Same as ../GetCalls.sh but for internal (multi-cell) branches.
# Branch names are groups of cells that share a common ancestor,
# e.g. "H2H1D2D1" = the internal branch subtending cells H2, H1, D2, D1.
#
# Output: internal/bed/expanded_hetSNPs/${branch}.expanded.hetSNPs.bed

module load bedtools

MYDIR="/gpfs/commons/groups/landau_lab/skao/pta/invitro_comparison/Illumina/internal"
HETSNP_VCF="/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/Invitro/RESULTS/HaplotypeCaller.Invitro.vqsr.heterozygous.vcf.gz"
CHROMSIZES="/gpfs/commons/groups/landau_lab/skao/pta/hg38.chrom.sizes"

for branch in B2A2C2B1G1G2F1F2 B2C2B1G1G2F1F2 B2G1G2F1F2 B2G2F1F2 G2F1F2 G2F1 C2B1 \
              E1H2H1D2D1E2A1 E1H2H1D2D1E2 H2H1D2D1E2 H2H1D2D1 H2H1 D2D1
do
    bedtools slop \
        -i ${MYDIR}/bed/${branch}.bed \
        -g ${CHROMSIZES} \
        -header \
        -b 100 \
        > ${MYDIR}/bed/expanded/${branch}.expanded.bed

    bedtools intersect \
        -wa \
        -wb \
        -a ${MYDIR}/bed/expanded/${branch}.expanded.bed \
        -b ${HETSNP_VCF} \
        -header \
        > ${MYDIR}/bed/expanded_hetSNPs/${branch}.expanded.hetSNPs.bed
done
