#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 05:00:00
#SBATCH --mem 50G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//' | sed 's/E$//')

##############################
# Filter the calls even more #
##############################
# Use distance to indels, mean average vaf, mean variance, missing data and more
module purge
source deactivate
module load anaconda3/10.19
source deactivate
conda deactivate
conda activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
Rscript ~/myscripts/DNAphy/src/R/ExtraFilterExome.R $PATIENT
conda deactivate

# Create a new vcf file containing only the kept mutations
module purge
module load bcftools/1.15.1
VCF="${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs-indels.hg38_multianno.vaf.vcf"
bgzip -c ${VCF} > \
        ${VCF}.gz
tabix -p vcf ${VCF}.gz

bcftools view \
        --regions-file ${WORKDIR}/../SomaticCalls/Mutect${PATIENT}.filtered.bed \
        --output-file ${WORKDIR}/../SomaticCalls/Mutect2.exome.${PATIENT}.SNVshardfiltering.vcf \
        ${VCF}.gz
