#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 10:00:00
#SBATCH --mem 40G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
module purge
module load bedtools/2.29.0
module load bcftools/1.9
module load vcftools/0.1.17

VCF=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesCellPhy.treemutkept.vcf

bcftools query \
        -f '%CHROM:%POS\t%ABMC\n' $VCF > \
	${WORKDIR}/MutationalContext.txt	

bcftools query \
        -f '%CHROM:%POS\t%num_alt\t%NS\t%ABMC\n' $VCF > \
        ${WORKDIR}/MutationalContextDirectGenotypesVCF.nomapping.txt
	

# For obtaining the germline signatures
VCF=${WORKDIR}/recall/VariantSitesRecall.vcf
bcftools query \
        -f '%CHROM:%POS\t%ABMC\n' $VCF > \
        ${WORKDIR}/MutationalContextGermline.txt
