#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 05:00:00
#SBATCH --mem 100G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

module load bcftools/1.9

# old list of SNPs passing filters
#VCF2=/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/4295E_Results_TP/0829-Tranche99.9_4295_Exome_SNP_ScHet.recode.vcf
#VCF2=/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/4295E_Results_TP/4295_somatic.unidentified.recode.vcf
# new list of SNPs passing filters(Yakun 10 Feb 2022)

VCF2="/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/4295E_Results_TP/4295_SNP_Indel_filtered.34.vcf"
bgzip -c ${VCF2} > ${VCF2}.gz
tabix -p vcf ${VCF2}.gz

VCF1=${WORKDIR}/Mutect2.${PATIENT}.SNVs.PASS.vcf.gz
VCF1=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.vcf

bcftools isec -p ${WORKDIR}/isec_WGS-WES \
	-Oz ${VCF1} \
	${VCF2}.gz

module load bedtools
# Intersect mutect private calls with the exome capture bed file
EXOME_PANEL_BED=/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/xgen-exome-research-panel-targets_grch38.bed
bedtools intersect -header -a ${WORKDIR}/isec_WGS-WES/0000.vcf.gz \
	-b $EXOME_PANEL_BED > \
	${WORKDIR}/isec_WGS-WES/0000-exome.vcf

bedtools intersect -header -a ${WORKDIR}/isec_WGS-WES/0001.vcf.gz \
        -b $EXOME_PANEL_BED > \
        ${WORKDIR}/isec_WGS-WES/0001-exome.vcf

bedtools intersect -header -a ${WORKDIR}/isec_WGS-WES/0002.vcf.gz \
        -b $EXOME_PANEL_BED > \
        ${WORKDIR}/isec_WGS-WES/0002-exome.vcf

bedtools intersect -header -a ${WORKDIR}/isec_WGS-WES/0003.vcf.gz \
        -b $EXOME_PANEL_BED > \
        ${WORKDIR}/isec_WGS-WES/0003-exome.vcf
