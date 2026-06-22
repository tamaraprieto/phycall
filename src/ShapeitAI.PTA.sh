#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type END
#SBATCH --cpus-per-task 1
#SBATCH -t 10:00:00
#SBATCH --mem 40G

source ReadConfig.sh $1

# Create a vcf with a single chromosome and just the healthy sample
CHR=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${RESDIR}/${REF}.fasta.fai | awk '{print $1}')

# Split VCF files by chromosome
module purge
#module load bcftools/1.9
#echo $CHR
#bcftools view ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.vcf.gz \
#        --regions ${CHR} \ |
#	--samples ${HEALTHY} | bgzip -c \
#        > ${WORKDIR}/HaplotypeCaller.${PATIENT}.${HEALTHY}.snps.${CHR}.vcf.gz



# Phase the SNPs
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
#PATIENT+=".hg38"
module load shapeit/4.2.2

SHAPEITDIR="/gpfs/commons/groups/landau_lab/tprieto/jlaval/SHAPEIT/"


#if using hg38
refname=$(basename $RESDIR)
REFERENCE="${SHAPEITDIR}${refname}/1kGP_high_coverage_Illumina.${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
#map b38
#region chr${CHR}
#if using hg19
#REFERENCE=/gpfs/commons/home/jlaval/SHAPEIT/${refname}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#map b37
#region ${CHR}

mkdir -p ${WORKDIR}/SHAPEIT
shapeit4 --input ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.vcf.gz \
         --map ${SHAPEITDIR}${refname}/${CHR}.b38.gmap.gz \
         --reference ${REFERENCE} \
	 --region ${CHR} \
	 --output ${WORKDIR}/SHAPEIT/Phased.${PATIENT}.${CHR}.vcf.gz \
	 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m \
	 --pbwt-depth 8

#--mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m
#--mcmc-iterations 5b,1p,1b,1p,1b,1p,5m
#default: pbwt 4
