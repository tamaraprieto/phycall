#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 01:00:00
#SBATCH --mem 120G

source ReadConfig.sh $1
module purge
module load bcftools/1.9
module load tabix/1.1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

# SELECT ONLY THE HEALTHY SAMPLE
HEALTHY=$(head -n 1 ${ORIDIR}/${CONTROL})

tag=$(head -n 1 ${RESDIR}/${REF}.fasta.fai | \
        awk '{print $1}')
HEALTHY_ID=$(zcat ${WORKDIR}/SHAPEIT/Phased.${PATIENT}.${tag}.vcf | grep "^#CHR" | tr -s '\t' '\n' | grep $HEALTHY)

VCFs=$(head -n 22 ${RESDIR}/${REF}.fasta.fai | \
	awk '{print $1}' | \
	awk -v workdir=$WORKDIR \
	-v patient=$PATIENT \
	'{print workdir"/SHAPEIT/Phased."patient"."$0".vcf.gz"}' | \
	tr -s '\n' ' ')

bcftools concat -o ${WORKDIR}/SHAPEIT/Phased.${PATIENT}.combined.vcf $VCFs

#remove phases 0|0 and 1|1
echo "remove phases 0|0 and 1|1"
bcftools view -s $HEALTHY_ID \
	-o ${WORKDIR}/SHAPEIT/Phased.${PATIENT}.healthy.vcf \
	${WORKDIR}/SHAPEIT/Phased.${PATIENT}.combined.vcf

bcftools view -e "GT=='0|0' || GT=='1|1'" \
	${WORKDIR}/SHAPEIT/Phased.${PATIENT}.healthy.vcf > \
	${WORKDIR}/SHAPEIT/Phased.${PATIENT}.healthy.hetSNPs.vcf
	
