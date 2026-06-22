#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 02:00:00
#SBATCH --mem 10G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//'  | sed 's/E//')
echo "Patient: "$PATIENT

interval=${SLURM_ARRAY_TASK_ID}

module purge
module load gatk/4.1.8.1
gatk LearnReadOrientationModel \
        --input ${WORKDIR}/f1r2.mutect2.exome.${PATIENT}.${interval}.tar.gz \
        -O ${WORKDIR}/Read-orientation-model.exome.${PATIENT}.${interval}.tar.gz

## It seems I can use more than one contamination table (1 per sample)
CONTAMINATION_TABLES=$( cat ${ORIDIR}/${SAMPLELIST} | \
	awk -v dir=$WORKDIR \
	'{printf "-contamination-table "dir"$0".unmatched.contamination.table "}' | \
	tr '\n' ' ')

# APPLY THE FILTER
gatk \
	FilterMutectCalls \
	-R ${RESDIR}/${REF}.fasta \
	-V ${WORKDIR}/Mutect2.exome.${PATIENT}.${interval}.vcf \
	${CONTAMINATION_TABLES} \
	-stats ${WORKDIR}/Mutect2.exome.${PATIENT}.${interval}.vcf.stats \
	--ob-priors ${WORKDIR}/Read-orientation-model.exome.${PATIENT}.${interval}.tar.gz \
	-O ${WORKDIR}/Mutect2.exome.${PATIENT}.${interval}.FILTERS.vcf


# REMOVE INDELS
gatk SelectVariants \
     -R ${RESDIR}/${REF}.fasta \
     -V ${WORKDIR}/Mutect2.exome.${PATIENT}.${interval}.FILTERS.vcf \
     -O ${WORKDIR}/Mutect2.exome.${PATIENT}.${interval}.SNVs.vcf \
     --select-type-to-include SNP     

module load bcftools/1.9
bgzip -c ${WORKDIR}/Mutect2.exome.${PATIENT}.${interval}.FILTERS.vcf > ${WORKDIR}/Mutect2.exome.${PATIENT}.${interval}.FILTERS.vcf.gz
tabix -p vcf ${WORKDIR}/Mutect2.exome.${PATIENT}.${interval}.FILTERS.vcf.gz
bgzip -c ${WORKDIR}/Mutect2.exome.${PATIENT}.${interval}.SNVs.vcf > ${WORKDIR}/Mutect2.exome.${PATIENT}.${interval}.SNVs.vcf.gz
tabix -p vcf ${WORKDIR}/Mutect2.exome.${PATIENT}.${interval}.SNVs.vcf.gz
