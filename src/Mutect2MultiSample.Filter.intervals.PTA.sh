#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 02:00:00
#SBATCH --mem 30G


source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
echo "Patient: "$PATIENT

number_intervals=$(wc -l ${RESDIR}/hg38.even.handcurated.20k.intervals | awk '{print $1}')
interval=$(seq 1 $number_intervals | sed "${SLURM_ARRAY_TASK_ID}q;d")

module purge
module load gatk/4.1.8.1
gatk LearnReadOrientationModel \
        --input ${WORKDIR}/f1r2.mutect2.${PATIENT}.${interval}.tar.gz \
        -O ${WORKDIR}/Read-orientation-model.${PATIENT}.${interval}.tar.gz


## It seems I can use more than one contamination table (1 per sample)
CONTAMINATION_TABLES=$( cat ${ORIDIR}/${SAMPLELIST} | \
	awk -v dir=$WORKDIR \
	'{printf "-contamination-table "dir"$0".unmatched.contamination.table "}' | \
	tr '\n' ' ')

# APPLY THE FILTER
gatk \
	FilterMutectCalls \
	-R ${RESDIR}/${REF}.fasta \
	-V ${WORKDIR}/Mutect2.${PATIENT}.${interval}.vcf \
	${CONTAMINATION_TABLES} \
	-stats ${WORKDIR}/Mutect2.${PATIENT}.${interval}.vcf.stats \
	--ob-priors ${WORKDIR}/Read-orientation-model.${PATIENT}.${interval}.tar.gz \
	-O ${WORKDIR}/Mutect2.${PATIENT}.${interval}.FILTERS.vcf


# REMOVE INDELS
gatk SelectVariants \
     -R ${RESDIR}/${REF}.fasta \
     -V ${WORKDIR}/Mutect2.${PATIENT}.${interval}.FILTERS.vcf \
     -O ${WORKDIR}/Mutect2.${PATIENT}.${interval}.SNVs.vcf \
     --select-type-to-include SNP     

