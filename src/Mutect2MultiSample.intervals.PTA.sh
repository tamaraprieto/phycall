#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 150:00:00
#SBATCH --mem-per-cpu 50G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
echo "Patient: "$PATIENT

interval=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${RESDIR}/hg38.even.handcurated.20k.intervals)
CHR=$(echo $interval | sed 's/:.*//')
GVCFs=$(awk -v dir=$WORKDIR \
        -v chr=$CHR \
        '{print "-V "dir"/HaplotypeCaller."$0"."chr".g.vcf.gz"}' \
        ${ORIDIR}/${SAMPLELIST} | tr '\n' ' ')

# SELECT NORMAL SAMPLE
NORMAL_NAME=$(head -n 1 ${ORIDIR}/${CONTROL})
#NORMAL_NAME="4295_UndeterminedBarcode_TTAGGATAGA_CACCTTAATC_CACCTTAATC_TTAGGATAGA"

module purge
module load bedtools/2.29.0
module load gatk/4.1.8.1

echo "> Run Mutect for interval "${interval}", which corresponds to number "${SLURM_ARRAY_TASK_ID}
SAMPLES=$(awk -v dir=$WORKDIR '{print "-I "dir"/"$0".recal.bam"}' ${ORIDIR}/${SAMPLELIST} | tr '\n' ' ')
echo $SAMPLES
gatk --java-options "-Xmx20G" Mutect2 \
	-R ${RESDIR}/${REF}.fasta \
	${SAMPLES} \
        -normal ${NORMAL_NAME} \
	--germline-resource ${RESDIR}/bestpractices/af-only-gnomad.hg38.vcf.gz \
	-L $interval \
	-O ${WORKDIR}/Mutect2.${PATIENT}.${SLURM_ARRAY_TASK_ID}.vcf \
	--f1r2-tar-gz ${WORKDIR}/f1r2.mutect2.${PATIENT}.${SLURM_ARRAY_TASK_ID}.tar.gz
