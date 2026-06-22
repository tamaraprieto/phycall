#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 10:00:00
#SBATCH --mem 20G

source ReadConfig.sh $1
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/../${IDLIST})
echo $SAMPLE

healthyremissionlist=$(echo $IDLIST | sed 's/Samples./HealthyOrRemission./')
HEALTHY=$(head -n 1 ${ORIDIR}/../${healthyremissionlist})

SUFFIX=".bqsr.marked"
module load gatk/4.1.8.1
# ${ORIDIR}/Bams/${SAMPLE}${SUFFIX}.bam # 4295E
# gatk GetPileupSummaries \
#       -I ${ORIDIR}/${SAMPLE}${SUFFIX}.bam \
#       -V ${RESDIR}/bestpractices/small_exac_common_3.hg38.vcf.gz \
#       -L ${RESDIR}/bestpractices/small_exac_common_3.hg38.vcf.gz \
#       -O ${WORKDIR}/Getpileupsummaries.${SAMPLE}.table \
#       --disable-read-filter NotDuplicateReadFilter

#exit # first part

# Estimate contamination with CalculateContamination.
gatk CalculateContamination \
    -I ${WORKDIR}/Getpileupsummaries.${SAMPLE}.table \
    -matched ${WORKDIR}/Getpileupsummaries.${HEALTHY}.table \
    -O ${WORKDIR}/${SAMPLE}.matched.contamination.table
