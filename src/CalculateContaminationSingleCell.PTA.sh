#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 10:00:00
#SBATCH --mem 2G

source ReadConfig.sh $1
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})
echo $SAMPLE

#healthyremissionlist=$(echo $IDLIST | sed 's/Samples./HealthyOrRemission./')
#HEALTHY=$(head -n 1 ${ORIDIR}/../${healthyremissionlist})

SUFFIX=".dedup"
module load gatk/4.1.8.1
 gatk GetPileupSummaries \
       -I ${WORKDIR}/${SAMPLE}${SUFFIX}.bam \
       -V ${RESDIR}/bestpractices/small_exac_common_3.hg38.vcf.gz \
       -L ${RESDIR}/bestpractices/small_exac_common_3.hg38.vcf.gz \
       -O ${WORKDIR}/Getpileupsummaries.${SAMPLE}.table \
       --disable-read-filter NotDuplicateReadFilter

#exit # first part

# Estimate contamination with CalculateContamination
# For the matched sample I think would be possible to use the remission bulk 
gatk CalculateContamination \
    -I ${WORKDIR}/Getpileupsummaries.${SAMPLE}.table \
    -O ${WORKDIR}/${SAMPLE}.unmatched.contamination.table

#    -matched ${WORKDIR}/Getpileupsummaries.${HEALTHY}.table \
