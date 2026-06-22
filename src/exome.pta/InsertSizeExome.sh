#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 01:00:00
#SBATCH --mem 30G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//')

################
# LOAD MODULES #
################

module purge
module load picard/2.16.0

#################
# SELECT SAMPLE #
#################

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/../${IDLIST})
echo $SAMPLE

##############
# RUN PICARD #
############## 

mkdir ${WORKDIR}/InsertSizeMetrics
java -jar /nfs/sw/picard-tools/picard-tools-2.8.0/picard.jar \
      CollectInsertSizeMetrics \
      I=${ORIDIR}/${SAMPLE}.bqsr.marked.bam \
      O=${WORKDIR}InsertSizeMetrics/${SAMPLE}.insert_size_metrics.txt \
      H=${WORKDIR}InsertSizeMetrics/${SAMPLE}.insert_size_histogram.pdf
