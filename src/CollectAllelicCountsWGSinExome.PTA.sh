#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH --mem 80G
#SBATCH -t 05:00:00

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})
echo $SAMPLE

module purge
module load gatk/4.1.8.1
mkdir -p ${WORKDIR}/CollectAllelicCountsExomeCalls
gatk CollectAllelicCounts \
          -I ${WORKDIR}/${SAMPLE}.recal.bam \
          -R ${RESDIR}/${REF}.fasta \
          -L  ${WORKDIR}/../../../2021-07-19_BioSkryb/SomaticCalls/${PATIENT}_SNP_filtered.bed \
          -O ${WORKDIR}/CollectAllelicCountsExomeCalls/${SAMPLE}.alleliccounts.tsv
