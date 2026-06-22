#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH --mem 20G
#SBATCH -t 01:00:00

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//' | sed 's/E$//')
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/../${IDLIST})
echo $SAMPLE

module purge
module load gatk/4.1.8.1
mkdir -p ${WORKDIR}/CollectAllelicCountsExomeCalls
gatk CollectAllelicCounts \
          -I ${ORIDIR}/${SAMPLE}.bqsr.marked.bam \
          -R ${RESDIR}/${REF}.fasta \
          -L  ${WORKDIR}/../SomaticCalls/Mutect${PATIENT}.filtered.bed \
          -O ${WORKDIR}/../CollectAllelicCountsExomeCalls/${SAMPLE}.WES.alleliccounts.tsv

#          -L  ${WORKDIR}/../SomaticCalls/${PATIENT}_SNP_filtered.bed \ 
