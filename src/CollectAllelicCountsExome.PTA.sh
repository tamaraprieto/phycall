#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH --mem 80G
#SBATCH -t 05:00:00


source ReadConfig.sh $1
PATIENT=$(echo $1 | sed 's/.txt//' | sed 's/Config.//')

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})
echo $SAMPLE

module purge
module load gatk/4.1.8.1
EXOME_BAM_DIR=/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/
#RESDIR=
#REF=

gatk CollectAllelicCounts \
          -I ${EXOME_BAM_DIR}/4295E_Results/Bams/${SAMPLE}.bqsr.marked.bam \
          -R ${RESDIR}/${REF}.fasta \
          -L ${EXOME_BAM_DIR}/SomaticCalls/${PATIENT}_SNP_filtered.bed \
          -O ${WORKDIR}/CollectAllelicCounts/${PATIENT}_${SAMPLE}.exome.alleliccounts.tsv
