#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 02:00:00
#SBATCH --mem 30G

source ReadConfig.sh $1
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})
echo $SAMPLE
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

# Create directory if it does not exist
GINKGO_INST_DIR=/gpfs/commons/groups/landau_lab/tprieto/apps/ginkgo/uploads/
mkdir -p ${GINKGO_INST_DIR}${PATIENT}

# Prepare input
module purge
module load samtools/1.9
module load bedtools/2.27.1
samtools view -bq 40 ${WORKDIR}/${SAMPLE}.dedup.bam | \
	bedtools bamtobed -i stdin | gzip > \
	${GINKGO_INST_DIR}${PATIENT}/${SAMPLE}.bed.gz
