#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem 5G

source ReadConfig.sh $1
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})
CHR=$2

module load gatk/4.1.8.1
gatk HaplotypeCaller \
	-R ${RESDIR}/${REF}.fasta \
	-I ${WORKDIR}/${SAMPLE}.recal.bam \
	-O ${WORKDIR}/HaplotypeCaller.${SAMPLE}.${CHR}.g.vcf.gz \
	-ERC GVCF \
	-L ${CHR}
