#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH --mem 2G
#SBATCH -t 05:00:00

source ReadConfig.sh $1
module purge
module load gatk/4.1.8.1

VCF=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesCellPhy.treemutkept.vcf
# Pretreatment
#SAMPLE="4295_4272T1"
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})

gatk CollectAllelicCounts \
          -I ${WORKDIR}/${SAMPLE}.recal.bam \
          -R ${RESDIR}/${REF}.fasta \
          -L ${VCF} \
          -O ${WORKDIR}/${SAMPLE}.somatic.alleliccounts.tsv
