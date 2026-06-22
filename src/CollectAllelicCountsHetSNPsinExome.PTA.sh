#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH --mem 80G
#SBATCH -t 05:00:00

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//' | sed 's/E$//')
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/../${IDLIST})
echo $SAMPLE

module purge
module load gatk/4.1.8.1
mkdir -p ${WORKDIR}/../CollectAllelicCountsHetGermlineSNPs
gatk CollectAllelicCounts \
          -I ${ORIDIR}/${SAMPLE}.bqsr.marked.bam \
          -R ${RESDIR}/${REF}.fasta \
          -L  /gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/${PATIENT}/RESULTS//HaplotypeCaller.${PATIENT}.vqsr.heterozygous.vcf.gz \
          -O ${WORKDIR}/../CollectAllelicCountsHetGermlineSNPs/${SAMPLE}.WES.alleliccounts.tsv


