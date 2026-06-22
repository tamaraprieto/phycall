#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 01:00:00
#SBATCH --mem 120G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
SAMPLE=$(cat ${ORIDIR}/${SAMPLELIST} ${ORIDIR}/${CONTROL} | uniq -u | \
	grep -v "T1" | \
	sed "${SLURM_ARRAY_TASK_ID}q;d")
echo "SAMPLE="$SAMPLE

HEALTHY=$(head -n 1 ${ORIDIR}/${CONTROL})
TUMOR=$(grep "T1" ${ORIDIR}${SAMPLELIST})

module purge
module load gatk/4.1.8.1
module load bcftools
NEWDIR=${WORKDIR}/isec_somatic-germline/FILTERED/CallingBenchmark/
mkdir -p $NEWDIR

# Remove sites with not enough coverage in healthy bulk
gatk CollectAllelicCounts \
          -I ${WORKDIR}/${HEALTHY}.recal.bam \
          -R ${RESDIR}/${REF}.fasta \
          -L chrX \
          -O ${NEWDIR}/TN4.${HEALTHY}.tsv
awk '{if ($3+$4>10 && $4==0){print $0}}' \
        ${NEWDIR}/TN4.${HEALTHY}.tsv | awk '{print $1"\t"$2-1"\t"$2}' > \
        ${NEWDIR}/TN4.${HEALTHY}.morethan10.bed
rm ${NEWDIR}/TN4.${HEALTHY}.tsv

# Remove sites with not enough coverage in tumor bulk
gatk CollectAllelicCounts \
          -I ${WORKDIR}/${TUMOR}.recal.bam \
          -R ${RESDIR}/${REF}.fasta \
          -L  ${NEWDIR}/TN4.${HEALTHY}.morethan10.bed \
          -O ${NEWDIR}/TN4.${TUMOR}.tsv
awk '{if ($3+$4>5){print $0}}' \
        ${NEWDIR}/TN4.${TUMOR}.tsv | awk '{print $1"\t"$2-1"\t"$2}' > \
        ${NEWDIR}/TN4.potential.bed
rm ${NEWDIR}/TN4.${TUMOR}.tsv

