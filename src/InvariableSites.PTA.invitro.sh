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
HEALTHY="Invitro_Bulk"

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
          -O ${NEWDIR}/TN4.tsv
awk '{if ($3+$4>15 && $4==0){print $0}}' \
        ${NEWDIR}/TN4.tsv | awk '{print $1"\t"$2-1"\t"$2}' > \
        ${NEWDIR}/TN4.potential.bed
rm ${NEWDIR}/TN4.tsv


