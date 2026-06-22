#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 05:00:00
#SBATCH --mem 10G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

interval=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${RESDIR}/hg38.even.handcurated.20k.intervals)
CHR=$(echo $interval | sed 's/:.*//')
GVCFs=$(awk -v dir=$WORKDIR \
	-v chr=$CHR \
	'{print "-V "dir"/HaplotypeCaller."$0"."chr".g.vcf.gz"}' \
	${ORIDIR}/${SAMPLELIST} | tr '\n' ' ')

rm -r ${PATIENT}.${SLURM_ARRAY_TASK_ID}
module load gatk/4.1.8.1
gatk --java-options "-Xmx4G" GenomicsDBImport \
	${GVCFs} \
	-R ${RESDIR}/${REF}.fasta \
	--genomicsdb-workspace-path ${PATIENT}.${SLURM_ARRAY_TASK_ID} \
	--tmp-dir ${WORKDIR} \
	-L $interval

gatk GenotypeGVCFs \
    -V gendb://${PATIENT}.${SLURM_ARRAY_TASK_ID} \
    -R ${RESDIR}/${REF}.fasta \
    -O ${WORKDIR}/GenotypedGVCFs.${PATIENT}.${SLURM_ARRAY_TASK_ID}.vcf \
    -L $interval

