#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 150:00:00
#SBATCH --mem 200G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

CHR=$(awk '{print $1}' ${RESDIR}/${REF}.fasta.fai | \
	grep -e "^chr[0-9|X|Y][0-9]*$" | \
	sed "${SLURM_ARRAY_TASK_ID}q;d")
GVCFs=$(awk -v dir=$WORKDIR \
	-v chr=$CHR \
	'{print "-V "dir"/HaplotypeCaller."$0"."chr".g.vcf.gz"}' \
	${ORIDIR}/${SAMPLELIST} | tr '\n' ' ')

rm -r ${PATIENT}.${CHR}
module load gatk/4.1.8.1
gatk --java-options "-Xmx4G" GenomicsDBImport \
	${GVCFs} \
	-R ${RESDIR}/${REF}.fasta \
	--genomicsdb-workspace-path ${PATIENT}.${CHR} \
	--tmp-dir ${WORKDIR} \
	-L $CHR

gatk GenotypeGVCFs \
	-V gendb://${PATIENT}.${CHR} \
	-R ${RESDIR}/${REF}.fasta \
	-O ${WORKDIR}/GenotypedGVCFs.${PATIENT}.${CHR}.vcf \
	-L $CHR
