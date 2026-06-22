#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem 50G

source ReadConfig.sh $1

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})

samples=$(grep -E -e $SAMPLE"-"  -e $SAMPLE"_" -e $SAMPLE"$" -e $SAMPLE"\." ${ORIDIR}${IDLIST} | awk -v dir=$WORKDIR '{print "I="dir$0".sorted.bam"}' | tr '\n' ' ')

module purge
module load picard/2.18.17

java -jar /nfs/sw/picard-tools/picard-tools-2.18.17/picard.jar \
	MarkDuplicates \
	${samples} \
	OUTPUT=${WORKDIR}/${SAMPLE}.dedup.bam \
	CREATE_INDEX=true \
	REMOVE_DUPLICATES=false \
	TMP_DIR=${WORKDIR} \
	M=${WORKDIR}/Duplicates_${SAMPLE}.txt \
	VALIDATION_STRINGENCY=LENIENT \
	TAGGING_POLICY=All

