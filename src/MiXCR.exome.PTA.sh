#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 50:00:00
#SBATCH --mem 30G


source ReadConfig.sh $1

module purge
module load mixcr/3.0.13

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})
echo $SAMPLE

# READ THE EXOME BAM FILES INSTEAD OF THE FASTQ FILES!!!!! TO DOOOOO

FASTQ_1=$(grep -E  -e $SAMPLE"_" -e $SAMPLE"-" -e $SAMPLE"$" ${ORIDIR}${IDLIST} | awk -v dir=$ORIDIR '{print dir$0".R1.fastq.gz "}' | tr '\n' ' ')
FASTQ_2=$(grep -E  -e $SAMPLE"_" -e $SAMPLE"-" -e $SAMPLE"$" ${ORIDIR}${IDLIST} | awk -v dir=$ORIDIR '{print dir$0".R2.fastq.gz "}' | tr '\n' ' ')

echo "> Merging the FASTQ files"
zcat ${FASTQ_1} | gzip -c > ${WORKDIR}/${SAMPLE}_1.fastq.gz 
zcat ${FASTQ_2} | gzip -c > ${WORKDIR}/${SAMPLE}_2.fastq.gz

echo "> Running mixcr"
# Non-enriched fragmented (shotgun) Exome-Seq data. Preset exome-cdr3 is used to assemble CDR3 clonotypes, while exome-full-length additionally runs consensus contig assembly to reconstruct all available parts of V-D-J receptor rearrangement sequence
mixcr analyze exome-full-length \
    --species hsa \
    ${WORKDIR}/${SAMPLE}_1.fastq.gz ${WORKDIR}/${SAMPLE}_2.fastq.gz
      result

#mixcr analyze shotgun \
#	--species HomoSapiens \
#	--starting-material dna \
#	--report ${WORKDIR}/mixcr.${SAMPLE}.report \
#	--receptor-type xcr \
#	${WORKDIR}/${SAMPLE}_1.fastq.gz ${WORKDIR}/${SAMPLE}_2.fastq.gz \
#        ${WORKDIR}/mixcr.xCR.${SAMPLE}

rm ${WORKDIR}/${SAMPLE}_1.fastq.gz
rm ${WORKDIR}/${SAMPLE}_2.fastq.gz
