#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem 40G

source ReadConfig.sh $1

CHR="chr21"

# BULK
SAMPLE="Invitro_Bulk"
INPUT_FILE=${WORKDIR}/${SAMPLE}.recal.bam
DEPTH=""

# SINGLE CELL
#SAMPLE="Invitro_H1" # PTA
SAMPLE="SCMDA.C6" # MDA
DEPTH=".30X"
INPUT_FILE=${WORKDIR}/${SAMPLE}.recal.${DEPTH}.bam

#module purge
#module load gatk/4.1.8.1
#gatk CollectAllelicCounts \
#        -I ${INPUT_FILE} \
#        -R ${RESDIR}/${REF}.fasta \
#	-L ${CHR} \
#        -O ${WORKDIR}/AllelicCounts.${SAMPLE}${DEPTH}.${CHR}.tsv \
#        --disable-read-filter NotDuplicateReadFilter

# -L chr20:3762925-3898969
#         -L ${RESDIR}/hg38.even.handcurated.20k.chr20_test.intervals \
#        #-L ${WORKDIR}/InvariableSites.${cell_line}.${CHR}.vcf.gz \

#grep -v -e "^@" -e ^CONTIG ${WORKDIR}/AllelicCounts.${SAMPLE}${DEPTH}.${CHR}.tsv | \
#	awk '{if ($3>0 || $4>0) {print $0}}' >  ${WORKDIR}/AllelicCountsCoverage.${SAMPLE}${DEPTH}.${CHR}.tsv

awk '{print $4/($3+$4)}' ${WORKDIR}/AllelicCountsCoverage.${SAMPLE}${DEPTH}.${CHR}.tsv |\
	 sort | uniq -c > ${WORKDIR}/AllelicCountsCoverageHistogram.${SAMPLE}${DEPTH}.${CHR}.tsv

cat  ${WORKDIR}/AllelicCountsCoverage.${SAMPLE}${DEPTH}.${CHR}.tsv | \
	awk -v chr=$CHR -v sample=$SAMPLE \
	'{sumREF+=$3; sumALT+=$4; if ($4>0){alt+=1}; if ($3>0){ref+=1}}END{print sample"\t"chr"\t"sumREF"\t"sumALT"\t"sumALT/(sumREF+sumALT)*100"\t"alt"\t"ref"\t"alt/(ref+alt)*100}' > ${WORKDIR}/AltPerc.${SAMPLE}${DEPTH}.${CHR}.txt
