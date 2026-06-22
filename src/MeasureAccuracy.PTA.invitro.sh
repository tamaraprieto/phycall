#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 01:00:00
#SBATCH --mem 120G

source ReadConfig.sh $1
module purge
module load gatk/4.1.8.1
module load bcftools/1.9

PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
SAMPLE=$(cat ${ORIDIR}/${SAMPLELIST} | \
        sed "${SLURM_ARRAY_TASK_ID}q;d")
echo "SAMPLE="$SAMPLE
NEWDIR=${WORKDIR}/isec_somatic-germline/FILTERED/CallingBenchmark/
mkdir -p $NEWDIR

vcf=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.treemutkept.vcf

# Select SAMPLE from VCF and just variants in chromosome X
bcftools view --samples $SAMPLE $vcf | \
bcftools filter \
        --include "CHROM=='chrX'" \
        - | bgzip -c \
        > ${NEWDIR}/TotalRaw.${SAMPLE}.vcf.gz
gatk IndexFeatureFile \
     -I ${NEWDIR}/TotalRaw.${SAMPLE}.vcf.gz

SC_COV=15
# Select variants overlapping with the sites which have more than 15 reads in bulk and no alternative
# and that show more than 15 reads in the cell
# I created ${NEWDIR}/TN4.potential.bed in previous script (contains positions with >15 depth, and no alternative reads in the bulk)
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalRaw.${SAMPLE}.vcf.gz \
	-L ${NEWDIR}/TN4.potential.bed \
        -O ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$SAMPLE\").getDP()  > ${SC_COV}'" | bash -

# Reduce TNs 
# Remove TNs which overlap with vcf positions, we don't want to work on the same site for 2 categories
# Select the real TNs for the downstream analysis
# I created ${NEWDIR}/TN4.potential.bed in previous script (contains positions with >15 depth, and no alternative reads in the bulk)
bedtools subtract \
        -a ${NEWDIR}/TN4.potential.bed \
        -b ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        > ${NEWDIR}/PotentialTNs.${SAMPLE}.bed
# Now select only those potential TNs which are not in vcf and show enough coverage in the single cell
module purge
module load gatk/4.1.8.1
gatk CollectAllelicCounts \
          -I ${WORKDIR}/${SAMPLE}.recal.bam \
          -R ${RESDIR}/${REF}.fasta \
          -L ${NEWDIR}/PotentialTNs.${SAMPLE}.bed \
          -O ${NEWDIR}/TN.${SAMPLE}.tsv
awk -v sc_cov=$SC_COV '{if ($3+$4>sc_cov){print $0}}' \
        ${NEWDIR}/TN.${SAMPLE}.tsv > \
        ${NEWDIR}/TN.${SAMPLE}.enoughcov.tsv
rm ${NEWDIR}/TN.${SAMPLE}.tsv


### DEFINE FP and TP and more TNs #####
# TUMOR VAF=0
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/FP.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$SAMPLE\").isHet()'" | bash -
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/TP.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$SAMPLE\").isHomVar()'" | bash -
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/TN.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$SAMPLE\").isHomRef()'" | bash -


grep -v "^#" ${NEWDIR}/FP.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV  '{print "FP "$1":"$2" "sample}END{print "FPt "NR" "sample" "sc_cov}' > ${NEWDIR}/Evaluation.${SAMPLE}.txt 
grep -v "^#" ${NEWDIR}/TN.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "TN "$1":"$2" "sample}END{print "TNt "NR" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt
grep -v "^#" ${NEWDIR}/TP.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "TP "$1":"$2" "sample}END{print "TPt "NR" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt
wc -l ${NEWDIR}/TN.${SAMPLE}.enoughcov.tsv | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "TNextra "$1" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt


FP=$(grep "FPt" ${NEWDIR}/Evaluation.${SAMPLE}.txt | awk '{print $2}')
TP=$(grep "TPt" ${NEWDIR}/Evaluation.${SAMPLE}.txt | awk '{print $2}')
TN=$(grep "TNt" ${NEWDIR}/Evaluation.${SAMPLE}.txt | awk '{print $2}')
TNextra=$(grep "TNextra" ${NEWDIR}/Evaluation.${SAMPLE}.txt | awk '{print $2}')

awk -v tn1=$TN \
	-v tn4=$TNextra \
	-v fp1=$FP \
	-v sample=$SAMPLE \
	-v sc_cov=$SC_COV \
	'BEGIN{fpr=(fp1)/(fp1+tn1+tn4);specif=1-fpr; print "FPR "fpr" "sample" "sc_cov"\nspecif "specif" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt

