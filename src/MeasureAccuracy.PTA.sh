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
SAMPLE=$(cat ${ORIDIR}/${SAMPLELIST} ${ORIDIR}/${CONTROL} | uniq -u | \
        grep -v "T1" | \
        sed "${SLURM_ARRAY_TASK_ID}q;d")
echo "SAMPLE="$SAMPLE
HEALTHY=$(head -n 1 ${ORIDIR}/${CONTROL})
TUMOR=$(grep "T1" ${ORIDIR}${SAMPLELIST})
NEWDIR=${WORKDIR}/isec_somatic-germline/FILTERED/CallingBenchmark/
mkdir -p $NEWDIR

#vcf=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered9.vcf
vcf=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.treemutkeptwithchrX.vcf

bcftools view --samples $HEALTHY,$TUMOR,$SAMPLE $vcf | \
bcftools filter \
        --include "CHROM=='chrX'" \
        - \
        > ${NEWDIR}/TotalRaw.${SAMPLE}.vcf

SC_COV=7
# I had already filtered the healthy bulk. It had to be 0/0 and show more than 10 reads
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalRaw.${SAMPLE}.vcf \
        -O ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$HEALTHY\").getAD().1 == 0 && vc.getGenotype(\"$TUMOR\").getDP()  > 5 && vc.getGenotype(\"$SAMPLE\").getDP()  > ${SC_COV}'" | bash -

### DEFINE TN ####
# Make sure TN4 are not present in ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf
# Select the real TNs, so they should not have a variant site
# I created ${NEWDIR}/TN4.potential.bed in previous script (contains positions with >10 in healthy, no alternative and >5 in bulk)
bedtools subtract \
        -a ${NEWDIR}/TN4.potential.bed \
        -b ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        > ${NEWDIR}/PotentialTNs.${SAMPLE}.bed
# Select sites with enough reads in the single cell which do not carry a variant site
module purge
module load gatk/4.1.8.1
gatk CollectAllelicCounts \
          -I ${WORKDIR}/${SAMPLE}.recal.bam \
          -R ${RESDIR}/${REF}.fasta \
          -L ${NEWDIR}/PotentialTNs.${SAMPLE}.bed \
          -O ${NEWDIR}/TN4.${SAMPLE}.tsv
awk -v sc_cov=$SC_COV '{if ($3+$4>sc_cov}){print $0}}' \
        ${NEWDIR}/TN4.${SAMPLE}.tsv > \
        ${NEWDIR}/TN4.${SAMPLE}.enoughcov.tsv
rm ${NEWDIR}/TN4.${SAMPLE}.tsv

### DEFINE FP and TP and more TNs #####
# TUMOR VAF=0
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/FP1.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$TUMOR\").getAD().1 == 0 && vc.getGenotype(\"$SAMPLE\").isHet()'" | bash -
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/TP1.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$TUMOR\").getAD().1 == 0 && vc.getGenotype(\"$SAMPLE\").isHomVar()'" | bash -
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/TN1.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$TUMOR\").getAD().1 == 0 && vc.getGenotype(\"$SAMPLE\").isHomRef()'" | bash -

# TUMOR 0>VAF<1
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/FP2.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$TUMOR\").getAD().0 != 0 && vc.getGenotype(\"$TUMOR\").getAD().1 != 0 && vc.getGenotype(\"$SAMPLE\").isHet()'" | bash -
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/TN2.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$TUMOR\").getAD().0 != 0 && vc.getGenotype(\"$TUMOR\").getAD().1 != 0 && vc.getGenotype(\"$SAMPLE\").isHomRef()'" | bash -
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/TP2.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$TUMOR\").getAD().0 != 0 && vc.getGenotype(\"$TUMOR\").getAD().1 != 0 && vc.getGenotype(\"$SAMPLE\").isHomVar()'" | bash -

#TUMOR VAF=1
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/FP3.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$TUMOR\").getAD().0 == 0 && vc.getGenotype(\"$SAMPLE\").isHet()'" | bash -
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/FN3.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$TUMOR\").getAD().0 == 0 && vc.getGenotype(\"$SAMPLE\").isHomRef()'" | bash -
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${NEWDIR}/TotalEnoughDepth.${SAMPLE}.vcf \
        -O ${NEWDIR}/TP3.${SAMPLE}.vcf \
        -select 'vc.getGenotype(\"$TUMOR\").getAD().0 == 0 && vc.getGenotype(\"$SAMPLE\").isHomVar()'" | bash -


grep -v "^#" ${NEWDIR}/FP1.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV  '{print "FP1 "$1":"$2" "sample}END{print "FP1t "NR" "sample" "sc_cov}' > ${NEWDIR}/Evaluation.${SAMPLE}.txt 
grep -v "^#" ${NEWDIR}/FP2.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "FP2 "$1":"$2" "sample}END{print "FP2t "NR" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt
grep -v "^#" ${NEWDIR}/FP3.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "FP3 "$1":"$2" "sample}END{print "FP3t "NR" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt
grep -v "^#" ${NEWDIR}/TN1.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "TN1 "$1":"$2" "sample}END{print "TN1t "NR" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt
grep -v "^#" ${NEWDIR}/TN2.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "TN2 "$1":"$2" "sample}END{print "TN2t "NR" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt
grep -v "^#" ${NEWDIR}/TP1.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "TP1 "$1":"$2" "sample}END{print "TP1t "NR" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt
grep -v "^#" ${NEWDIR}/TP2.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "TP2 "$1":"$2" "sample}END{print "TP2t "NR" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt
grep -v "^#" ${NEWDIR}/TP3.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "TP3 "$1":"$2" "sample}END{print "TP3t "NR" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt
grep -v "^#" ${NEWDIR}/FN3.${SAMPLE}.vcf  | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "FN3 "$1":"$2" "sample}END{print "FN3t "NR" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt
wc -l ${NEWDIR}/TN4.${SAMPLE}.enoughcov.tsv | awk -v sample=$SAMPLE -v sc_cov=$SC_COV '{print "TN4 "$1" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt


FP1=$(grep "FP1t" ${NEWDIR}/Evaluation.${SAMPLE}.txt | awk '{print $2}')
FP2=$(grep "FP2t" ${NEWDIR}/Evaluation.${SAMPLE}.txt | awk '{print $2}')
TN1=$(grep "TN1t" ${NEWDIR}/Evaluation.${SAMPLE}.txt | awk '{print $2}')
TN2=$(grep "TN2t" ${NEWDIR}/Evaluation.${SAMPLE}.txt | awk '{print $2}')
TN4=$(grep "TN4" ${NEWDIR}/Evaluation.${SAMPLE}.txt | awk '{print $2}')

awk -v tn1=$TN1 -v tn2=$TN2 \
	-v tn4=$TN4 \
	-v fp1=$FP1 -v fp2=$FP2 \
	-v sample=$SAMPLE \
	-v sc_cov=$SC_COV \
	'BEGIN{fpr=(fp1+fp2)/(fp1+fp2+tn1+tn2+tn4);specif=1-fpr; print "FPR "fpr" "sample" "sc_cov"\nspecif "specif" "sample" "sc_cov}' >> ${NEWDIR}/Evaluation.${SAMPLE}.txt

