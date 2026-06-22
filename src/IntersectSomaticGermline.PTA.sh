#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 4
#SBATCH -t 20:00:00
#SBATCH --mem 100G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

module purge
module load bcftools/1.15.1
module load vcftools/0.1.17

echo "> Compress Mutect2 files"
bgzip -c ${WORKDIR}/Mutect2.${PATIENT}.SNVs.PASS.vcf > ${WORKDIR}/Mutect2.${PATIENT}.SNVs.PASS.vcf.gz
tabix -p vcf ${WORKDIR}/Mutect2.${PATIENT}.SNVs.PASS.vcf.gz
echo "> Combine Mutect2 with HaplotypeCaller"

# For patients
bcftools isec -p ${WORKDIR}/isec_somatic-germline \
        -Oz ${WORKDIR}/Mutect2.${PATIENT}.SNVs.PASS.vcf.gz \
        ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.reheaded.vcf.gz
# For old 
#bcftools isec -p ${WORKDIR}/isec_somatic-germline \
#	-Oz ${WORKDIR}/Mutect2.${PATIENT}.SNVs.PASS.vcf.gz \
#	${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.somatic.vcf.gz


mkdir -p ${WORKDIR}/isec_somatic-germline
# For patients
VCF=${WORKDIR}/isec_somatic-germline/0003.vcf.gz
# FOR Invitro
rm -r ${WORKDIR}/isec_somatic-germline/FILTERED/; mkdir ${WORKDIR}/isec_somatic-germline/FILTERED/

echo "> Annotating the VCF with NS and VAF per sample"
# By default the values are calculated across all samples, but also per-population values can be calculated. For this, provide a file with the list of samples in the first column and a comma-separated list of populations in the second column. The file can look for example like this
# The order is important
# The FORMAT/VAF will be 0  for those samples with no reads. It would be better if it was NA 
# NS: Number of samples with data
bcftools +fill-tags $VCF \
	--threads 4 \
	-Oz -o ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered.vcf.gz \
	-- \
	--samples-file ${ORIDIR}/Samples.${PATIENT}.group.txt \
	--tags FORMAT/VAF,NS,AF
        #--samples-file ${ORIDIR}/Samples.${PATIENT}.group.txt \

# Filter the file using the new annotations
echo "> Remove sites not genotyped in half of the samples"
# half of the cells
num_half_samples=$(wc -l ${ORIDIR}/${SAMPLELIST} | awk '{print $1/2}')
bcftools filter --include "NS>=${num_half_samples}"  \
        ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered.vcf.gz | \
	awk '{if (length($5)==1){print $0}else if ($0~/^#/){print $0}}' | \
	bgzip -c \
        > ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered2.vcf.gz

echo "> Create a table with average VAF in the cells with alternative alleles"

#if [ "$PATIENT" = "molm13" ]
#then
#        cell_names=$(cat ${ORIDIR}/${SAMPLELIST} | tr -s '\n' ',' | sed 's/,$//')
#else
cell_names=$(grep "_[A-Z][0-9]" ${ORIDIR}/${SAMPLELIST} | tr -s '\n' ',' | sed 's/,$//')
#fi
BED=SC.annotations.bed
# GQ of 0 means normally that genotype has not been estimated
bcftools query \
        --samples $cell_names -f '%CHROM %POS [ %VAF,%DP,%GQ]\n' \
        ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered2.vcf.gz | \
                awk '{printf $1"\t"$2"\t";sum=0;nt=0;total_dp=0;total_GQ=0;dp=0;non_nt=0;gq=0;
                for(i=3; i<=NF; i++) {
                split($i,a,",");
                if(a[1]==0) {dp+=a[2];non_nt+=1;gq+=a[3]}        
                if(a[1]>0) {printf a[1]",";sum+=a[1];nt+=1;
                        printf a[2]",";total_dp+=a[2];
                        printf a[3]"|";total_GQ+=a[3]}}
#                if (nt>0) {printf "\t"sum/nt"\t"total_dp/nt"\t"total_GQ/nt"\t"nt"\t"dp/non_nt"\t"gq/non_nt"\n" }
                if (nt>0 && non_nt>0) {printf "\t"sum/nt"\t"total_dp/nt"\t"total_GQ/nt"\t"nt"\t"dp/non_nt"\t"gq/non_nt"\n" }
                else if (nt>0 && non_nt==0) {printf "\t"sum/nt"\t"total_dp/nt"\t"total_GQ/nt"\t"nt"\t.\t.\n" }
                else {printf "\t.\t0\t.\t.\t"nt"\t.\t.\n"}}'  | \
        sed 's/|\t/\t/g' | \
	awk 'BEGIN{print "chrom\tstart\tend\tinfo_var\tsc_vaf\tmean_dp\tmean_gq\tnum_alt\tdp_invar\tgq_invar"}{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' \
        > \
        ${WORKDIR}/isec_somatic-germline/FILTERED/$BED

# Obtain the nucleotide context
module load bedtools/2.29.0
# Obtain context
	CONTEXT_LEN=1
	bcftools view ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered2.vcf.gz | \
		grep  -v "#" | awk -v var=$CONTEXT_LEN \
        	'{OFS="\t"; print $1, $2-1-var, $2+var}' | \
        	bedtools getfasta -fi ${RESDIR}/${REF}.fasta -bed - | \
        	grep -v "^>" | awk '{print toupper($0)}' > $PATIENT".tmp"
# Obtain the alternative allele
bcftools view ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered2.vcf.gz | grep  -v "#" | \
        awk '{print $5}' | \
        paste -d"\t" $PATIENT".tmp" - | \
	awk 'BEGIN{print "ABC\tM"}{print $0}' | \
        paste -d"\t" ${WORKDIR}/isec_somatic-germline/FILTERED/$BED - | \
        awk '{for (i = 1; i <= 10; i++){printf $i"\t"};A=substr($11,1,1);B=substr($11,2,1);C=substr($11,3,1); printf A"\t"B"\t"C"\t"$12"\n"}' | \
        awk -v OFS='\t' \
        '{if ($12=="A" || $12=="G")
        {   for (i = 1; i <= 10; i++){printf $i"\t"};
	    for (i = 11; i <= 14; i++){
               if ($i=="A"){$i="T"} else if ($i=="T"){$i="A"} else if ($i=="C"){$i="G"} else if ($i=="G"){$i="C"}}
        printf $13$12">"$14$11"\n"
        }
        else{for (i = 1; i <= 10; i++){printf $i"\t"};printf $11$12">"$14$13"\n"}}' | \
	sed 's/AB>MC/ABMC/' \
	> ${WORKDIR}/isec_somatic-germline/FILTERED/${BED}2
rm $PATIENT".tmp"

mv ${WORKDIR}/isec_somatic-germline/FILTERED/${BED}2 ${WORKDIR}/isec_somatic-germline/FILTERED/$BED

# Used the bed file to annotate the vcf
echo "> Create the tsv tab with annotations"
awk 'BEGIN{printf "#"}{print $0}' ${WORKDIR}/isec_somatic-germline/FILTERED/${BED} > \
         ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab
bgzip -c ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab > \
	${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz
tabix -c 1 -f -b 2 -0 -e 3 \
	--comment '#' \
	${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz

grep "^chrom" ${WORKDIR}/isec_somatic-germline/FILTERED/${BED} | sed 's/#//' | sed 's/\t/\n/g' > newfilters.txt
#1:chrom,2:start,3:end,4:info_var,5:sc_vaf,6:mean_dp,7:mean_gq,8:num_alt,9:dp_invar,10:gq_invar,11:ABMC
a=$(head -n 5 newfilters.txt | tail -n 1)
b=$(head -n 6 newfilters.txt | tail -n 1)
c=$(head -n 7 newfilters.txt | tail -n 1)
d=$(head -n 8 newfilters.txt | tail -n 1)
e=$(head -n 9 newfilters.txt | tail -n 1)
f=$(head -n 10 newfilters.txt | tail -n 1)
g=$(head -n 11 newfilters.txt | tail -n 1)
rm newfilters.txt
echo "> Annotate the vcf files"
bcftools annotate \
       --annotations ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz \
       -h <(echo '##INFO=<ID='${a}',Number=1,Type=Float,Description="Mean VAF of cells showing at least one alternative read">\n') \
       -c CHROM,FROM,TO,-,${a},-,-,-,-,-,- \
       --output ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered3.vcf \
       ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered2.vcf.gz

bcftools annotate \
       --annotations ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz \
       -h <(echo '##INFO=<ID='${b}',Number=1,Type=Float,Description="Mean DP of cells showing at least one alternative read">\n') \
       -c CHROM,FROM,TO,-,-,${b},-,-,-,-,- \
       --output ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered4.vcf \
       ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered3.vcf

bcftools annotate \
       --annotations ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz \
       -h <(echo '##INFO=<ID='${c}',Number=1,Type=Float,Description="Mean GQ of cells showing at least one alternative read">\n') \
       -c CHROM,FROM,TO,-,-,-,${c},-,-,-,- \
       --output ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered5.vcf \
       ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered4.vcf

bcftools annotate \
       --annotations ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz \
       -h <(echo '##INFO=<ID='${d}',Number=1,Type=Integer,Description="Number of cells with alternative alleles">\n') \
       -c CHROM,FROM,TO,-,-,-,-,${d},-,-,- \
       --output ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered6.vcf \
       ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered5.vcf

bcftools annotate \
       --annotations ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz \
       -h <(echo '##INFO=<ID='${e}',Number=1,Type=Float,Description="Mean DP of cells without alternative reads">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,${e},-,- \
       --output ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered7.vcf \
       ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered6.vcf

bcftools annotate \
       --annotations ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz \
       -h <(echo '##INFO=<ID='${f}',Number=1,Type=Float,Description="Mean GQ of cells without alternative reads">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,-,${f},- \
       --output ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered8.vcf \
       ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered7.vcf

bcftools annotate \
       --annotations ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz \
       -h <(echo '##INFO=<ID='${g}',Number=1,Type=String,Description="Context of the mutation for fitting or finding de novo signatures. The string contains the 5 prime base, the reference base, the alternative allele or mutation, and the 3 prime base">\n') \
       -c CHROM,FROM,TO,-,-,-,-,-,-,-,${g} \
       --output ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered9.vcf \
       ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered8.vcf


# before May 24
#echo "> Remove sites not showing an average VAF higher than 0.35 in the cells carrying the mutation"
bcftools filter --include "${a}>0.35"  \
	${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered9.vcf \
        > ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered10.vcf


# This should be the new way but I don't want to rerun it now
#echo "> Remove sites not showing an average VAF higher than 0.35 in the cells carrying the mutation"
#bcftools filter --include "QD>=10"  \
#        ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered9.vcf \
#        > ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered10.vcf


ln -s ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered10.vcf ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.vcf


#module load bedtools
## Intersect mutect private calls with the exome capture bed file
#EXOME_PANEL_BED=/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/xgen-exome-research-panel-targets_grch38.bed
#bedtools intersect -header -a ${WORKDIR}/isec_somatic-germline/0000.vcf.gz \
#        -b $EXOME_PANEL_BED > \
#        ${WORKDIR}/isec_somatic-germline/0000-exome.vcf
#
#bedtools intersect -header -a ${WORKDIR}/isec_somatic-germline/0001.vcf.gz \
#        -b $EXOME_PANEL_BED > \
#        ${WORKDIR}/isec_somatic-germline/0001-exome.vcf
#
#bedtools intersect -header -a ${WORKDIR}/isec_somatic-germline/0002.vcf.gz \
#        -b $EXOME_PANEL_BED > \
#        ${WORKDIR}/isec_somatic-germline/0002-exome.vcf
#
#bedtools intersect -header -a ${WORKDIR}/isec_somatic-germline/0003.vcf.gz \
#        -b $EXOME_PANEL_BED > \
#        ${WORKDIR}/isec_somatic-germline/0003-exome.vcf
