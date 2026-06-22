#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 4
#SBATCH -t 10:00:00
#SBATCH --mem 60G

source ReadConfig.sh $1

module purge
module load bcftools/1.15.1

PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
mkdir -p ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping

# File containing all mutect calls (filtered and unfiltered)
VCF="${WORKDIR}/Mutect2.${PATIENT}.SNVs"

# Select variants which have not been used for phylogeny building and might actually be real
awk '{if ($7=="clustered_events" || $7=="weak_evidence" || $1~/^#/){print $0}}' ${VCF}.vcf > \
	${VCF}.filtered.vcf
bgzip -c ${VCF}.filtered.vcf > ${VCF}.filtered.vcf.gz
tabix -p vcf ${VCF}.filtered.vcf.gz

echo "> Combine Mutect2 with HaplotypeCaller"
bcftools isec -p ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/ \
        -Oz ${VCF}.filtered.vcf.gz \
        ${WORKDIR}/isec_somatic-germline/0001.vcf.gz

bcftools +fill-tags ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/0003.vcf.gz \
        --threads 4 \
        -Oz -o ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/0003.vaf.vcf \
        -- \
        --tags FORMAT/VAF

bcftools query -H -f '%CHROM:%POS[ %VAF]\n' ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/0003.vaf.vcf | \
        sed 's/:VAF//g' | sed 's/\[[0-9]*\]//g' | \
        sed 's/^# //' \
        > ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/Mutect.AD.bed


module load anaconda3/10.19
source deactivate
conda deactivate
conda activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
Rscript ${SCRIPTDIR}/R/CreateBedFileHeritableVariants.R $PATIENT
conda deactivate

bcftools view -R ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/FilteredToRescue.bed \
	${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/0003.vcf.gz > \
	${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants.vcf	


################################
# Obtain the nucleotide context
###############################
module load bedtools/2.29.0
# Obtain context
        CONTEXT_LEN=1
        bcftools view ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants.vcf | \
                grep  -v "#" | awk -v var=$CONTEXT_LEN \
                '{OFS="\t"; print $1, $2-1-var, $2+var}' | \
                bedtools getfasta -fi ${RESDIR}/${REF}.fasta -bed - | \
                grep -v "^>" | awk '{print toupper($0)}' > $PATIENT".tmp"
# Obtain the alternative allele
bcftools view ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants.vcf | grep  -v "#" | \

bcftools query \
        -f '%CHROM\t%POS\t%ALT\n' ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants.vcf | \
        #awk '{print $5}' | \
        paste -d"\t" $PATIENT".tmp" - | \
	awk '{print $2"\t"$3-1"\t"$3"\t"$1"\t"$4}' |\
        awk 'BEGIN{print "chrom\tstart\tend\tABC\tM"}{print $0}' | \
	awk '{A=substr($4,1,1);B=substr($4,2,1);C=substr($4,3,1); printf $1"\t"$2"\t"$3"\t"A"\t"B"\t"C"\t"$5"\n"}' | \
        awk -v OFS='\t' \
        '{if ($5=="A" || $5=="G")
        { for (i = 4; i <= 7; i++){
               if ($i=="A"){$i="T"} else if ($i=="T"){$i="A"} else if ($i=="C"){$i="G"} else if ($i=="G"){$i="C"}}
        printf $1"\t"$2"\t"$3"\t"$6$5">"$7$4"\n"
        }
        else{printf $1"\t"$2"\t"$3"\t"$4$5">"$7$6"\n"}}' | \
        sed 's/AB>MC/ABMC/' \
        > ${WORKDIR}/isec_somatic-germline/FILTERED/${PATIENT}.tmp2

mv ${WORKDIR}/isec_somatic-germline/FILTERED/${PATIENT}.tmp2 ${WORKDIR}/isec_somatic-germline/FILTERED/${PATIENT}.tmp

# Add chrom start end to tmp

# Used the bed file to annotate the vcf
echo "> Create the tsv tab with annotations"
awk 'BEGIN{printf "#"}{print $0}' ${WORKDIR}/isec_somatic-germline/FILTERED/${PATIENT}.tmp > \
         ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab
bgzip -c ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab > \
        ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz
tabix -c 1 -f -b 2 -0 -e 3 \
        --comment '#' \
        ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz

grep "^chrom" ${WORKDIR}/isec_somatic-germline/FILTERED/${PATIENT}.tmp | sed 's/#//' | sed 's/\t/\n/g' > newfilters.txt
#1:chrom,2:start,3:end,4:ABMC
a=$(head -n 4 newfilters.txt | tail -n 1)
echo "> Annotate the vcf files"
bcftools annotate \
       --annotations ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.tab.gz \
       -h <(echo '##INFO=<ID='${a}',Number=1,Type=String,Description="Context of the mutation for fitting or finding de novo signatures. The string contains the 5 prime base, the reference base, the alternative allele or mutation, and the 3 prime base">\n') \
       -c CHROM,FROM,TO,${a} \
       --output  ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants.mutcontext.vcf \
       ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants.vcf


################################
# Annotate the rescued variants
################################
cd /gpfs/commons/groups/landau_lab/tprieto/apps/annovar
# the output file will be named based on the prefix in out followed by "hg38_multianno.vcf"
./table_annovar.pl \
        ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants.mutcontext.vcf \
        humandb/ \
        -buildver hg38 \
        -out ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants \
        -protocol refGene,dbnsfp42c,cosmic70,avsnp150,exac03,clinvar_20220320 \
        -remove \
        -operation g,f,f,f,f,f \
        -nastring . --vcfinput


################################
# Combine all variants
################################
bgzip -c ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf > ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf.gz
tabix -p vcf ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf.gz
bgzip -c ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants.hg38_multianno.vcf > \
	${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants.hg38_multianno.vcf.gz
tabix -p vcf ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants.hg38_multianno.vcf.gz

bcftools concat --allow-overlaps \
	-o ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhyPlusRescued.vcf \
	 ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf.gz \
        ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMapping/RescuedVariants.hg38_multianno.vcf.gz 


