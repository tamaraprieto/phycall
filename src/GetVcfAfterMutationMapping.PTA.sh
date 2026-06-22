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

##########################################
# For obtaining VCF with mapped variants #
##########################################

PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
# Run it only once
## We are going to run the mutation calling on the variants that we have in the tree (after forcing the mutation mapping)

if [ $PATIENT = "Invitro" ]
then
        VCF="${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhyPlusRescued.vcf"
else
        VCF="${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf"
fi
#VCF=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf
bgzip -c ${VCF} > \
        ${VCF}.gz
tabix -p vcf ${VCF}.gz

bcftools view \
        --regions-file ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesCellPhyAfterMutationMapping.genome.bed \
        --output-file ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesCellPhy.treemutkept.vcf \
        ${VCF}.gz

###############################
# For obtaining chrX on males #
###############################

VCF="${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.vcf"
bgzip -c ${VCF} > \
        ${VCF}.gz
tabix -p vcf ${VCF}.gz

bcftools view \
        --regions-file ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.genome.withchrX.bed \
        --output-file ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.treemutkept.withchrX.vcf \
        ${VCF}.gz

VCF=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.treemutkept.withchrX.vcf
OUT=$(echo $VCF | sed 's/.vcf//' | sed 's/.gz//')
module load gatk/4.1.8.1
module load bcftools
info_fields=$(cat $VCF | grep "INFO=<" |  grep "Number=1" | sed 's/##INFO=<ID=//' | sed 's/,.*//' | tr -t "\n" " " | sed 's/ / -F /g' | awk 'BEGIN{printf "-F "}{print $0}' | sed 's/ -F $//')
genotype_fields=$(cat $VCF | grep "FORMAT=<" | sed 's/##FORMAT=<ID=//' | sed 's/,.*//' | \
        tr -t "\n" " " | sed 's/ / -GF /g' | awk 'BEGIN{printf "-GF "}{print $0}' | sed 's/ -GF $//')
#info_fields_perallele=$(cat $VCF | grep "INFO=<" | grep "Number=A" | sed 's/##INFO=<ID=//' | sed 's/,.*//' | \
#        tr -t "\n" " " | sed 's/ / -ASF /g' | awk 'BEGIN{printf "-ASF "}{print $0}' | sed 's/ -ASF $//')
gatk VariantsToTable \
        --add-output-vcf-command-line \
        --show-filtered \
        -V ${OUT}.vcf \
        -F CHROM -F POS -F REF -F ALT -F FILTER -F QUAL -F INFO \
        -F EVENTLENGTH \
        -F MULTI-ALLELIC \
        -F TRANSITION \
        -F HOM-REF \
        -F HET \
        -F HOM-VAR \
        -F VAR \
        -F NO-CALL \
        -F NCALLED \
        -F NSAMPLES \
        -F TYPE \
        ${info_fields} ${genotype_fields} \
        -O ${OUT}.table \
        -LE
grep -e "^CHROM" -e "^chrX" ${OUT}.table > ${OUT}.chrX.table

