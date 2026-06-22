#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 8
#SBATCH -t 100:00:00
#SBATCH --mem 10G

source ReadConfig.sh $1
module purge
module load bcftools/1.9
module load bedtools/2.29.0
module load anaconda3/10.19
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
export PATH=/gpfs/commons/groups/landau_lab/tprieto/apps/cellphy-2022:$PATH
export PATH=/gpfs/commons/groups/landau_lab/tprieto/apps/angsd:$PATH
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

filename=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.vcf
# Annovar annotated file
filename=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.vcf

# Remove samples from the tree
bcftools view \
	-s ^4295_4272T1,4295_UndeterminedBarcode,4295_H3,417_417T,417_368B,445_445T,445_380T,4084_4072NonB,4084_4072T,InvitroBulk --force-samples \
	$filename > ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.vcf

# by default the PLs are used instead of the GTs. BE CAREFUL because GATK does not generate raw PLs but normalized PLs
MODEL=GTGTR4+G+FO
MODEL=GT16+FO+E
# On unphased genotype input data, which in fact only contains 10 states, the GT10 model appears to be as accurate as GT16 but requires only half of the time
MODEL=GT10+FO+E
mkdir -p ${WORKDIR}/CellPhy-03
rm ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf

if [ "$PATIENT" = "Invitro" ] || [ "$PATIENT" = "molm13" ]
	then
		rm ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.vcf
		ln -s ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.vcf ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf
else
	# Remove SNVs in CNVs
	echo "> Create the tsv tab with annotations"
	awk 'BEGIN{printf "#"}{print $0}' ${WORKDIR}/AllelicImbalance25.bed > \
        	 ${WORKDIR}/AI.annotations.tab
	bgzip -c ${WORKDIR}/AI.annotations.tab > \
        	${WORKDIR}/AI.annotations.tab.gz
	tabix -c 1 -f -b 2 -e 3 -0 \
        	--comment '#' \
        	${WORKDIR}/AI.annotations.tab.gz
	bcftools annotate \
       		--annotations ${WORKDIR}/AI.annotations.tab.gz \
       		-h <(echo '##INFO=<ID=allelicimbalance,Number=.,Type=String,Description="Allelic imbalance accross samples">\n') \
       		-c CHROM,FROM,TO,allelicimbalance \
       		--output ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.AI.vcf \
       		${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.vcf

	bcftools filter \
       		--include "allelicimbalance=='diploid' || allelicimbalance=='diploidAI2cells' || allelicimbalance=='diploidAI1cell'"  \
        	${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.AI.vcf  \
        	> ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.vcf
		ln -s ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.vcf ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf
fi


# Remove potential false positives in patient males
# If I remove duplicates they are not going to be there because chrX is not diploid
#if [ "$GENDER" = "XX" ]
#	then
#		ln -s ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.vcf ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf
#	else
#		a=$(grep "^chrom" ${WORKDIR}/isec_somatic-germline/FILTERED/SC.annotations.bed | sed 's/#//' | \
#			sed 's/\t/\n/g' | head -n 5 | tail -n 1)
#	bcftools filter \
#        	--exclude "CHROM=='chrX' || CHROM=='chrY' && ${a}<1" \
#        	${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.vcf \
#        	> ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf
#fi

#exit

# The PL provided by HaplotypeCaller is not the PL expected by CellPhy. It is the SCcaller after modification
# If I don't do the prob-msa off then the program is going to use PL
cellphy.sh RAXML \
	--msa ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy.vcf \
	--model ${MODEL} \
	--threads 8 \
	--all \
	--redo \
	--prefix ${WORKDIR}/CellPhy-03/CellPhy.${MODEL}.nobulks.noCNVs


#        --msa ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.vcf \

#	--prefix ${WORKDIR}/CellPhy-02/CellPhy.msaoff.${MODEL}.nobulks \
#	--prob-msa off \
#	--bs-trees 100


#        --bs-trees 100 \ # by default cellphy will determine the number of bootstraps, max 1000. I can add this line to specify it
#        --prefix ${WORKDIR}/CellPhy-02/CellPhy.msaoff.${MODEL}.nobulks \
#        --prob-msa off
#         --prefix ${WORKDIR}/CellPhy-02/CellPhy.${MODEL}.nobulks

