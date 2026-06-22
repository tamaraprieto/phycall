#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem 40G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//')

################
# LOAD MODULES #
################

module purge
module load cellphy # my own module
module load R/4.0.2
module load bcftools/1.9
module load vcftools/0.1.14 # 1.17 version not properly installed I think 
module load anaconda3/10.19
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4

##############################################
# CREATE A FASTA MULTIPLE SEQUENCE ALIGNMENT #
##############################################

VCF=sc_covered.recode
VCF=4295_somatic.unidentified.recode # 51 somatic SNVs
VCF=0829-Tranche99.9_4295_Exome_SNP_ScHet.recode # 12 somatic SNVs email 

#  Previous file: ${ORIDIR}/filtered/${VCF}.vcf
bcftools view --max-alleles 2 \
	--samples-file ${ORIDIR}/../SC.${PATIENT} \
	--types snps \
	-e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
	${ORIDIR}/../SomaticCalls/${VCF}.vcf > ${WORKDIR}/${VCF}.vcf

bgzip -c ${WORKDIR}/${VCF}.vcf > ${WORKDIR}/${VCF}.vcf.gz 
tabix -p vcf ${WORKDIR}/${VCF}.vcf.gz
# vcf-to-tab is a vcftools function
zcat ${WORKDIR}/${VCF}.vcf.gz | vcf-to-tab > ${WORKDIR}/${VCF}.tab

###############
# RUN CELLPHY #
############### 

MODEL=GTGTR4+G+FO
mkdir ${WORKDIR}/CellPhy

/gpfs/commons/groups/landau_lab/tprieto/apps/cellphy/bin/raxml-ng-cellphy-linux \
                --all \
                --threads 1 \
                --redo \
		--msa ${WORKDIR}/${VCF}.vcf \
		--msa-format VCF \
                --model ${MODEL} \
                --prob-msa off \
                --prefix ${WORKDIR}/CellPhy/CellPhy.${MODEL}.GT.$VCF
