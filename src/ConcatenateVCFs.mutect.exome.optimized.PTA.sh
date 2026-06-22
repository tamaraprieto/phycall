#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 05:00:00
#SBATCH --mem 15G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//' | sed 's/E//')

number_intervals=$2
module load bcftools/1.9
ulimit -s 65536

# SNVs and indels
VCFs=$(seq 1 $number_intervals | \
        awk -v dir=$WORKDIR \
        -v patient=$PATIENT \
        '{print " "dir"/Mutect2.exome."patient"."$0".FILTERS.vcf.gz"}' \
        | tr '\n' ' ')
bcftools concat -a -o ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs-indels.vcf $VCFs
bcftools view -f PASS ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs-indels.vcf > ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs-indels.PASS.vcf

# SNVs
VCFs=$(seq 1 $number_intervals | \
        awk -v dir=$WORKDIR \
        -v patient=$PATIENT \
        '{print " "dir"/Mutect2.exome."patient"."$0".SNVs.vcf.gz"}' \
        | tr '\n' ' ')
bcftools concat -a -o ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs.vcf $VCFs

module load bcftools/1.9
cat ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs.vcf | \
        awk '{if ($1~/^#CHR/){gsub(/_[A|T|C|G][A|C|G|T]+_\S+/,""); print}
                else{print $0}}' | \
        bgzip -c > ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs.reheaded.vcf.gz
tabix -p vcf ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs.reheaded.vcf.gz
bcftools view -f PASS ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs.reheaded.vcf.gz > ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs.PASS.vcf

