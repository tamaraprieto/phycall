#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 05:00:00
#SBATCH --mem 200G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

number_intervals=$(wc -l ${RESDIR}/hg38.even.handcurated.20k.intervals | awk '{print $1}')
VCFs=$(seq 1 $number_intervals | \
        awk -v dir=$WORKDIR \
        -v patient=$PATIENT \
        '{print " "dir"/GenotypedGVCFs."patient"."$0".vcf"}' \
        | tr '\n' ' ')

module load bcftools/1.9
ulimit -s 65536
bcftools concat -o ${WORKDIR}/HaplotypeCaller.${PATIENT}.vcf $VCFs

## This code below was to see if it was faster than merging so many VCF file but
## I haven't try it in the end
#VCFs_gatk=$(seq 1 $number_intervals | \
#        awk -v dir=$WORKDIR \
#        -v patient=$PATIENT \
#        '{print -I" "dir"/GenotypedGVCFs."patient"."$0".vcf"}' \
#        | tr '\n' ' ')
#module purge
#module load gatk/4.1.8.1
#gatk MergeVcfs \
#   $VCFs_gatk \
#   -O ${WORKDIR}/HaplotypeCaller.${PATIENT}.vcf
