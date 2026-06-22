#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 20:00:00
#SBATCH --mem 100G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

HEALTHY=$(head -n 1 ${ORIDIR}/${CONTROL})
TUMOR=$(grep "T1" ${ORIDIR}${SAMPLELIST})

module purge
module load gatk/4.1.8.1
module load bcftools

NEWDIR=${WORKDIR}/isec_somatic-germline/groundtruth/
rm -r $NEWDIR
mkdir $NEWDIR

# TP somatic autosomes
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${WORKDIR}/isec_somatic-germline/0003.vcf.gz \
        -O ${NEWDIR}/TP.vcf.gz \
        -select 'vc.getGenotype(\"$HEALTHY\").isHomRef() && vc.getGenotype(\"${TUMOR}\").isHet() && vc.getGenotype(\"$TUMOR\").getDP()  > 9'" | bash -


zgrep -v "^chrX" ${NEWDIR}/TP.vcf.gz | grep -v "^chrY" | gzip -c > ${NEWDIR}/TP.autosomes.vcf.gz
tabix -p vcf ${NEWDIR}/TP.autosomes.vcf.gz
rm ${NEWDIR}/TP.vcf.gz*

# TPs somatic X chromosome
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${WORKDIR}/isec_somatic-germline/0003.vcf.gz \
        -O ${NEWDIR}/TP.X.vcf.gz \
	-L chrX \
	-select 'vc.getGenotype(\"$HEALTHY\").isHomRef() && vc.getGenotype(\"${TUMOR}\").isHomVar() && vc.getGenotype(\"$TUMOR\").getDP()  > 4'" | bash -

# FPs
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${WORKDIR}/isec_somatic-germline/0003.vcf.gz \
        -O ${NEWDIR}/FP.X.vcf.gz \
        -L chrX \
        -select 'vc.getGenotype(\"$HEALTHY\").isHomRef() && vc.getGenotype(\"${TUMOR}\").isHomRef() && vc.getGenotype(\"$TUMOR\").getDP()  > 4 && vc.getHetCount() > 2'" | bash -

# Create the bed files and a combined file
zgrep -v "^#" ${NEWDIR}/TP.autosomes.vcf.gz | \
	awk '{print $1"\t"$2-1"\t"$2"\tTP"}' > ${NEWDIR}/TP.autosomes.bed

zgrep -v "^#" ${NEWDIR}/FP.X.vcf.gz | \
        awk '{print $1"\t"$2-1"\t"$2"\tFP"}' > ${NEWDIR}/FP.X.bed

cat ${NEWDIR}/TP.autosomes.bed ${NEWDIR}/FP.X.bed | sort -k 1,1 -k2,2n > ${NEWDIR}/TP-FP.bed 

# Intersect the bed file with the HC VCF file to obtain the informations

num_fields_vcf=$(zgrep "^#CHR" ${WORKDIR}/isec_somatic-germline/0003.vcf.gz | awk '{print NF}')
module load bedtools
bedtools intersect -wao -header \
	-b ${NEWDIR}/TP-FP.bed \
	-a ${WORKDIR}/isec_somatic-germline/0003.vcf.gz | \
	awk -v numf=$num_fields_vcf '{if ($NF==0){}else if($0~/^#/){print $0} else{for (i=1;i<=6;i++){printf $i"\t"}; printf $(NF-1)"\t"; for (i=8;i<=numf;i++){printf $i"\t"}; printf "\n"}}' | sed 's/\t\n/\n/' | bgzip -c > ${NEWDIR}/TP-FP.vcf.gz
tabix -p vcf ${NEWDIR}/TP-FP.vcf.gz

# Create a file with AD counts and INFO for the complex heatmap
module load vcftools 
vcftools --gzvcf ${NEWDIR}/TP-FP.vcf.gz --extract-FORMAT-info AD --out ${NEWDIR}/TP-FP
info_params=$(zcat ${NEWDIR}/TP-FP.vcf.gz | grep INFO= | sed 's/,.*//' | sed 's/.*=//' | awk '{print "--get-INFO "$1}' | tr -s "\n" "\t")
zgrep -v '^##' ${NEWDIR}/TP-FP.vcf.gz | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' | sed 's/^#//' > ${NEWDIR}/TP-FP.INFO
vcftools --gzvcf ${NEWDIR}/TP-FP.vcf.gz $info_params --out ${NEWDIR}/TP-FP.ALL.HC

# Intersect the bed file with the Mutect2 VCF file to obtain the informations

num_fields_vcf=$(zgrep "^#CHR" ${WORKDIR}/isec_somatic-germline/0002.vcf.gz | awk '{print NF}')
module load bedtools
bedtools intersect -wao -header \
        -b ${NEWDIR}/TP-FP.bed \
        -a ${WORKDIR}/isec_somatic-germline/0002.vcf.gz | \
        awk -v numf=$num_fields_vcf '{if ($NF==0){}else if($0~/^#/){print $0} else{for (i=1;i<=numf;i++){printf $i"\t"}; printf "\n"}}' | sed 's/\t\n/\n/' | bgzip -c > ${NEWDIR}/TP-FP.mutect.vcf.gz
tabix -p vcf ${NEWDIR}/TP-FP.mutect.vcf.gz

# Create a file with AD counts and INFO for the complex heatmap
module load vcftools
info_params=$(zcat ${NEWDIR}/TP-FP.mutect.vcf.gz | grep INFO= | sed 's/,.*//' | sed 's/.*=//' | awk '{print "--get-INFO "$1}' | tr -s "\n" "\t")
vcftools --gzvcf ${NEWDIR}/TP-FP.mutect.vcf.gz $info_params --out ${NEWDIR}/TP-FP.ALL.mutect

