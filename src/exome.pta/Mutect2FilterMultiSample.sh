#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 01:00:00
#SBATCH --mem 40G


source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
echo "Patient: "$PATIENT

# Chromosomes: use all independently of the ploidy 
num=24

# COMBINE ALL VCFs FROM MUTECT
mutect_files=$(grep -E "^chr([0-9]*|X|Y)\s" ${RESDIR}/${REF}.fasta.fai | \
        awk '{print $1}' | sort --version-sort | head -n $num | \
	awk -v dir=$WORKDIR -v patient=$PATIENT \
	'{print "--variant "dir"/Mutect2."patient"."$1".vcf"}' | tr '\n' ' ')

module purge
module load gatk/3.8.1
java -cp /nfs/sw/gatk/gatk-3.8.1/GenomeAnalysisTK.jar \
	org.broadinstitute.gatk.tools.CatVariants \
	-R ${RESDIR}/${REF}.fasta \
	${mutect_files} \
	-out ${WORKDIR}/Mutect2.${PATIENT}.unfiltered.vcf \
	--assumeSorted

module purge
module load gatk/4.1.8.1
# COMBINE ALL f1r2 FILES FROM MUTECT (one per chromosome)
all_f1r2_input=$( grep -E "^chr([0-9]*|X|Y)\s" ${RESDIR}/${REF}.fasta.fai | \
        awk '{print $1}' | sort --version-sort | head -n $num | \
        awk -v dir=$WORKDIR -v patient=$PATIENT \
        '{print "--input "dir"/f1r2.mutect2."patient"."$1".tar.gz"}' | tr '\n' ' ')
gatk LearnReadOrientationModel \
        $all_f1r2_input \
        -O ${WORKDIR}/Read-orientation-model.${PATIENT}.tar.gz



# COMBINE ALL STATS FILES FROM MUTECT
STATS=$( grep -E "^chr([0-9]*|X|Y)\s" ${RESDIR}/${REF}.fasta.fai | \
        awk '{print $1}' | sort --version-sort | head -n $num | \
	awk -v dir=$WORKDIR -v patient=$PATIENT \
	'{print "-stats "dir"/Mutect2."patient"."$1".vcf.stats"}' | tr '\n' ' ')
gatk MergeMutectStats \
        ${STATS} \
        -O ${WORKDIR}/Mutect2.${PATIENT}.merged.stats

# CONTAMINATION FILES
CONTAMS= SAMPLES=$(grep -v "NonB"  ${WORKDIR}/../${IDLIST} |  awk -v dir=$WORKDIR '{print "--contamination-table "dir"/"$0".matched.contamination.table"}' | tr '\n' ' ')

# APPLY THE FILTER
gatk \
	FilterMutectCalls \
	-R ${RESDIR}/${REF}.fasta \
	-V ${WORKDIR}/Mutect2.${PATIENT}.unfiltered.vcf \
	-stats ${WORKDIR}/Mutect2.${PATIENT}.merged.stats \
	--ob-priors ${WORKDIR}/Read-orientation-model.${PATIENT}.tar.gz \
	${CONTAMS} \
	-O ${WORKDIR}/Mutect2.${PATIENT}.FILTERS.vcf


# REMOVE INDELS AND NON PASS
module purge
module load gatk/4.1.8.1
gatk SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${WORKDIR}/Mutect2.${PATIENT}.FILTERS.vcf \
        -O ${WORKDIR}/Mutect2.${PATIENT}.SNVs.vcf \
        --exclude-filtered true \
        --select-type-to-include SNP

