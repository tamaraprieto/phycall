#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 05:00:00
#SBATCH --mem 50G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//' | sed 's/E$//')

cd /gpfs/commons/groups/landau_lab/tprieto/apps/annovar

VCF=${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs-indels.vcf

# the output file will be named based on the prefix in out followed by "hg38_multianno.vcf"
./table_annovar.pl \
	${VCF} \
	humandb/ \
	-buildver hg38 \
	-out ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs-indels \
	-protocol refGene,dbnsfp42c,cosmic70,avsnp150,exac03,clinvar_20220320 \
	-remove \
	-operation g,f,f,f,f,f \
	-nastring . --vcfinput


# The -operation argument tells ANNOVAR which operations to use for each of the protocols: g means gene-based, gx means gene-based with cross-reference annotation (from -xref argument), r means region-based and f means filter-based.

module purge
module load bcftools/1.15.1
bcftools +fill-tags ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs-indels.hg38_multianno.vcf \
        --threads 4 \
        -Oz -o ${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs-indels.hg38_multianno.vaf.vcf \
        -- \
        --tags FORMAT/VAF

bcftools query -H -f '%CHROM:%POS %Gene.refGene %FILTER %REF %ALT %DP %AAChange.refGene %ExonicFunc.refGene[ %VAF][ %DP]\n' \
	${WORKDIR}/Mutect2.exome.${PATIENT}.SNVs-indels.hg38_multianno.vaf.vcf | \
	sed 's/\[[0-9]*\]//g' | \
	sed 's/_L003//g' | \
	sed 's/MRD-BALL-PTA-NEXTERA-EXOME-//g' | \
        sed 's/^# //' \
        > ${WORKDIR}/../SomaticCalls/Mutect${PATIENT}.annotations.bed

