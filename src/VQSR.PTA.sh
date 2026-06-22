#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 10:00:00
#SBATCH --mem 40G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

module purge
module load gatk/4.1.8.1
module load R/3.2.1
cd ${RESDIR}


# Separate indels and SNVs because SNVs should be recalibrated separately
gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${WORKDIR}/HaplotypeCaller.${PATIENT}.vcf \
         -O ${WORKDIR}/HaplotypeCaller.${PATIENT}.SNVs.vcf.gz \
         --select-type-to-include SNP

gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${WORKDIR}/HaplotypeCaller.${PATIENT}.vcf \
         -O ${WORKDIR}/HaplotypeCaller.${PATIENT}.indels.vcf.gz \
         --select-type-to-include INDEL



# The maximum DP (depth) filter only applies to whole genome data. But I am not sure it can be used for scWGA data
gatk VariantRecalibrator \
   -R ${RESDIR}/${REF}.fasta \
   -V  ${WORKDIR}/HaplotypeCaller.${PATIENT}.SNVs.vcf.gz  \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz \
   -an QD -an MQ \
   -an MQRankSum -an ReadPosRankSum \
   -an FS -an SOR \
   -mode SNP \
   -O ${WORKDIR}/VQSR.recal \
   --tranches-file ${WORKDIR}/VQSR.tranches \
   --rscript-file ${WORKDIR}/VQSR.plots.R 
#   -an DP \ Might not be a good idea using this parameter for WGA data Mar2022


gatk ApplyVQSR \
   -R ${RESDIR}/${REF}.fasta \
   -V ${WORKDIR}/HaplotypeCaller.${PATIENT}.SNVs.vcf.gz \
   -O ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.vcf.gz \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file ${WORKDIR}/VQSR.tranches \
   --recal-file ${WORKDIR}/VQSR.recal \
   -mode SNP


module load bcftools/1.9
zcat ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.vcf.gz | \
	awk '{if ($1~/^#CHR/){gsub(/_[A|T|C|G][A|C|G|T]+_\S+/,""); print}
		else{print $0}}' | \
	sed 's/417_B9_AGCGCTGTGT/417_B9/' | \
	bgzip -c > ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.reheaded.vcf.gz
tabix -p vcf ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.reheaded.vcf.gz


HEALTHY=$(head -n 1 ${ORIDIR}/${CONTROL})

# The step below I think I don't need it anymore after having called somatic variants with Mutect2
#When it fails the select commmand, it is ussually because the sample does not exists on the header
echo "gatk \
        SelectVariants \
        -R ${RESDIR}/${REF}.fasta \
        -V ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.reheaded.vcf.gz \
	 -O ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.somatic.vcf.gz \
	 --select-type-to-include SNP \
	 --remove-unused-alternates \
	 -select 'vc.getGenotype(\"$HEALTHY\").isHomRef() && vc.getGenotype(\"${HEALTHY}\").getDP() > 10'" | bash -

TUMOR=$(grep -E "(T1|4072T|380T|368B)" ${ORIDIR}${SAMPLELIST})
# Obtain a set of heterozygous SNPs to calculate allelic imbalance and dropout rates
echo "gatk \
         SelectVariants \
         -R ${RESDIR}/${REF}.fasta \
         -V ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.reheaded.vcf.gz \
         -O ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.heterozygous.vcf.gz \
         --select-type-to-include SNP \
         --remove-unused-alternates \
         --restrict-alleles-to BIALLELIC \
         -select 'vc.getGenotype(\"$HEALTHY\").isHet() && vc.getGenotype(\"${HEALTHY}\").getDP() > 15 && vc.getGenotype(\"$TUMOR\").isHet()'" | bash -


echo "gatk \
         SelectVariants \
         -R ${RESDIR}/${REF}.fasta \
	 -L chrX \
         -V ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.reheaded.vcf.gz \
         -O ${WORKDIR}/HaplotypeCaller.${PATIENT}.vqsr.hom.chrX.vcf.gz \
         --select-type-to-include SNP \
         --remove-unused-alternates \
         --restrict-alleles-to BIALLELIC \
         -select 'vc.getGenotype(\"$HEALTHY\").isHomVar() && vc.getGenotype(\"${HEALTHY}\").getDP() > 15 && vc.getGenotype(\"$HEALTHY\").getAD().0 == 0 && vc.getGenotype(\"$TUMOR\").isHomVar() && vc.getGenotype(\"$TUMOR\").getAD().0 == 0'" | bash -

