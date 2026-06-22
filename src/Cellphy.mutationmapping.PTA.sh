#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 01:00:00
#SBATCH --mem 40G


source ReadConfig.sh $1
module purge
module load bcftools/1.9
module load bedtools/2.29.0
module load anaconda3/10.19
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
export PATH=/gpfs/commons/groups/landau_lab/tprieto/apps/cellphy-2022:$PATH
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

MODEL=GTGTR4+G+FO
MODEL=GT16+FO+E
# On unphased genotype input data, which in fact only contains 10 states, the GT10 model appears to be as accurate as GT16 but requires only half of the time
MODEL="GT10+FO+E"


####################
# MAP GENOME CALLS #
####################

# All mutations, even those present in CNVs (but mapping to the tree built without CNVs)
cellphy.sh RAXML \
        --mutmap \
        --msa  ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.vcf \
        --model ${WORKDIR}/CellPhy-02/CellPhy.${MODEL}.nobulks.noCNVs.raxml.bestModel \
        --tree ${WORKDIR}/CellPhy-02/CellPhy.${MODEL}.nobulks.noCNVs.raxml.bestTree \
        --opt-branches off \
        --prefix  ${WORKDIR}/CellPhy-02/CellPhy.${MODEL}.MappedGenome \
        --threads 1 \
        --redo
# Create the bed file with Gene name
bcftools query -f '%CHROM %POS %POS %Func.refGene %Gene.refGene\n' \
	${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.vcf | \
	awk '{if ($4!="exonic"){$5="NA"}; print $0}' | \
	awk '{gensub("\\.*","",$5);print $0}' | \
	#sed 's/\\.*//' |
	awk '{$2=$2-1;print $0}' > \
	${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.bed

# Create a bed file with mutational pattern (only those not present in CNVs)
bcftools view ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.vcf | \
            grep  -v "#" | awk -v var=$CONTEXT_LEN \
            '{OFS="\t"; print $1, $2-1-var, $2+var}' | \
                bedtools getfasta -fi ${RESDIR}/${REF}.fasta -bed - | \
                grep -v "^>" | awk '{print toupper($0)}' > \
	${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.mutationalpattern.bed.2
bcftools query -f '%CHROM:%POS\t%ALT\n' \
        ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.vcf | \
        paste ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.mutationalpattern.bed.2 - | \
        awk '{A=substr($1,1,1);B=substr($1,2,1);C=substr($1,3,1); print A"\t"B"\t"C"\t"$3"\t"$2}' | \
        awk -v OFS='\t' \
        '{if ($2=="A" || $2=="G")
        {   for (i = 1; i <= 4; i++){
               if ($i=="A"){$i="T"} else if ($i=="T"){$i="A"} else if ($i=="C"){$i="G"} else if ($i=="G"){$i="C"}}
        print $5"\t"$3$2">"$4$1
        }
        else{print $5"\t"$1$2">"$4$3}}' > \
        ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.mutationalpattern.bed
rm ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.mutationalpattern.bed.2

#Obtain VAF from the VCF (compare with the mapping)
grep "^#CHR" ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.vcf | cut -f1,2,10- | \
        sed 's/#//' >  ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.VAF-DP.txt
bcftools query \
        -f '%CHROM %POS[\t%VAF;%DP]\n' ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.vcf >> \
        ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.noCNV.VAF-DP.txt

###################
# MAP EXOME CALLS #
###################

# I want to rescue as many mutations as possible even if they are present in CNVs
#module load bedtools/2.25.0
# Create the annotation file
bedtools intersect -header \
        -a ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.nobulks.vcf \
        -b ${WORKDIR}/../../../2021-07-19_BioSkryb/SomaticCalls/${PATIENT}_SNP_filtered.bed \
        > ${WORKDIR}/MutationsToMap.${PATIENT}.vcf

# Annotate the mutation on the tree
cellphy.sh RAXML \
	--mutmap \
	--msa ${WORKDIR}/MutationsToMap.${PATIENT}.vcf \
	--model ${WORKDIR}/CellPhy-02/CellPhy.${MODEL}.nobulks.noCNVs.raxml.bestModel \
	--tree ${WORKDIR}/CellPhy-02/CellPhy.${MODEL}.nobulks.noCNVs.raxml.bestTree \
	--opt-branches off \
	--prefix  ${WORKDIR}/CellPhy-02/CellPhy.${MODEL}.MappedExome \
	--threads 1 \
	--redo

#Obtain VAF from the VCF (compare with the mapping)
grep "^#CHR" ${WORKDIR}/MutationsToMap.${PATIENT}.vcf | cut -f1,2,10- | \
	sed 's/#//' >  ${WORKDIR}/MutationsToMap.${PATIENT}.VAF-DP.txt
bcftools query \
	-f '%CHROM %POS[\t%VAF;%DP]\n' ${WORKDIR}/MutationsToMap.${PATIENT}.vcf >> \
	${WORKDIR}/MutationsToMap.${PATIENT}.VAF-DP.txt 
