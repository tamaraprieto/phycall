#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 10:00:00
#SBATCH --mem 40G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
module purge
module load bedtools/2.29.0
module load bcftools/1.9
module load vcftools/0.1.17
CONTEXT_LEN="1"
REF_GENOME=${RESDIR}/${REF}.fasta
SUFFIX=".internalbranch"
SUFFIX=".externalbranch"

VCF=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.hg38_multianno.treemutkept${SUFFIX}.vcf

# check the chromosome format in REF_GENOME
chr_name=`head -1 $REF_GENOME | awk '{print $1}'`

mkdir -p ${WORKDIR}/isec_somatic-germline/persample${SUFFIX}
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})
bcftools view -s $SAMPLE \
	 ${VCF} | \
	 grep -v '0/0:' | \
	 grep -v '0|0:' | \
	 grep -vE "\.\/\.:" | \
	grep -vE "\.\|\.:" \
	 > \
	 ${WORKDIR}/isec_somatic-germline/persample${SUFFIX}/${SAMPLE}.vcf

VCF=${WORKDIR}/isec_somatic-germline/persample${SUFFIX}/${SAMPLE}.vcf

echo $SAMPLE

    if [[ $chr_name == *"chr"* ]]
    then
    	    bcftools view -f 'PASS' --types snps --max-alleles 2 $VCF | \
            		    grep  -v "#" | awk -v var=$CONTEXT_LEN \
                        '{OFS="\t"; print $1, $2-1-var, $2+var}' | \
                            bedtools getfasta -fi $REF_GENOME -bed - | \
                            grep -v "^>" | awk '{print toupper($0)}' > $SAMPLE".tmp"
    else
            bcftools view -f 'PASS' --types snps --max-alleles 2 $VCF | \
			grep  -v "#" | sed "s/^chr//g" | awk -v var=$CONTEXT_LEN \
                        '{OFS="\t"; print $1, $2-1-var, $2+var}' | \
                            bedtools getfasta -fi $REF_GENOME -bed - | \
                            grep -v "^>" | awk '{print toupper($0)}' > $SAMPLE".tmp"

    fi
        bcftools view -f 'PASS' --types snps --max-alleles 2 $VCF | grep  -v "#" | \
	awk -v var=$SAMPLE '{print $5"\t"var}' | \
	paste -d"\t" $SAMPLE".tmp" - | \
	awk '{A=substr($1,1,1);B=substr($1,2,1);C=substr($1,3,1); print A"\t"B"\t"C"\t"$2"\t"$3}' | \
	awk -v OFS='\t' \
	'{if ($2=="A" || $2=="G")
	{   for (i = 1; i <= 4; i++){
	       if ($i=="A"){$i="T"} else if ($i=="T"){$i="A"} else if ($i=="C"){$i="G"} else if ($i=="G"){$i="C"}}
	print $3$2">"$4$1"\t"$5 
	}
	else{print $1$2">"$4$3"\t"$5}}' | \
	sort -k2,2 -k1,1 | \
	uniq -c | awk '{print $2"\t"$3"\t"$1}' > ${WORKDIR}/96MutPat.${SAMPLE}.Unspread${SUFFIX}.txt

rm $SAMPLE".tmp"

