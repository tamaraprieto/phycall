#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 8
#SBATCH -t 100:00:00
#SBATCH --mem 120G

source ReadConfig.sh $1

module purge
module load bcftools/1.9

PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
HEALTHY=$(head -n 1 ${ORIDIR}/${CONTROL})

mkdir -p ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMappingBootstrap

VCF="${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForCellPhy"

# Create the bed file with gene name for Signatures and everything else
bcftools query -H -f '%CHROM:%POS[ %AD]\n' \
         ${VCF}.vcf  | \
        sed 's/:AD//g' | sed 's/\[[0-9]*\]//g' | \
        sed 's/^# //' \
        > ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMappingBootstrap/AD.bed
# Create the bed file with Gene name
bcftools query -f '%CHROM %POS %POS %Func.refGene %Gene.refGene\n' \
        ${VCF}.vcf | \
        awk '{if ($4!="exonic"){$5="NA"}; print $0}' | \
        awk '{gensub("\\.*","",$5);print $0}' | \
        #sed 's/\\.*//' |
        awk '{$2=$2-1;print $0}' > \
        ${WORKDIR}/isec_somatic-germline/FILTERED/MutationMappingBootstrap/gene.bed

module purge
source deactivate
module load anaconda3/10.19
source deactivate
conda deactivate
conda activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
Rscript ${SCRIPTDIR}/R/MutationMappingCellPhyTreeMutBootstrap.R $PATIENT ${SLURM_ARRAY_TASK_ID}
conda deactivate
