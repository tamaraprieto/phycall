#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 10:00:00
#SBATCH --mem 40G

source ReadConfig.sh $1
module purge
module load samtools/1.9
module load bcftools/1.9
module load miniconda2/4.4.10
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/vep

VCF=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.vcf
VCF=${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.unfiltered4.vcf

mkdir -p ${WORKDIR}/isec_somatic-germline/vcf2maf
echo $SAMPLE
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})
bcftools view -s $SAMPLE \
	 ${VCF} | \
        sed 's/^chr//' | \
        sed 's/ID=chr/ID=/' \
	> \
	${WORKDIR}/isec_somatic-germline/vcf2maf/${SAMPLE}.vcf

VCF=${WORKDIR}/isec_somatic-germline/vf2maf/${SAMPLE}.vcf

perl /gpfs/commons/groups/landau_lab/tprieto/apps/vcf2maf-1.6.21/vcf2maf.pl \
	--input-vcf ${WORKDIR}/isec_somatic-germline/vcf2maf/${SAMPLE}.vcf \
	--output-maf ${WORKDIR}/isec_somatic-germline/vcf2maf/${SAMPLE}.maf \
	--ncbi-build GRCh38 \
	--vep-overwrite \
	--ref-fasta /gpfs/commons/home/tprieto/.vep/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
	--vep-path /gpfs/commons/groups/landau_lab/tprieto/conda/vep/bin \
	--tumor-id $SAMPLE
