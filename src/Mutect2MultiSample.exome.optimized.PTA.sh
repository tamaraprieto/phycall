#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 40G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//' | sed 's/E//')
echo "Patient: "$PATIENT

#CHR=$(awk '{print $1}' ${RESDIR}/${REF}.fasta.fai | grep -E "chr[0-9|X|Y]+$" | sed "${SLURM_ARRAY_TASK_ID}q;d")

#grep "^${CHR}\s" \
#	/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/Illumina_Exome_TargetedRegions_v1.2.hg38.bed > \
#	/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/Illumina_Exome_TargetedRegions_v1.2.hg38.${CHR}.bed 


int_end=$(awk -v jobid=${SLURM_ARRAY_TASK_ID} 'BEGIN{print jobid*50}')
int_start=$(awk -v end=${int_end} 'BEGIN{print end-49}')

awk -v A=$int_start -v B=$int_end '{if (NR>=A && NR<=B){print $0}}' /gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/Illumina_Exome_TargetedRegions_v1.2.hg38.bed > /gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/Illumina_Exome_TargetedRegions_v1.2.hg38.${SLURM_ARRAY_TASK_ID}.bed


# SELECT NORMAL SAMPLE
NORMAL_NAME=$(head -n 1 ${ORIDIR}/${CONTROL})

module purge
module load bedtools/2.29.0
module load gatk/4.1.8.1

echo "> Run Mutect for chr "${SLURM_ARRAY_TASK_ID}
SAMPLES=$(awk -v dir=$WORKDIR '{print "-I "dir"/"$0".bqsr.marked.bam"}' ${ORIDIR}/${SAMPLELIST} | tr '\n' ' ')
#echo $SAMPLES
echo ""
echo ""
gatk --java-options "-Xmx40G" Mutect2 \
	-R ${RESDIR}/${REF}.fasta \
	${SAMPLES} \
        -normal ${NORMAL_NAME} \
	--germline-resource ${RESDIR}/bestpractices/af-only-gnomad.hg38.vcf.gz \
	-L /gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/Illumina_Exome_TargetedRegions_v1.2.hg38.${SLURM_ARRAY_TASK_ID}.bed \
	-O ${WORKDIR}/Mutect2.exome.${PATIENT}.${SLURM_ARRAY_TASK_ID}.vcf \
	--f1r2-tar-gz ${WORKDIR}/f1r2.mutect2.exome.${PATIENT}.${SLURM_ARRAY_TASK_ID}.tar.gz

rm /gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/Illumina_Exome_TargetedRegions_v1.2.hg38.${SLURM_ARRAY_TASK_ID}.bed

# -L /gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/2021-07-19_BioSkryb/Illumina_Exome_TargetedRegions_v1.2.hg38.${CHR}.bed \
