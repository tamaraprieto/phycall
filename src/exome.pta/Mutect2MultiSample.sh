#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 80G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//')
echo "Patient: "$PATIENT

CHR=$(grep -E "^chr([0-9]*|X|Y)\s" ${RESDIR}/${REF}.fasta.fai | \
        awk '{print $1}' | sort --version-sort | \
	sed "${SLURM_ARRAY_TASK_ID}q;d")

# SELECT NORMAL SAMPLE
NORMAL_NAME=$(grep "NonB" ${WORKDIR}/../${IDLIST})

module purge
module load bedtools/2.29.0
module load gatk/4.1.8.1

#for CHR in ${CHRS[@]}
#do
echo "> Run Mutect for chromosome "${CHR}
SAMPLES=$(awk -v dir=$ORIDIR '{print "-I "dir"/Bams/"$0".bqsr.marked.bam"}' ${WORKDIR}/../${IDLIST} | tr '\n' ' ')
echo $SAMPLES
gatk --java-options "-Xmx4G" Mutect2 \
	-R ${RESDIR}/${REF}.fasta \
	${SAMPLES} \
        -normal ${NORMAL_NAME} \
	--germline-resource ${RESDIR}/bestpractices/af-only-gnomad.hg38.vcf.gz \
	-L $CHR \
	-O ${WORKDIR}/Mutect2.${PATIENT}.${CHR}.vcf \
	--f1r2-tar-gz ${WORKDIR}/f1r2.mutect2.${PATIENT}.${CHR}.tar.gz
#done
