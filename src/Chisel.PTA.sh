#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 8
#SBATCH --mem-per-cpu 20G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
module purge
#module load chisel
module load anaconda2/4.3.1
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/chisel

#variables
window=$2
hap_window=$3

#for 4295
#hg38
#BARCODED='/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/4295/RESULTS/CHISEL/input/'
#hg19
#BARCODED="/gpfs/commons/groups/landau_lab/tprieto/BioSkryb/data/WGS/4295/RESULTS_hg19/"
#for the others
#BARCODED="/gpfs/commons/home/jlaval/PreprocessingResults/BARCODED/${PATIENT}/"

HEALTHY=$(sed "1q;d" $ORIDIR/$CONTROL)
#HEALTHY=$(sed "1q;d" ${ORIDIR}/HealthyCellBadQual.txt)

#Use chisel command
chisel  -t ${WORKDIR}/Chisel/Chisel.BarcodedCells.${PATIENT}.bam \
        -n ${WORKDIR}/${HEALTHY}.recal.bam \
	--chromosomes "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22" \
        -r ${RESDIR}/${REF}.fasta \
	-l ${WORKDIR}/SHAPEIT/Phased.${PATIENT}.healthy.hetSNPs.vcf \
	--rundir ${WORKDIR}/Chisel \
        --seed 12 \
        --maxploidy 3 \
        --jobs 8 \
        --size $window \
        --blocksize $hap_window # We use blocks of length 50kb in this work as phasing errors are unlikely at this scale

#--chromosomes "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22" \ hg38
#--chromosomes "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22" \ hg19
#phased VCF for 417, 445, 4084
#-l ${WORKDIR}/Phased.${PATIENT}.hg38.edited.vcf \
#hg19
#-l ${WORKDIR}/chrs.dose.edited.vcf \

#mv baf ${WORKDIR}/Chisel
#mv clones ${WORKDIR}/Chisel
#mv combo ${WORKDIR}/Chisel
#mv plots ${WORKDIR}/Chisel
#mv rdr ${WORKDIR}/Chisel
#mv calls ${WORKDIR}/Chisel
