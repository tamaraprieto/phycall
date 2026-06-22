#!/bin/sh

source ReadConfig.sh $1
CONFIG=$1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
window=$2
hap_window=$3

REF="hg38"
MAINDIR="${WORKDIR}/BAF/"
COUNTS_FOLDER="${MAINDIR}counts_5k/"
mkdir $MAINDIR
mkdir $COUNTS_FOLDER

CHISEL_FOLDER=${WORKDIR}/Chisel/
PreprocessingResults="/gpfs/commons/groups/landau_lab/tprieto/jlaval/PreprocessingResults/"
mkdir $PreprocessingResults

module purge
module load samtools
NAMES="cell_names_indexes_${PATIENT}_${REF}.tsv"
samtools view -H ${CHISEL_FOLDER}Chisel.BarcodedCells.${PATIENT}.bam | \
	grep ^@RG | sed 's/ID:CB:Z://' | sed 's/SM://' | sed 's/.dedup.bam//' | \
	sed 's/-0//' > ${PreprocessingResults}${NAMES}
COUNTS="counts_${PATIENT}_${REF}.tsv"
ln -s ${CHISEL_FOLDER}/rdr/rdr.tsv ${PreprocessingResults}${PATIENT}_${REF}_rdr_5kb.tsv
ln -s ${CHISEL_FOLDER}baf/baf.tsv ${PreprocessingResults}${COUNTS}
ln -s ${CHISEL_FOLDER}rdr/total.tsv ${PreprocessingResults}${PATIENT}_${REF}_total_reads.tsv

module load anaconda3/10.19
source activate  /gpfs/commons/home/jquentin/.conda/envs/py27 # jeans environment with the required packages
#source activate /gpfs/commons/groups/landau_lab/tprieto/conda/phython2.7.8
SCRIPTS_FOLDER=/gpfs/commons/home/jquentin/tamara_project/jlaval/scripts/
sbatch --array=1-22 $SCRIPTS_FOLDER/block_count.sh $CONFIG \
	5000 $COUNTS $NAMES $COUNTS_FOLDER $REF $PreprocessingResults
source deactivate
