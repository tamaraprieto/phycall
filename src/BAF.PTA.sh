#!/bin/sh

source ReadConfig.sh $1
CONFIG=$1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
window=$2

REF="hg38"
MAINDIR="${WORKDIR}/BAF/"
PreprocessingResults="/gpfs/commons/groups/landau_lab/tprieto/jlaval/PreprocessingResults/"

# Make sure that nb_blocks = bin_size // block_size is an integer and that nb_blocks // (sliding_window-1) is an integer
#From the counts, generate the BAF profiles 
# with/without a sliding window and standard (no switches) or switches allowed. 
# The reasoning of not using sliding window and no switches is that the data is still good

## Option 1: No slide no switch (closest version to the data)
# Blocksize must divide binsize
#BLOCKSIZE=5000
#COUNTS_FOLDER="${MAINDIR}counts_5k/"
#SLIDING_NUMBER=2 #2 is for 0 slide
#MAX_SWITCH=0

# Option 2: Slide and switch
SLIDING_NUMBER=4  # 6 is the number was set as default
BLOCKSIZE=5000
COUNTS_FOLDER="${MAINDIR}counts_5k/"
MAX_SWITCH=1 # jean annotated as 1 but I think is zero?

#BAF_FOLDER="${MAINDIR}BAF_${window}/"
BAF_FOLDER="${MAINDIR}BAF_${window}_sliding${SLIDING_NUMBER}_switch${MAX_SWITCH}/"
mkdir $BAF_FOLDER


module purge
source deactivate
NAMES=${PreprocessingResults}cell_names_indexes_${PATIENT}_${REF}.tsv
SCRIPTS_FOLDER=/gpfs/commons/home/jquentin/tamara_project/jlaval/scripts/
#SCRIPTS_FOLDER=/gpfs/commons/groups/landau_lab/tprieto/jlaval/scripts/jean_scripts/
module load anaconda3/10.19
#source activate  /gpfs/commons/home/jquentin/.conda/envs/py27 # jeans
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/jeans_copy # jeans cloned
#source activate /gpfs/commons/groups/landau_lab/tprieto/conda/phython2.7.8
sbatch --array=1-22 ${SCRIPTS_FOLDER}BAF.sh $COUNTS_FOLDER $BAF_FOLDER $NAMES $window $BLOCKSIZE $SLIDING_NUMBER $MAX_SWITCH
source deactivate
