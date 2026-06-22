#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu 30G

source ReadConfig.sh $1
CONFIG=$1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')
window=$2
REF="hg38"
MAINDIR="${WORKDIR}/BAF/"

## Option 1: No slide no switch (closest version to the data)
#SLIDING_NUMBER=2 #2 is for 0 slide
#MAX_SWITCH=0
# Option 2: Slide and switch
SLIDING_NUMBER=4  # 6 is the number was set as default
MAX_SWITCH=1 # jean annotated as 1 but I think is zero?


#BAF_FOLDER="${MAINDIR}BAF_${window}/"
BAF_FOLDER="${MAINDIR}BAF_${window}_sliding${SLIDING_NUMBER}_switch${MAX_SWITCH}/"
PreprocessingResults="/gpfs/commons/groups/landau_lab/tprieto/jlaval/PreprocessingResults/"

module purge
module load anaconda3/10.19
#source activate  /gpfs/commons/home/jquentin/.conda/envs/py27 # jeans
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/jeans_copy
#source activate /gpfs/commons/groups/landau_lab/tprieto/conda/phython2.7.8
#/gpfs/commons/groups/landau_lab/tprieto/conda/jeans_copy/bin/python ${SCRIPTDIR}/BAF.LOH.py $PATIENT $REF $MAINDIR $PreprocessingResults $window $SLIDING_NUMBER $MAX_SWITCH
/gpfs/commons/groups/landau_lab/tprieto/conda/jeans_copy/bin/python ${SCRIPTDIR}/BAF.LOH.reduceswitches.py $PATIENT $REF $MAINDIR $PreprocessingResults $window $SLIDING_NUMBER $MAX_SWITCH
source deactivate
