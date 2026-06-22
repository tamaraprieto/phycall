#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 15:00:00
#SBATCH --mem 30G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

GINKGO_INST_DIR=/gpfs/commons/groups/landau_lab/tprieto/apps/ginkgo/
# Create sample list file
ls ${GINKGO_INST_DIR}uploads/${PATIENT} | grep .bed.gz$ > ${GINKGO_INST_DIR}uploads/${PATIENT}/list
# Create configuration file (default) 500kb
# cp $POSADALAB/APPS/ginkgo/config.example $POSADALAB/APPS/ginkgo/uploads/${PATIENT}/config
cp ${GINKGO_INST_DIR}config.10kb.example ${GINKGO_INST_DIR}uploads/${PATIENT}/config

module purge
module load anaconda3/10.19
module load php/8.1.0 
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4

# the first argument passed to analyze.sh has to be the name of the folder in ginkgo uploads
bash ${GINKGO_INST_DIR}scripts/analyze.sh ${PATIENT}

source deactivate
