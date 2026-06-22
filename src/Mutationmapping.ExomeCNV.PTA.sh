#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 01:00:00
#SBATCH --mem-per-cpu 1G


source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

module purge
source deactivate
module load anaconda3/10.19
source deactivate
conda deactivate
conda activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
Rscript ${SCRIPTDIR}/R/ExtraFilterExome.R $PATIENT
conda deactivate
