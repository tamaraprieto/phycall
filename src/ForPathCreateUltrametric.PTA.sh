#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 4
#SBATCH -t 100:00:00
#SBATCH --mem 1G

source ReadConfig.sh $1
indir=${WORKDIR}/isec_somatic-germline/FILTERED/MutationMappingBootstrap
input=${indir}/TreeMutWithZeros.${SLURM_ARRAY_TASK_ID}.rds
output=${indir}/TreeMutWithZeros.${SLURM_ARRAY_TASK_ID}.ultrametric.rds

module purge
source deactivate
module load anaconda3/10.19
source deactivate
conda deactivate
conda activate /gpfs/commons/groups/landau_lab/tprieto/conda/myR4
# original script from Sri located at /gpfs/commons/home/srajagopalan/run_treemut.R
Rscript /gpfs/commons/home/tprieto/myscripts/DNAphy/src/R/ForPathCreateUltrametric.R $input $output
conda deactivate
