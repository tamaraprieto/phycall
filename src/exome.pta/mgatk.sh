#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 12
#SBATCH -t 30:00:00
#SBATCH --mem 30G

source ReadConfig.sh $1

PATIENT=$(basename $1 | sed 's/Config.//')

module purge
module load java/1.9
module load miniconda2/4.4.10
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/mgatk

#ORIDIR=${ORIDIR}Bams/ # 4295 only

mgatk call -i ${ORIDIR} \
	--mito-genome 'hg38' \
	--output ${WORKDIR}/mgatk \
	--ncores 12 \
	--alignment-quality 20 \
	--base-qual 20 \
	--emit-base-qualities

source deactivate
