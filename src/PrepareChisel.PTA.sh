#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user jlaval@nygenome.org
#SBATCH --mail-type END
#SBATCH --cpus-per-task 8
#SBATCH -t 100:00:00
#SBATCH --mem 160G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

module load anaconda2/4.3.1
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/chisel


#*.dedup.bam
chisel_prep ${WORKDIR}/*.dedup.bam \
        --output ${WORKDIR}/Chisel.BarcodedCells.${PATIENT}.bam \
        --noduplicates \
	-j 8

rm -r ${WORKDIR}/Chisel
mkdir -p ${WORKDIR}/Chisel
mv Chisel.BarcodedCells.${PATIENT}.*  ${WORKDIR}/Chisel
