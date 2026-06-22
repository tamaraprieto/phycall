#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 30:00:00
#SBATCH --mem 60G


source ReadConfig.sh $1

SAMPLE="Invitro_H1"
SAMPLETAG="H1"
SAMPLE="SCMDA.C6"
SAMPLETAG="C6"
DEPTH=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${WORKDIR}/ADO-${SAMPLETAG}/Depths.txt)

echo $SAMPLE
suffix=".recal"
myoptions=" -c 'chr21' -d "$DEPTH" -r"${RESDIR}/${REF}.fasta
export PATH=/gpfs/commons/groups/landau_lab/tprieto/apps/SingleCheck:$PATH

SingleCheckH1b \
        $myoptions \
        ${WORKDIR}/${SAMPLE}${suffix}.bam
