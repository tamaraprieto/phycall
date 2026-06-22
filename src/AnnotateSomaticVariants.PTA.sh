#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tamara.prieto.fernandez@gmail.com
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 1
#SBATCH -t 05:00:00
#SBATCH --mem 50G

source ReadConfig.sh $1
PATIENT=$(basename $1 | sed 's/Config.//' | sed 's/.txt//')

cd /gpfs/commons/groups/landau_lab/tprieto/apps/annovar

# the output file will be named based on the prefix in out followed by "hg38_multianno.vcf"
./table_annovar.pl \
	${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny.vcf \
	humandb/ \
	-buildver hg38 \
	-out ${WORKDIR}/isec_somatic-germline/FILTERED/VariantSitesForPhylogeny \
	-protocol refGene,dbnsfp42c,cosmic70,avsnp150,exac03,clinvar_20220320 \
	-remove \
	-operation g,f,f,f,f,f \
	-nastring . --vcfinput

# The -operation argument tells ANNOVAR which operations to use for each of the protocols: g means gene-based, gx means gene-based with cross-reference annotation (from -xref argument), r means region-based and f means filter-based.
