#!/bin/bash
#SBATCH --mem 16G -c 24 -p short -J fetch --out fetch.%A_%a.log
CPU=24
module load sratoolkit
module load parallel-fastq-dump

N=1
if [ ! -z ${SLURM_ARRAY_TASK_ID} ]; then
 N=${SLURM_ARRAY_TASK_ID}
elif [ ! -z $1 ]; then
 N=$1
fi

SRA=sra.txt
sed -n ${N}p $SRA | while read SRARUN
do
 if [ ! -f ${SRARUN}_1.fastq.gz ]; then
	parallel-fastq-dump --tmpdir $SCRATCH --gzip  --sra-id $SRARUN --threads $CPU -O $SRARUN --split-files
 fi
done

