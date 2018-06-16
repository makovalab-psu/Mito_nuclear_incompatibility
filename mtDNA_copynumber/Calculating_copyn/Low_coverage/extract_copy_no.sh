#!/bin/bash

#SBATCH -C new
##SBATCH --cpus-per-task=8
#SBATCH -J DNAdamage
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=saz5078@psu.edu
#SBATCH -e slurm-%j.err
#SBATCH -D /nfs/brubeck.bx.psu.edu/scratch3/arslan/mtnuc_project/copy_number


source /nfs/brubeck.bx.psu.edu/scratch4/arslan/.bash_profile
source activate cenv

location=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data
iid=${1}
bam=${2}

alignment=${location}/${iid}/alignment/${bam}
index=${location}/${iid}/alignment/${bam}.bai

read_length=`samtools view ${alignment} |head -n 1000 |awk '{sum+=length($10)} END{print sum/NR}'`
echo $read_length
echo "running idxstats"

cd /nfs/brubeck.bx.psu.edu/scratch3/arslan/mtnuc_project/copy_number

samtools idxstats ${alignment} > ${iid}.idxstats

echo "calculating depth of sequencing"

Rscript cal_mtdna_copy.R ${iid}.idxstats ${read_length}



