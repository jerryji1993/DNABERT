#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH -job-name="PromID"
#SBATCH -p buyin
#SBATCH -A b1017
#SBATCH -t 168:00:00
#SBATCH --output=/projects/b1017/Jerry/DNABERT/outputs/PromID.out


module load python/anaconda3.6
source activate deeplearning
cd /projects/b1017/Jerry/PromID/promid
python3 __main__.py -I /projects/b1017/Jerry/DNABERT/data/genome/GRCh38.chr1to22.fa -O /projects/b1017/Jerry/DNABERT/results/benchmark/promid/GRCh38.promoters.promid.bed
