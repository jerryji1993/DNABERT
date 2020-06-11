#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=512G
#SBATCH -job-name="prepare_data_noTATAprom"
#SBATCH -p genomics-himem
#SBATCH -A b1042
#SBATCH -t 48:00:00
#SBATCH --output=/projects/b1017/Jerry/DNABERT/outputs/prepare_data_noTATAprom.out


module load python/anaconda3.6
source activate deeplearning
cd /projects/b1017/Jerry/DNABERT/scripts
python3 ./prepare_data_prom.py