#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=120G
#SBATCH -job-name="additional_noTATAprom_data"
#SBATCH -p buyin
#SBATCH -A b1017
#SBATCH -t 168:00:00
#SBATCH --output=/projects/b1017/Jerry/DNABERT/outputs/additonal_noTATAprom_data.out


module load python/anaconda3.6
source activate deeplearning
cd /projects/b1017/Jerry/DNABERT/scripts
python3 ./additional_prom_data.py