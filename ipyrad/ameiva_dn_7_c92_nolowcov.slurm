#!/bin/bash

#SBATCH --job-name ameiva_denovo
#SBATCH -A inbreh
#SBATCH -t 7-00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/inbreh/ameiva/scripts/err_out/ameiva_%A_%a.err
#SBATCH -o /project/inbreh/ameiva/scripts/err_out/ameiva_%A_%a.out

module load miniconda3/24.3.0

## run ipyrad
conda activate ipyrad
ipyrad -p params-ameiva_dn_c92_nolowcov.txt -s 7 -c 24
