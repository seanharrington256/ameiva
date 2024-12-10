#!/bin/bash

#SBATCH --job-name iqtree
#SBATCH -A inbreh
#SBATCH -t 0-24:00
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=32
#SBATCH --mem=36G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/inbreh/ameiva/scripts/err_out/err_iqtree%A.err
#SBATCH -o /project/inbreh/ameiva/scripts/err_out/std_iqtree%A.out

# Load module necessary
module load iq-tree/2.3.6

# set up the input file, outfile prefix, and output directory
INFILE="/project/inbreh/ameiva/ipyrad_out/ameiva_dn_c92_nolowcov_outfiles/ameiva_dn_c92_nolowcov.phy"
OUTFIX="ameiva_dn_c92_nolowcov"
OUTDIR="/project/inbreh/ameiva/iqtree_out"

# change into the output directory
cd $OUTDIR

# Run iqtree
iqtree2 -s $INFILE -m MFP -T auto --prefix $OUTFIX -B 1000 -alrt 1000 -nt AUTO -ntmax 32

