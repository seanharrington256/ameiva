#!/bin/bash

#SBATCH --job-name SVDQbootstrap
#SBATCH -A inbreh
#SBATCH -t 0-24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24   
#SBATCH --mem=36G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sharrin2@uwyo.edu
#SBATCH -e /project/inbreh/ameiva/scripts/err_out/err_svdqboot_%A_%a.err
#SBATCH -o /project/inbreh/ameiva/scripts/err_out/std_svdqboot_%A_%a.out


# Load R module to combine bootstraps with R script below
module load gcc/14.2.0 r/4.4.0

# Had to install PAUP in my software directory first and add it to my path
# variable for PAUP version installed, in path - or full path if not in path
PAUP=paup4a168_centos64

# set up some input & output file paths
OUTDIR="/project/inbreh/ameiva/svdq_out"

# Path to R script to combine bootstraps onto best tree
BOOTSCRIPT="/project/inbreh/software/bootanno.R"

#define slurm variables for the PAUPblock
filebname="ameiva_dn_c92_nolowcov_tips" #basename for all produced files
infile="/project/inbreh/ameiva/ipyrad_out/ameiva_dn_c92_nolowcov_outfiles/ameiva_dn_c92_nolowcovPAUP.nex" #name of input nexus file; can give a path so the input files don't have to be part of the working directory
nthreads=$SLURM_CPUS_PER_TASK #number of threads in the slurm script
nreps=200 #number of reps for bootstrapping


# Note that to make the ameiva_dn_c92_nolowcovPAUP.nex file, I took the regular .nex output from ipyrad and just deleted the charsets



################################################################################################################################################################
################################################################################################################################################################
####    Run based on the parameters set above
################################################################################################################################################################
################################################################################################################################################################


#change working directory to where your output files will go
cd $OUTDIR


cat <<EOF > $filebname.paup.txt
Begin paup;
set autoclose=yes warntree=no warnreset=no flock=no;
log start file=$filebname.log ;
execute $infile;
svdQuartets evalQuartets=all showScores=no ambigs=distribute bootstrap=no nthreads=$nthreads;
savetrees file=$filebname.besttree.tre;
svdQuartets evalQuartets=all showScores=no ambigs=distribute bootstrap=standard nreps=$nreps nthreads=$nthreads treefile=$filebname.svdqboots.tre;  
quit; 
end;
EOF

$PAUP $filebname.paup.txt #execute your new paup block file

Rscript $BOOTSCRIPT $filebname.besttree.tre $filebname.svdqboots.tre  #run R script with the best tree and bootstrap values 
