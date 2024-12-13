## Ameiva RAD data assembly and analysis

Sean Harrington, Nicole Angeli, Andrew Gottscho

<br>
<br>

This repository contains code to assemble and analyze RADseq data for Carribean *Ameiva*.

<br>
<br>



### ipyrad data assembly

Directory `ipyrad`

1. `params-ameiva_denovo_c92.txt` contains ipyrad parameters for initial de novo assembly at 0.92 clustering threshold.

2. `ameiva_denovo_1_7_c92.slurm` uses the params file from `1` to execute ipyrad steps 1-7 (all steps).

3. Branch the assembly to remove individuals with low coverage, do this interactively then edit the line of the params file to allow less missing data (require locus present in 51 of 73 individuals)

```
# Load miniconda and the ipyrad environment
module load miniconda3/24.3.0
conda activate ipyrad

# branch the assembly
ipyrad -p params-ameiva_denovo_c92.txt -b ameiva_dn_c92_nolowcov ../../metadata/ameiva_names_rm_low_cov.txt
```

4. `ameiva_dn_7_c92_nolowcov.slurm` runs step 7 on the new branch with fewer individuals.


5. Branch the assembly to remove Cuba outgroup individual, do this interactively

```
# Load miniconda and the ipyrad environment
module load miniconda3/24.3.0
conda activate ipyrad

# branch the assembly
ipyrad -p params-ameiva_dn_c92_nolowcov.txt -b ameiva_dn_c92_no_outgroup ../../metadata/ameiva_names_no_out.txt
```

- `ameiva_dn_7_c92_no_outgroup.slurm` runs step 7 on the new branch with no Cuba outgroup.

4. Branch the assembly to get only the *A. exsul* clade


```
# Load miniconda and the ipyrad environment
module load miniconda3/24.3.0
conda activate ipyrad

# branch the assembly
ipyrad -p params-ameiva_dn_c92_no_outgroup.txt -b ameiva_dn_c92_exsul ../../metadata/ameiva_names_exsul.txt
```

- `ameiva_dn_7_c92_exsul.slurm` runs step 7 on the new branch with *A. exsul* clade only.

<br>
<br>


## Phylogenetic trees

Directory `trees`

1. `iqtree_ameiva.slurm` runs IQtree on the `ameiva_dn_c92_nolowcov_outfiles` dataset.

2. `svdq_as_tips.slurm` runs SVDQuartets with all individuals as tips on the same `ameiva_dn_c92_nolowcov_outfiles` dataset.





<br>
<br>

## Admixture analysis using program Admixture

Admixture binary from [here](https://dalexander.github.io/admixture/download.html)


1. Convert files to plink format


- All samples, except outgroup `ameiva_dn_c92_no_outgroup` assembly:

```
cd /project/inbreh/ameiva/ipyrad_out/ameiva_dn_c92_no_outgroup_outfiles
module load miniconda3/24.3.0
conda activate plink

# First make a thinned vcf file with only independent SNPs
vcftools --vcf ameiva_dn_c92_no_outgroup.vcf --thin 150 --recode --out ameiva_dn_c92_no_out_thin
mv ameiva_dn_c92_no_out_thin.recode.vcf ameiva_dn_c92_no_out_thin.vcf

# then make plink file
plink --vcf ameiva_dn_c92_no_out_thin.vcf --recode12 --out ameiva_dn_c92_no_out_thin --allow-extra-chr --double-id
```



- *A. exsul* only `ameiva_dn_c92_exsul` assembly:

```
cd /project/inbreh/ameiva/ipyrad_out/ameiva_dn_c92_exsul_outfiles
module load miniconda3/24.3.0
conda activate plink

# First make a thinned vcf file with only independent SNPs
vcftools --vcf ameiva_dn_c92_exsul.vcf --thin 150 --recode --out ameiva_dn_c92_exsul_thinned
mv ameiva_dn_c92_exsul_thinned.recode.vcf ameiva_dn_c92_exsul_thinned.vcf

# then make plink file
plink --vcf ameiva_dn_c92_exsul_thinned.vcf --recode12 --out ameiva_dn_c92_exsul_thin --allow-extra-chr --double-id
```

2. Run Admixture. Just running this from an interactive session because it's fast and easy:


- Run on all samples, except outgroup `ameiva_dn_c92_no_outgroup` assembly:


```
salloc -A inbreh -t 0-03:00 --mem=8G --cpus-per-task=8

cd /project/inbreh/ameiva/admix_out/ameiva_dn_c92_no_outgroup
INFILE="/project/inbreh/ameiva/ipyrad_out/ameiva_dn_c92_no_outgroup_outfiles/ameiva_dn_c92_no_out_thin.ped"
BNAME=$(basename "$INFILE" | cut -d "." -f 1)

for K in 1 2 3 4 5 6 7 8 9 10; \
	do admixture --cv $INFILE $K -j8 | tee ${BNAME}_log${K}.out; done
	
# Take a look at cross validationL
grep -h CV *log*.out
grep -h CV *log*.out > CV_summ.txt  # write to a file
```




- Run on *A. exsul* only `ameiva_dn_c92_exsul` assembly:

```
salloc -A inbreh -t 0-03:00 --mem=8G --cpus-per-task=8

cd /project/inbreh/ameiva/admix_out/ameiva_dn_c92_exsul
INFILE="/project/inbreh/ameiva/ipyrad_out/ameiva_dn_c92_exsul_outfiles/ameiva_dn_c92_exsul_thin.ped"
BNAME=$(basename "$INFILE" | cut -d "." -f 1)

for K in 1 2 3 4 5 6 7 8 9 10; \
	do admixture --cv $INFILE $K -j8 | tee ${BNAME}_log${K}.out; done
	
# Take a look at cross validationL
grep -h CV *log*.out
grep -h CV *log*.out > CV_summ.txt  # write to a file
```



<br>
<br>


## Population structure in R

Directory `R_popstruct`

1. `R_popstruct.R` runs sNMF clustering.

2. `Plot_admixture.R` plots results from Admixture. 




