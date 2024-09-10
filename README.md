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


