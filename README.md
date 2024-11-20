# modelling_hpp

## Installing and Using PLINK 1.9
PLINK 1.9 is a powerful and efficient tool for whole-genome association studies. This guide provides step-by-step instructions for downloading and installing PLINK 1.9 on a Linux system.

---

### **1. Download PLINK 1.9**

Download [Plink 1.9](https://www.cog-genomics.org/plink/)
```bash
unzip plink_linux_*.zip
mv plink ~/bin
chomd 755 ~/bin/plink
export PATH=~/bin:$PATH
source ~/.bashrc
plink --version
```

## Run the preprocessing Script

This script processes multiple `.vcf.gz` files in a directory, applies quality control filters using PLINK, and merges the filtered data into a single dataset.
### Quality Control Thresholds:

GENO: Missingness rate (default: 0.05).
MAF: Minor allele frequency threshold (default: 0.05).
HWE: Hardy-Weinberg equilibrium threshold (default: 1e-6).

### **Running the Script**
```bash
   bash preprocessing_vcf.sh
```

Temporary files are stored in the ../results/temp_file directory.
The final merged dataset is saved as:
../results/temp_file/merged_data.bed
../results/temp_file/merged_data.bim
../results/temp_file/merged_data.fam

## GWAS Pipeline Script (TODO)

This script automates the key steps of a Genome-Wide Association Study (GWAS) pipeline, from genotype data quality control to GWAS analysis and visualization. It leverages PLINK for data processing and R for generating plots.
