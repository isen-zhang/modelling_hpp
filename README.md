# modelling_hpp

## Installing PLINK 1.9
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

This [preprocess_vcf.sh](https://github.com/isen-zhang/modelling_hpp/blob/main/preprocess_vcf.sh) processes multiple `.vcf.gz` files in a directory, applies quality control filters using PLINK, and merges the filtered data into a single dataset.

### Quality Control Thresholds:

GENO: Missingness rate (default: 0.05).

MAF: Minor allele frequency threshold (default: 0.05).

HWE: Hardy-Weinberg equilibrium threshold (default: 1e-6).

### **Running the Script**
Change the input_dir varible with the actual directory.
```bash
   bash preprocess_vcf.sh
```

Temporary files are stored in the ../results/temp_file directory.
The final merged dataset is saved as:

results/temp_file/merged_data.bed

results/temp_file/merged_data.bim

results/temp_file/merged_data.fam

## GWAS Pipeline Script (TODO)

An example scripts can be found in [gwas_v1.sh](https://github.com/isen-zhang/modelling_hpp/blob/main/gwas_v1.sh) automates the key steps of a Genome-Wide Association Study (GWAS) pipeline, from genotype data quality control to GWAS analysis and visualization. It leverages PLINK for data processing and R for generating plots.
Create Phenotype File as create Phenotype.txt
```bash
plink --bfile temp_file/merged_data_qc_filtered --pheno phenotype.txt --assoc --out gwas_results
```

Set up covariates (Optional) as covariates.txt
```bash
# for Binary Traits
plink --bfile temp_file/merged_data_qc_filtered --pheno phenotype.txt --covar covariates.txt --logistic --out gwas_results
# for Continuous Traits
plink --bfile temp_file/merged_data_qc_filtered --pheno phenotype.txt --covar covariates.txt --linear --out gwas_results
```

### Key Features
1. **Quality Control**:
   - Filters genotype data based on missingness (`GENO`), minor allele frequency (`MAF`), and Hardy-Weinberg equilibrium (`HWE`).
   - Removes samples with high missing data (`MIND`).

2. **Population Stratification Correction**:
   - Performs PCA to account for population stratification, saving results for downstream analysis.

3. **GWAS Analysis**:
   - Conducts association testing using logistic regression (with covariates) or basic association testing (without covariates).
   - Requires `phenotype.txt` and optionally `covariates.txt`.

4. **Post-GWAS Visualization**:
   - Generates Manhattan and QQ plots of GWAS results using R.

5. **Significant SNP Filtering**:
   - Extracts SNPs with a p-value < 5e-8 and saves them to a file.

### Output Files
- **Quality-Controlled Genotype Data**: `gwas_results/merged_data_qc_filtered.*`
- **PCA Results**: `gwas_results/pca_results.eigenvec`
- **GWAS Results**: `gwas_results/gwas_results.*`
- **Plots**:
  - Manhattan Plot: `gwas_results/manhattan_plot.png`
  - QQ Plot: `gwas_results/qq_plot.png`
- **Significant SNPs**: `gwas_results/significant_snps.txt`
