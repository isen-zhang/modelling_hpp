# modelling_hpp

## Installing PLINK 1.9
PLINK 1.9 is a powerful and efficient tool for whole-genome association studies. This guide provides step-by-step instructions for downloading and installing PLINK 1.9 on a Linux system.

---

### **1. Download PLINK 1.9**

Download [Plink 1.9](https://www.cog-genomics.org/plink/)
```bash
unzip plink_linux_*.zip
cp plink_install/plink ~/bin/
chomd 755 ~/bin/plink
export PATH=~/bin:$PATH
source ~/.bashrc
plink --version
```


### **Running the Script**
An example scripts can be found in [gwas_v1.sh](https://github.com/isen-zhang/modelling_hpp/blob/main/gwas_v1.sh) , which streamlines genotype data quality control, GWAS analysis, and significant SNP extraction.



```bash
   bash run_gwas.sh
```

---

## **Input Files**

### 1. **`merge_list.txt`**
   - A text file listing the PLINK binary files (`.bed`, `.bim`, `.fam`) to be merged.
     - `.bed`: Binary genotype file.
     - `.bim`: SNP information file.
     - `.fam`: Sample information file.

   - Format: Each line should specify the path to the binary files without the file extensions. Example:
     ```
     dataset1
     dataset2
     dataset3
     ```

### 2. **`phenotype.txt`**
   - A tab-delimited file containing phenotype data.
   - Format:
     - The first two columns must match the `FID` (Family ID) and `IID` (Individual ID) in the `.fam` file.
     - Subsequent columns represent phenotypes. Each column header should indicate the phenotype name.
   - Example:
     ```
     FID IID Height Disease
     1   1   165.2    1
     1   2   170.8    0
     2   1   160.5    1
     ```
   - **Continuous Phenotypes**: Values like height or BMI.
   - **Binary Phenotypes**: Encoded as `1` (control) and `2` (case).

---
### Key Features
1. **File Merging**:
   - Combines genotype files based on a merge_list.txt.
     
2. **Quality Control**:
   - Filters genotype data based on missingness (`GENO`), minor allele frequency (`MAF`), and Hardy-Weinberg equilibrium (`HWE`).
   - GENO: Missingness rate (default: 0.05).
   - MAF: Minor allele frequency threshold (default: 0.05).
   - HWE: Hardy-Weinberg equilibrium threshold (default: 1e-6).
3. **GWAS Analysis**:
   - For each phenotype listed in the phenotype.txt file:
   - If the phenotype is continuous, linear regression is performed.
   - If the phenotype is binary, the pipeline ensures values are encoded as 1 (control) and 2 (case), converting from 0/1 if necessary. Logistic regression is then performed.

4. **Significant SNP Filtering**:
   - Extracts SNPs with a p-value < 5e-8 and saved in the significant_snps directory.
     
### Output Files
   - Quality-Controlled Genotype Data: qc_results/qc_data.*
   - GWAS Results: gwas_results/gwas_<phenotype_name>.*
   - Significant SNPs: significant_snps/<phenotype_name>_significant_snps.txt
