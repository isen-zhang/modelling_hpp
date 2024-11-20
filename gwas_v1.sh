#!/bin/bash

# Define quality control thresholds
GENO=0.05
MAF=0.01
HWE=1e-6

# Define directories
input_dir="temp_file"
output_dir="gwas_results"
mkdir -p "$output_dir"

# Step 1: Quality Control on Genotype Data
echo "Starting Quality Control on merged genotype data..."
plink --bfile "$input_dir/merged_data" --geno $GENO --maf $MAF --hwe $HWE --make-bed --out "$output_dir/merged_data_qc" --noweb
plink --bfile "$output_dir/merged_data_qc" --mind 0.1 --make-bed --out "$output_dir/merged_data_qc_filtered" --noweb
echo "Quality Control completed."

# Step 2: Population Stratification Correction (PCA)
echo "Performing PCA for population stratification correction..."
plink --bfile "$output_dir/merged_data_qc_filtered" --pca 10 --out "$output_dir/pca_results" --noweb
echo "PCA completed. Results saved to $output_dir/pca_results.eigenvec for downstream analysis."

# Step 3: Check Phenotype and Covariate Files
# Ensure phenotype.txt and optionally covariates.txt are available
if [[ ! -f "phenotype.txt" ]]; then
    echo "Error: phenotype.txt file not found. Please prepare the phenotype file."
    exit 1
fi

if [[ ! -f "covariates.txt" ]]; then
    echo "Warning: covariates.txt file not found. Proceeding without covariates."
fi

# Step 4: GWAS Analysis
echo "Starting GWAS analysis..."
if [[ -f "covariates.txt" ]]; then
    # Run logistic regression for binary traits with covariates
    plink --bfile "$output_dir/merged_data_qc_filtered" --pheno phenotype.txt --covar covariates.txt --logistic --out "$output_dir/gwas_results" --noweb
else
    # Run basic association test without covariates
    plink --bfile "$output_dir/merged_data_qc_filtered" --pheno phenotype.txt --assoc --out "$output_dir/gwas_results" --noweb
fi
echo "GWAS analysis completed."

# Step 5: Post-GWAS Visualization (Manhattan and QQ Plots)
echo "Generating Manhattan and QQ plots for GWAS results..."
Rscript - <<EOF
library(qqman)

# Load GWAS results
gwas_data <- read.table("$output_dir/gwas_results.assoc", header=TRUE)

# Generate Manhattan plot
png("$output_dir/manhattan_plot.png", width=1000, height=600)
manhattan(gwas_data, chr="CHR", bp="BP", p="P", snp="SNP", main="Manhattan Plot of GWAS Results")
dev.off()

# Generate QQ plot
png("$output_dir/qq_plot.png", width=600, height=600)
qq(gwas_data\$P, main="QQ Plot of GWAS Results")
dev.off()
EOF
echo "Manhattan and QQ plots saved to $output_dir."

# Step 6: Filter Significant SNPs
echo "Filtering significant SNPs with p-value < 5e-8..."
awk '$9 < 5e-8' "$output_dir/gwas_results.assoc" > "$output_dir/significant_snps.txt"
echo "List of significant SNPs saved to $output_dir/significant_snps.txt."

echo "GWAS pipeline completed."
