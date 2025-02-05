#!/bin/bash

# Reference files
g1000_REF="g1000_eas"
GENE_LOC="NCBI38.gene.loc"

# Variables
SUBSET_FILE="/home/ec2-user/studies/modeling-hpp/ning/data/height_testset.tsv"
GWAS_DIR="gwas_results"
GWAS_TYPE="logistic"
PHENOTYPE="height"
P_VALUE_THRESHOLD=5e-5

# Directories
WORK_DIR="magma_analysis/${PHENOTYPE}"
ANNOT_DIR="${WORK_DIR}/annotations"
RESULT_DIR="${WORK_DIR}/results"

# Create directories
mkdir -p $ANNOT_DIR $RESULT_DIR

# Process subset and run GWAS
awk '{print $1, $1}' "$SUBSET_FILE" > "${GWAS_DIR}/subset_ids.txt"
awk '{print $1, $1, $2}' "$SUBSET_FILE" > "${GWAS_DIR}/subset_phenotype.txt"

plink \
    --bfile "qc_results/qc_data" \
    --keep "${GWAS_DIR}/subset_ids.txt" \
    --pheno "${GWAS_DIR}/subset_phenotype.txt" \
    --${GWAS_TYPE} \
    --out "${GWAS_DIR}/gwas_${PHENOTYPE}"

# Determine GWAS results file
if [ -f "${GWAS_DIR}/gwas_${PHENOTYPE}.assoc.linear" ]; then
    GWAS_RESULTS="${GWAS_DIR}/gwas_${PHENOTYPE}.assoc.linear"
    echo "Using linear GWAS results: $GWAS_RESULTS"
elif [ -f "${GWAS_DIR}/gwas_${PHENOTYPE}.assoc.logistic" ]; then
    GWAS_RESULTS="${GWAS_DIR}/gwas_${PHENOTYPE}.assoc.logistic"
    echo "Using logistic GWAS results: $GWAS_RESULTS"
else
    echo "Error: No GWAS results file found for phenotype ${PHENOTYPE}"
    exit 1
fi

# Step 1: Prepare SNP Location File for Annotation
SNP_ANNOT_LOC_FILE="${WORK_DIR}/${PHENOTYPE}_snp_annot_loc.txt"
echo "Preparing SNP location file for annotation..."
awk -v pval="$P_VALUE_THRESHOLD" 'NR > 1 && $9 <= pval {print $2, $1, $3}' "$GWAS_RESULTS" > "$SNP_ANNOT_LOC_FILE"
if [ $? -ne 0 ]; then
    echo "Error: Failed to prepare SNP location file for annotation."
    exit 1
fi
echo "SNP annotation location file created: $SNP_ANNOT_LOC_FILE"

# Step 2: Prepare genes.annot File
ANNOT_FILE="${ANNOT_DIR}/${PHENOTYPE}.genes.annot"
echo "Preparing genes.annot file..."
magma --annotate \
    --snp-loc $SNP_ANNOT_LOC_FILE \
    --gene-loc $GENE_LOC \
    --out $ANNOT_FILE

if [ $? -ne 0 ]; then
    echo "Error: Failed to prepare genes.annot file."
    exit 1
fi

# Step 3: Prepare SNP Location File for P-value Analysis
SNP_PVAL_FILE="${WORK_DIR}/${PHENOTYPE}_snp_pval.txt"
echo "Preparing SNP location file for P-value analysis..."
awk -v pval="$P_VALUE_THRESHOLD" 'NR > 1 && $9 <= pval {print $2, $9}' "$GWAS_RESULTS" > "$SNP_PVAL_FILE"

# Step 4: Perform MAGMA SNP-to-Gene Analysis
RESULT_FILE="${RESULT_DIR}/${PHENOTYPE}_${P_VALUE_THRESHOLD}_gene_analysis"
echo "Running MAGMA SNP-to-Gene Analysis..."
magma --bfile $g1000_REF \
    --pval $SNP_PVAL_FILE N=10000 \
    --gene-annot ${ANNOT_FILE}.genes.annot \
    --out $RESULT_FILE
