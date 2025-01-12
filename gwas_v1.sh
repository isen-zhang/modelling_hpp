#!/bin/bash

# Step 0: Set file paths
MERGE_LIST="merge_list.txt"  # List of PLINK binary files to merge
BFILE="merged_data"         # Output file after merging
PHENO_FILE="phenotype.txt" # Phenotype file
QC_DIR="qc_results"        # Directory to store QC results
GWAS_DIR="gwas_results"    # Directory to store GWAS results
SIGNIFICANT_DIR="significant_snps"  # Directory for significant SNPs
P_THRESHOLD=5e-8            # P-value threshold for significant SNPs

# Step 1: Merge PLINK binary files
echo "Merging PLINK files..."
plink --merge-list $MERGE_LIST --make-bed --out $BFILE
if [ $? -ne 0 ]; then
    echo "Error during merging. Check logs."
    exit 1
fi

# Step 2: Perform Quality Control (QC)
echo "Performing Quality Control (QC)..."
mkdir -p $QC_DIR
plink \
    --bfile $BFILE \
    --geno 0.05 \
    --mind 0.05 \
    --maf 0.05 \
    --hwe 1e-6 \
    --make-bed \
    --out $QC_DIR/qc_data
if [ $? -ne 0 ]; then
    echo "Error during QC. Check logs."
    exit 1
fi

# Step 3: Perform GWAS for each phenotype
mkdir -p $GWAS_DIR
PHENOTYPES=$(head -1 $PHENO_FILE | cut -f3- | tr '\t' '\n')
echo "Performing GWAS for each phenotype..."
for PHENOTYPE in $PHENOTYPES; do
    PHENOTYPE_CLEAN=$(echo "$PHENOTYPE" | sed 's/ /_/g')
    echo "Processing phenotype: $PHENOTYPE"

    # Check if phenotype is binary or continuous
    FIRST_VALUE=$(awk -v col="$PHENOTYPE" 'NR==1 {for (i=1; i<=NF; i++) if ($i == col) col_idx=i} NR==2 {print $col_idx}' $PHENO_FILE)
    SECOND_VALUE=$(awk -v col="$PHENOTYPE" 'NR==1 {for (i=1; i<=NF; i++) if ($i == col) col_idx=i} NR==3 {print $col_idx}' $PHENO_FILE)

    if [[ "$FIRST_VALUE" =~ ^[0-9]+([.][0-9]+)?$ ]] && [[ "$SECOND_VALUE" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        # Continuous phenotype: perform linear GWAS
        plink \
            --bfile $QC_DIR/qc_data \
            --pheno $PHENO_FILE \
            --pheno-name "$PHENOTYPE" \
            --linear \
            --allow-no-sex \
            --out $GWAS_DIR/gwas_$PHENOTYPE_CLEAN
    else
        # Binary phenotype: check if values are 1,2, otherwise convert
        awk -v col="$PHENOTYPE" 'BEGIN {FS=OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i == col) col_idx=i; print; next} {
            if ($col_idx == 0) $col_idx=1;
            else if ($col_idx == 1) $col_idx=2;
            print
        }' $PHENO_FILE > temp_pheno.txt

        plink \
            --bfile $QC_DIR/qc_data \
            --pheno temp_pheno.txt \
            --pheno-name "$PHENOTYPE" \
            --logistic \
            --allow-no-sex \
            --out $GWAS_DIR/gwas_$PHENOTYPE_CLEAN

        rm temp_pheno.txt
    fi

    if [ $? -ne 0 ]; then
        echo "Error during GWAS for phenotype $PHENOTYPE. Check logs."
        exit 1
    fi

done

# Step 4: Extract significant SNPs
mkdir -p $SIGNIFICANT_DIR
echo "Extracting significant SNPs..."
for FILE in $GWAS_DIR/*.assoc.*; do
    BASENAME=$(basename $FILE)
    awk -v pval="$P_THRESHOLD" '$9 < pval' $FILE > $SIGNIFICANT_DIR/${BASENAME}_significant_snps.txt
    echo "Significant SNPs saved to $SIGNIFICANT_DIR/${BASENAME}_significant_snps.txt"
done

echo "GWAS pipeline completed successfully."
