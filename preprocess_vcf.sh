#!/bin/bash

# Define quality control thresholds
GENO=0.05
MAF=0.05
HWE=1e-6

# Define directories
input_dir="~/studies/human_genetics_bulk/reads/vcf"
temp_dir="../results/temp_file"
mkdir -p "$temp_dir"

# Define the starting file
start_file="gencove_vcf__5c89d2d0-418f-4177-8ee5-76020455da48.vcf.gz"
start_processing=false

# Loop through each .vcf.gz file in the specified directory
for vcf_file in "$input_dir"/*.vcf.gz; do
    # Check if we have reached the starting file
    if [[ "$(basename "$vcf_file")" == "$start_file" ]]; then
        start_processing=true
    fi

    # Skip files until we reach the specified starting file
    if ! $start_processing; then
        echo "Skipping $vcf_file"
        continue
    fi

    # Extract the base name (without .vcf.gz extension) for output naming
    base_name=$(basename "$vcf_file" .vcf.gz)

    # Convert VCF to PLINK format (.bed/.bim/.fam), saving in temp_dir
    echo "Processing $vcf_file..."
    plink --vcf "$vcf_file" --make-bed --out "$temp_dir/${base_name}_raw" --noweb

    # Perform quality control, saving filtered output in temp_dir
    plink --bfile "$temp_dir/${base_name}_raw" --geno $GENO --maf $MAF --hwe $HWE --make-bed --out "$temp_dir/${base_name}_filtered" --noweb
done

# Merge all the filtered .bed files into a single dataset
# Start by taking the first file as the base
first_file=$(ls "$temp_dir"/*_filtered.bed | head -n 1)
first_base=$(basename "$first_file" _filtered.bed)

# Create a list of the files to merge
merge_list="$temp_dir/merge_list.txt"
rm -f $merge_list
for bed_file in "$temp_dir"/*_filtered.bed; do
    if [[ "$bed_file" != "$temp_dir/${first_base}_filtered.bed" ]]; then
        base=$(basename "$bed_file" _filtered.bed)
        echo "$temp_dir/${base}_filtered" >> $merge_list
    fi
done

# Perform the merge, saving output in temp_dir
plink --bfile "$temp_dir/${first_base}_filtered" --merge-list $merge_list --make-bed --out "$temp_dir/merged_data" --noweb

echo "All VCF files processed and merged into $temp_dir/merged_data.bed, $temp_dir/merged_data.bim, and $temp_dir/merged_data.fam"
