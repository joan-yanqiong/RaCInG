#!/usr/bin/env bash
work_dir=/Users/joankant/Library/CloudStorage/OneDrive-UHN/Coding/RaCInG
input_file="${work_dir}/Data/GBM/GBM_rnaseqv2_counts.txt"
output_dir="${work_dir}/output"
cancer_type="GBM"


echo "Computing TPM"
Rscript "${work_dir}/Rscripts/10_preprocessing.R" \
    --input_file "${input_file}" \
    --output_dir "${output_dir}/10_preprocessing"

echo "Computing cell type abundances"
input_file="$output_dir/10_preprocessing/${cancer_type}_tpm.rds"

Rscript "${work_dir}/Rscripts/20_deconv.R" \
    --input_file "${input_file}" \
    --output_dir "${output_dir}/20_deconv" \
    --cancer_type "${cancer_type}"

echo "Post-processing deconvolution results"
input_dir="$output_dir/20_deconv"
Rscript "${work_dir}/Rscripts/21_deconv_post.R" \
    --input_dir "${input_dir}" \
    --output_dir "${output_dir}/21_deconv_post" \
    --cancer_type "${cancer_type}"