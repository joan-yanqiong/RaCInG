#!/usr/bin/env zsh
work_dir=/Users/joankant/Library/CloudStorage/OneDrive-UHN/Coding/RaCInG

input_dir="$work_dir/output/SKCM"
patient_id=0
n_patients=30
n_cells=1000
av=5
n_graphs=10
cancer_type=SKCM
weight_type=min
communication_type="W"
output_dir=$work_dir/output/SKCM/networks

source "/Users/joankant/miniforge-pypy3/bin/activate" "racing"

# for patient_id in {0..1}
# do
#     norm_dist=0
#     python3 "$work_dir/Python/1_simulate.py" \
#         --cancer_type $cancer_type \
#         --input_dir $input_dir \
#         --patient_id $patient_id \
#         --n_cells $n_cells \
#         --av $av \
#         --n_graphs $n_graphs \
#         --weight_type $weight_type \
#         --communication_type $communication_type \
#         --norm_dist $norm_dist \
#         --output_dir $output_dir/original

#     norm_dist=1
#     python3 $work_dir/Python/1_simulate.py \
#         --cancer_type $cancer_type \
#         --input_dir $input_dir \
#         --patient_id $patient_id \
#         --n_cells $n_cells \
#         --av $av \
#         --n_graphs $n_graphs \
#         --weight_type $weight_type \
#         --communication_type $communication_type \
#         --norm_dist $norm_dist \
#         --output_dir $output_dir/normalized

# done

norm_dist=0
python3 $work_dir/Python/2_combine_samples.py \
        --cancer_type $cancer_type \
        --input_dir $output_dir/original \
        --av $av \
        --communication_type $communication_type \
        --norm_dist $norm_dist \
        --output_dir $output_dir


norm_dist=1
python3 $work_dir/Python/2_combine_samples.py \
        --cancer_type $cancer_type \
        --input_dir $output_dir/normalized \
        --av $av \
        --communication_type $communication_type \
        --norm_dist $norm_dist \
        --output_dir $output_dir


# # echo $output_dir
# output_dir=/Users/joankant/Library/CloudStorage/OneDrive-UHN/Coding/RaCInG/output/SKCM/networks_v0

python3 $work_dir/Python/3_normalize.py \
    --cancer_type $cancer_type \
    --input_dir $input_dir \
    --n_cells $n_cells \
    --av $av \
    --n_graphs $n_graphs \
    --weight_type $weight_type \
    --communication_type $communication_type \
    --output_dir $output_dir