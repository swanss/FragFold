#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name process_colabfold_output
#SBATCH --exclusive
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
#SBATCH -o process_colabfold_output.%j.log

### -------------------------------------------------------------------- ###
# Set the following parameters

# Environments and installation paths
repo_dir=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold
conda_env_name=fragfold

# Inputs 
import_json=colabfold_output_ftsz.json # path to JSON file with colabfold output directories
experimental_data=$repo_dir/input/data/Savinov_2022_inhib_peptide_mapping.csv

### -------------------------------------------------------------------- ###

SECONDS=0

conda run --no-capture-output -n $conda_env_name python $repo_dir/fragfold/colabfold_process_output.py \
    --import_json $import_json \
    --experimental_data $experimental_data

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
