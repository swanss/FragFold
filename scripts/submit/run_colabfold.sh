#!/bin/bash

echo "Running colabfold with arguments: input = ${1}, colabfold_env_dir = {$2}, alphafold_params_dir = ${3}, cuda_version = ${4}"

# Set variables
input=$1
colabfold_env_dir=$2
result_dir=output
model_type=alphafold2_ptm #alt: alphafold2_multimer_v1 alphafold2_multimer_v2 alphafold2_multimer_v3
pair_mode=unpaired
alphafold_params_dir=$3
cuda_version=$4

module load "cuda/${cuda_version}"
export PATH="${colabfold_env_dir}/bin:$PATH"

# Create output directory
if [ ! -d data ]
then
        mkdir data
fi
cd data
name=${input##*/}
dirname=${name%.a3m}
mkdir $dirname
cd $dirname

if [ ! -d $result_dir ]
then
        mkdir $result_dir
fi

# Run ColabFold
SECONDS=0
echo colabfold_batch $input $result_dir --data $alphafold_params_dir --model-type $model_type --pair-mode $pair_mode
colabfold_batch $input $result_dir --data $alphafold_params_dir --model-type $model_type --pair-mode $pair_mode
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED