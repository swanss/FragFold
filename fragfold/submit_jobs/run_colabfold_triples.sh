#!/bin/bash

module load cuda/11.1
export PATH="/home/gridsan/sswanson/conda_envs/colabfold_batch/bin:$PATH"

input=$1
result_dir=output
model_type=AlphaFold2-ptm #alt: AlphaFold2-multimer-v1 AlphaFold2-multimer-v2
pair_mode=unpaired
data=/data1/groups/keatinglab/alphafold_ss/dataset

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

SECONDS=0

echo colabfold_batch $input $result_dir --data $data --model-type $model_type --pair-mode $pair_mode
colabfold_batch $input $result_dir --data $data --model-type $model_type --pair-mode $pair_mode

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED