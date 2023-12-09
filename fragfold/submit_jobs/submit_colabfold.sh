#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name colabfold
#SBATCH --gres=gpu:volta:1
#SBATCH --constraint=xeon-g6
#SBATCH --tasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH -o colabfold.%j.log

source $(conda info --base)/etc/profile.d/conda.sh
conda activate colabfold_batch

#GPU environment variables that frankly I don't understand
MAXRAM=$(echo `ulimit -m` '/ 1024.0'|bc)
GPUMEM=`nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits|tail -1`
export XLA_PYTHON_CLIENT_MEM_FRACTION=`echo "scale=3;$MAXRAM / $GPUMEM"|bc`
export TF_FORCE_UNIFIED_MEMORY='1'

input=$1 #may also be fasta, or dir with fastas/a3ms
result_dir=output
model_type=AlphaFold2-ptm #alt: AlphaFold2-multimer-v1 AlphaFold2-multimer-v2
pair_mode=unpaired
data=/data1/groups/keatinglab/alphafold_ss/dataset

if [ ! -d $result_dir ]
then
	mkdir $result_dir
fi

SECONDS=0

echo colabfold_batch $input $result_dir --data $data --model-type $model_type --pair-mode $pair_mode
colabfold_batch $input $result_dir --data $data --model-type $model_type --pair-mode $pair_mode

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
