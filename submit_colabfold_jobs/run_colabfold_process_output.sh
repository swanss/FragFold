#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name process_colabfold_output
#SBATCH --tasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
#SBATCH -o process_colabfold_output.%j.log

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate alphafold_fragment_prediction

REPO=/home/gridsan/sswanson/local_code_mirror/inhibitory_fragments_structure_prediction

python $REPO/src/colabfold_process_output.py --import_json $REPO/colabfold_output.json