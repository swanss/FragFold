#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name process_colabfold_output
#SBATCH --exclusive
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
#SBATCH -o process_colabfold_output.%j.log

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate alphafold_fragment_prediction

REPO=/data1/groups/keatinglab/swans/savinovCollaboration/inhibitory_fragments_structure_prediction

SECONDS=0

python $REPO/inhib_frag_pred/colabfold_process_output.py \
    --import_json $REPO/input/json/colabfold115_output.json \
    --experimental_data $REPO/input/data/Savinov_2022_inhib_peptide_mapping.csv

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED