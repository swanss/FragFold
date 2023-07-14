#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name predict_peaks
#SBATCH --tasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
#SBATCH -o predict_alphafold_peaks.%j.log

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate alphafold_fragment_prediction

REPO=/data1/groups/keatinglab/swans/savinovCollaboration/inhibitory_fragments_structure_prediction

SECONDS=0

python $REPO/predict_alphafold_peaks.py \
    --colabfold_data_csv $1 \
    --n_contacts 3 \
    --n_weighted_contacts 3 \
    --iptm 0.65 \
    --contact_distance 0.5 \
    --verbose

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED