#!/bin/bash
##SBATCH -N 1
##SBATCH --job-name predict_peaks_paramscan
##SBATCH --tasks-per-node=1
##SBATCH --time=1-00:00:00
##SBATCH --mem=8G
##SBATCH -o predict_alphafold_peaks_paramscan.%j.log

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate alphafold_fragment_prediction

REPO=/data1/groups/keatinglab/swans/savinovCollaboration/inhibitory_fragments_structure_prediction

SECONDS=0

python $REPO/predict_alphafold_peaks.py \
    --n_batches $LLSUB_SIZE \
    --batch_id $LLSUB_RANK \
    --colabfold_data_csv $1 \
    --verbose

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED