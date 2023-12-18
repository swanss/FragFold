#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name predict_peaks
#SBATCH --tasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
#SBATCH -o predict_alphafold_peaks.%j.log

CONDA_ENV=fragfold
REPO=/home/gridsan/sswanson/swans/savinovCollaboration/FragFold

SECONDS=0

conda run -n $CONDA_ENV --no-capture-output python $REPO/fragfold/predict_alphafold_peaks.py \
    --colabfold_data_csv ../processoutput/colabfold_predictions.csv \
    --n_contacts 3 \
    --n_weighted_contacts 3 \
    --iptm 0.65 \
    --contact_distance 0.5 \
    --verbose

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED