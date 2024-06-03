#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name contact_recovery_analysis
#SBATCH --tasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
#SBATCH -o contact_recovery_analysis.%j.log

ENV=fragfold3

REPO=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold

SECONDS=0

conda run -n $ENV --no-capture-output python $REPO/fragfold/contact_recovery_analysis.py \
    --import_json $REPO/input/json/contact_recovery/colabfold_contact_recovery_240521fix.json \
    --colabfold_data_csv $REPO/results/v1/nextflow_predictpeaks/main/colabfold115_v1_results_expmerge.csv \
    --native_pdbs $REPO/input/protein_structures/ \
    --contact_distance_cutoff 4.0

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED