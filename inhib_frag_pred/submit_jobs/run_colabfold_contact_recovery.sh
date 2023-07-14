#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name contact_recovery_analysis
#SBATCH --tasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=8G
#SBATCH -o contact_recovery_analysis.%j.log

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate alphafold_fragment_prediction

REPO=/data1/groups/keatinglab/swans/savinovCollaboration/inhibitory_fragments_structure_prediction

SECONDS=0

python $REPO/contact_recovery_analysis.py \
    --import_json $REPO/json/colabfold_contact_recovery_example.json \
    --colabfold_data_csv $REPO/analysis/process_colabfold_output_115/colabfold_predictions.csv \
    --native_pdbs $REPO/data/protein_structures/ \
    --contact_distance_cutoff 4.0

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED