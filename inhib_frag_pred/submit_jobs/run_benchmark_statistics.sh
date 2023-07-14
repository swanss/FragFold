#!/bin/bash
##SBATCH -N 1
##SBATCH --job-name benchmark_statistics
##SBATCH --tasks-per-node=1
##SBATCH --time=1-00:00:00
##SBATCH --mem=8G
##SBATCH -o benchmark_statistics.%j.log

CONDA_ROOT=/state/partition1/llgrid/pkg/anaconda/anaconda3-2022b/
source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda activate alphafold_fragment_prediction

REPO=/data1/groups/keatinglab/swans/savinovCollaboration/inhibitory_fragments_structure_prediction

echo $LLSUB_RANK

pred_peaks_csv=$REPO/analysis/predict_peaks_paramscan_3/predictalphafoldpeaks_paramscan_batch${LLSUB_RANK}.csv
exp_peaks_csv=$REPO/analysis/230223_expPeakFinding_5/230223_experimentalpeaks_Zscorecutoff2.5_modeabsolute_grouping25_maxgapdistance5.csv
exp_peaks_known_csv=$REPO/input_data/PPI_inhibitory_fragment_peaks_AS_230114.csv

SECONDS=0

python $REPO/calculate_benchmark_statistics.py \
    --batch_id $LLSUB_RANK \
    --pred_peaks_csv $pred_peaks_csv \
    --exp_peaks_csv $exp_peaks_csv \
    --exp_peaks_known_csv $exp_peaks_known_csv \

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED