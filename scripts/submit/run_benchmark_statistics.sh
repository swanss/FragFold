#!/bin/bash
##SBATCH -N 1
##SBATCH --job-name benchmark_statistics
##SBATCH --tasks-per-node=1
##SBATCH --time=1-00:00:00
##SBATCH --mem=8G
##SBATCH -o benchmark_statistics.%j.log

CONDA_ENV=fragfold2
REPO=/home/gridsan/sswanson/swans/savinovCollaboration/FragFold

echo $LLSUB_RANK

pred_peaks_csv=/data1/groups/keatinglab/swans/savinovCollaboration/inhibitory_fragments_structure_prediction/derived_data/peak_prediction/predict_peaks_paramscan/predictalphafoldpeaks_paramscan_batch${LLSUB_RANK}.csv
exp_peaks_csv=/data1/groups/keatinglab/swans/savinovCollaboration/analysis/param_scan_analysis/230714_experimentalpeaks_Zscorecutoff2.5_modeabsolute_grouping25_maxgapdistance5.csv
exp_peaks_known_csv=$REPO/input/data/PPI_inhibitory_fragment_peaks_AS_230114.csv

SECONDS=0

conda run -n $CONDA_ENV --no-capture-output python $REPO/fragfold/calculate_benchmark_statistics.py \
    --pred_peaks_csv $pred_peaks_csv \
    --exp_peaks_csv $exp_peaks_csv \
    --exp_peaks_known_csv $exp_peaks_known_csv \
    --by_gene \
    --min_cluster_size 6 \
    --cluster_peaks_frac_overlap 0.7 \
    --store_intermediate


ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED