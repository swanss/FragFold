import math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import time

from Bio.PDB import *

import glob, json, argparse, os, sys, re, functools, itertools, random
import string
import multiprocessing as mp

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from fragfold.src.analyze_predictions import *
from fragfold.src.peak_prediction import (clusterPeaksByOverlap,
                                            calculateOverlapBetweenPredandExp,
                                            calculateBenchmarkStatistics)
from fragfold.calculate_benchmark_statistics_paramscan import (benchmark_genes,
                                            benchmark_conditions,
                                            filterLengthsByGene)

def main(args):
    ### EXP peaks
    # Load experimental data and filter by benchmark genes
    path = args.exp_peaks_csv
    exp_df = pd.read_csv(path,index_col=0)
    # exp_df['gene'] = exp_df['gene'].str.replace('rpIL-coding-EcoliBL21DE3','rplL-coding-EcoliBL21DE3') # fix misnamed gene
    # filt_exp_df = filterLengthsByGene(exp_df[(exp_df['gene'].isin(benchmark_genes))]).copy(deep=True)
    filt_exp_df = exp_df[(exp_df['gene'].isin(benchmark_genes))]
    filt_exp_df = filterLengthsByGene(filt_exp_df).copy(deep=True)
    if (args.store_intermediate):
        filt_exp_df.to_csv("filtered_experimental_peaks.csv")
    print(f"There are {len(exp_df)} peaks in the original experimental data and {len(filt_exp_df)} after filtering for genes in the benchmark")

    # Label peaks with experimental structures as "known"
    path = args.exp_peaks_known_csv
    andrew_df = pd.read_csv(path)
    print(f"There are {len(andrew_df)} peaks in the curated set")

    # Go through the experimental peaks df and check if each peak is similar to one already in the curated set
    curated_found = set()
    peak_type = []
    for i,row in filt_exp_df.iterrows():
        gene = row['gene'].split('-')[0]
        peak_start = row['peak region first fragment center (aa)']
        peak_end = row['peak region last fragment center (aa)']
        filt_df = andrew_df[(andrew_df['protein-coding gene']==gene)&
                            (andrew_df['protein-protein interaction inhibitory peak center (aa)']>=peak_start)&
                            (andrew_df['protein-protein interaction inhibitory peak center (aa)']<=peak_end)]
        if len(filt_df) > 0:
            print("Fragment is in curated set")
            for x in filt_df.index:
                curated_found.add(int(x))
            peak_type.append('known')
        else:
            peak_type.append('novel')
    filt_exp_df['peak type'] = peak_type
    print(f"{filt_exp_df.groupby('peak type').size()} of the filtered experimental peaks are known")

    ### PRED peaks
    # Load predicted peaks and filter by benchmark genes and conditions
    path = args.pred_peaks_csv
    pred_df = pd.read_csv(path,index_col=0)
    filt_pred_df = filterLengthsByGene(pred_df[(pred_df['gene'].isin(benchmark_genes))&(pred_df['condition'].isin(benchmark_conditions))]).copy(deep=True)

    # Verify that we are considering a single set of peak prediction parameters
    params_list = ['n_contacts_cutoff','n_weighted_contacts_cutoff','iptm_cutoff','contact_distance_cutoff']
    assert filt_pred_df.groupby(params_list).ngroups == 1

    pred_df_list = []
    overlap_df_list = []
    stat_df_list = []

    for min_cluster_size in range(1,50):
        print(f"min_cluster_size = {min_cluster_size}")
        # Filter by cluster size and merge highly overlapping clusters
        print(f"Dataframe has {len(filt_pred_df)} peaks. Filtering out peaks with size < {min_cluster_size}")
        clus_filt_pred_df = filt_pred_df[filt_pred_df['cluster_n_fragments'] >= min_cluster_size].copy(deep=True)
        print(f"After filtering, {len(clus_filt_pred_df)} peaks remain")
        if clus_filt_pred_df.shape[0] == 0:
            continue

        if args.cluster_peaks_frac_overlap > 0.0 and args.cluster_peaks_frac_overlap < 1.0:
            print("Merging overlapping peaks...")
            clus_filt_pred_df = clusterPeaksByOverlap(clus_filt_pred_df,frac_overlap=args.cluster_peaks_frac_overlap,verbose=False)
            print(f"After merging, {len(clus_filt_pred_df)} peaks remain...")
        clus_filt_pred_df['cluster_n_fragments_threshold'] = min_cluster_size

        # Find overlap
        assert clus_filt_pred_df.groupby(params_list).ngroups == 1 # just in case something changed in the DF...
        res_overlap_frac_cutoffs = np.arange(15,31)/30 # convert to fractions to apply to different fragment lengths
        overlap_df = calculateOverlapBetweenPredandExp(clus_filt_pred_df,filt_exp_df,res_overlap_frac_cutoffs)
        overlap_df['cluster_n_fragments_threshold'] = min_cluster_size
        stat_df = calculateBenchmarkStatistics(overlap_df,filt_exp_df,30,(0,30),args.by_gene)
        stat_df['cluster_n_fragments_threshold'] = min_cluster_size

        # print("Check cluster_n_fragments_threshold in dfs")
        # for df in pred_df_list:
        #     if df.shape[0] > 0:
        #         print(df.iloc[0]['cluster_n_fragments_threshold'])

        pred_df_list.append(clus_filt_pred_df)
        overlap_df_list.append(overlap_df)
        stat_df_list.append(stat_df)

    filt_pred_df = pd.concat(pred_df_list,ignore_index=True)
    overlap_df = pd.concat(overlap_df_list,ignore_index=True)
    stat_df = pd.concat(stat_df_list,ignore_index=True)

    if (args.store_intermediate):
        filt_pred_df.to_csv("predicted_peaks_filteredbyclustersize_mergeoverlappingpeaks.csv")
        overlap_df.to_csv(f"overlap.csv")
    stat_df.to_csv(f"benchmark_statistics.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'calculateBenchmarkStatistics',
        description = 'Script for calculating benchmark statistics given clustering results for a single set of parameters with randomized predicted peak positions')
    parser.add_argument('--pred_peaks_csv',type=str,required=True)
    parser.add_argument('--exp_peaks_csv',type=str,required=True)
    parser.add_argument('--exp_peaks_known_csv',type=str,required=True)
    parser.add_argument('--store_intermediate',action='store_true')
    parser.add_argument('--by_gene',action='store_true')
    # parser.add_argument('--min_cluster_size_',type=int,default=6,
    #                     help='Predicted peaks with a smaller size than this will be filtered out prior to peak clustering (column == `cluster_n_fragments` in the table)')
    parser.add_argument('--cluster_peaks_frac_overlap',type=float,default=-1.0,
                        help='If provided (i.e. value is in the range (0,1)), will cluster predicted peaks before calculating overlap statistics')

    args = parser.parse_args()
    print(args)
    main(args)
    print('Done!')