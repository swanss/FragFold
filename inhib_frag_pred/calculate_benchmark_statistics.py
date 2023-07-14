import glob, json
from dataclasses import dataclass

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

from src.analyze_predictions import *
from src.peak_prediction import rangeOverlap,getResidueOverlapReq,calculateOverlapBetweenPredandExp,calculateBenchmarkStatistics

def main(args):
    ### EXP peaks
    # path = '/data1/groups/keatinglab/swans/savinovCollaboration/analysis/230223_expPeakFinding_5/230223_experimentalpeaks_Zscorecutoff2.5_modeabsolute_grouping25_maxgapdistance5.csv'
    path = args.exp_peaks_csv
    exp_df = pd.read_csv(path,index_col=0)
    # exp_df['gene'] = exp_df['gene'].str.replace('rpIL-coding-EcoliBL21DE3','rplL-coding-EcoliBL21DE3') # fix misnamed gene
    print(f"There are {len(exp_df)} peaks in the experimental data")

    # path = '/data1/groups/keatinglab/swans/savinovCollaboration/inhibitory_fragments_structure_prediction/data/PPI_inhibitory_fragment_peaks_AS_230114.csv'
    path = args.exp_peaks_known_csv
    andrew_df = pd.read_csv(path)
    print(f"There are {len(andrew_df)} peaks in the curated set")

    # Go through the experimental peaks df and check if each peak is similar to one already in the curated set
    curated_found = set()
    peak_type = []
    for i,row in exp_df.iterrows():
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
    exp_df['peak type'] = peak_type
    print(f"{exp_df.groupby('peak type').size()} of the experimental peaks are known")

    ### PRED peaks
    path = args.pred_peaks_csv
    pred_df = pd.read_csv(path,index_col=0)

    # Find overlap
    params_list = ['n_contacts_cutoff','n_weighted_contacts_cutoff','iptm_cutoff','contact_distance_cutoff']
    grouped_pred_df = pred_df.groupby(params_list)
    print(grouped_pred_df.ngroups)
    res_overlap_frac_cutoffs = np.arange(15,31)/30 # convert to fractions to apply to different fragment lengths
    df_list = []
    for (n_contacts,nweighted_contacts,iptm,contact_cutoff),group_df in grouped_pred_df:
        print(n_contacts,nweighted_contacts,iptm,contact_cutoff)
        overlap_df = calculateOverlapBetweenPredandExp(group_df,exp_df,res_overlap_frac_cutoffs)
        stat_df = calculateBenchmarkStatistics(overlap_df,30)
        stat_df['n_contacts_cutoff'] = n_contacts
        stat_df['n_weighted_contacts_cutoff'] = nweighted_contacts
        stat_df['iptm_cutoff'] = iptm
        stat_df['contact_distance_cutoff'] = contact_cutoff
        df_list.append(stat_df)
    comb_df = pd.concat(df_list,ignore_index=True)

    comb_df.to_csv(f"benchmark_statistics_batch{args.batch_id}.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'calculateBenchmarkStatistics',
        description = 'Script for calculating benchmark statistics given clustering results')
    parser.add_argument('--batch_id',type=int)
    parser.add_argument('--pred_peaks_csv',type=str,required=True)
    parser.add_argument('--exp_peaks_csv',type=str,required=True)
    parser.add_argument('--exp_peaks_known_csv',type=str,required=True)

    args = parser.parse_args()
    print(args)
    main(args)
    print('Done!')