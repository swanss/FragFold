from Bio.PDB import *

import glob, json, argparse, os, sys, re, functools, itertools, random, math
import string
import multiprocessing as mp

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from fragfold.src.analyze_predictions import *
from fragfold.src.peak_prediction import (clusterOverlap,
                                        calculateOverlapBetweenPredandExp,
                                        calculateBenchmarkStatistics,
                                        peakResLength)
from fragfold.calculate_benchmark_statistics_paramscan import (benchmark_genes,
                                                               benchmark_conditions,
                                                               filterLengthsByGene)

benchmark_gene_fragment_start_range = {
    'ftsZ-coding-EcoliBL21DE3':(1,383),
    'groL-coding-EcoliBL21DE3':(1,548),
    'gyrA-coding-EcoliBL21DE3':(1,872),
    'lptG-coding-EcoliBL21DE3':(1,360),
    'rplL-coding-EcoliBL21DE3':(1,121),
    'ssb-coding-EcoliBL21DE3' :(1,178),
}

def randomizePredictedPeakPositions(peak_df,avg_peak_widths=False):
    copy_peak_df = peak_df.copy(deep=True)

    # drop values that will no longer have meaning after randomizing predicted peak positions
    copy_peak_df.drop(columns=['fragment start (aa)', 'fragment center (aa)', 'fragment end (aa)', 'cluster_member_idx'])

    copy_peak_df['cluster width (aa)'] = copy_peak_df['cluster last fragment center (aa)'] - copy_peak_df['cluster first fragment center (aa)']
    copy_peak_df['cluster first fragment center (aa)'] = copy_peak_df.gene.apply(lambda x: random.randint(*benchmark_gene_fragment_start_range[x]))
    if not avg_peak_widths:
        # preserve the number of predicted peaks and their width, but randomize their positions
        copy_peak_df['cluster last fragment center (aa)'] = copy_peak_df['cluster first fragment center (aa)'] + copy_peak_df['cluster width (aa)']
    else:
        # average the peak widths
        copy_peak_df['cluster width (aa)'] = np.round(copy_peak_df['cluster width (aa)'].mean())
        copy_peak_df['cluster last fragment center (aa)'] = copy_peak_df['cluster first fragment center (aa)'] + copy_peak_df['cluster width (aa)']

    return copy_peak_df

def randomizePredictedPeakPositionsNoOverlap(peak_df,
                                             group_vars=['gene','condition'],
                                             avg_peak_widths=False,
                                             maxFracOverlapToExistingPeak=0.7):
    '''
    This version of the function verifies that a peak does not overlap any existing ones by more than maxFracOverlapToExistingPeak
    '''
    copy_peak_df = peak_df.copy(deep=True)

    # drop values that will no longer have meaning after randomizing predicted peak positions
    copy_peak_df.drop(columns=['fragment start (aa)', 'fragment center (aa)', 'fragment end (aa)',
                               'cluster first residue', 'cluster last residue'])
    copy_peak_df['cluster width (aa)'] = copy_peak_df['cluster last fragment center (aa)'] - copy_peak_df['cluster first fragment center (aa)']
    if avg_peak_widths:
        copy_peak_df['cluster width (aa)'] = np.round(copy_peak_df['cluster width (aa)'].mean())
        copy_peak_df['cluster last fragment center (aa)'] = copy_peak_df['cluster first fragment center (aa)'] + copy_peak_df['cluster width (aa)']

    all_resampled_peaks_df_list = []
    for (gene,condition),group_df in copy_peak_df.groupby(group_vars):
        resampled_peaks_list = []
        for i,peak in group_df.iterrows():
            overlapsExistingPeak = True
            first_center = -1
            while overlapsExistingPeak:
                # try sampling new positions for the peak
                peak['cluster first fragment center (aa)'] = random.randint(*benchmark_gene_fragment_start_range[peak['gene']])
                peak['cluster last fragment center (aa)'] = peak['cluster first fragment center (aa)'] + peak['cluster width (aa)']

                # check if it overlaps a peak that was already placed
                overlapsExistingPeak = False
                for placed_peak in resampled_peaks_list:
                    smaller_peak,larger_peak = (peak,placed_peak) if peakResLength(peak) <= peakResLength(placed_peak) else (placed_peak,peak)
                    if clusterOverlap(smaller_peak,larger_peak) >= maxFracOverlapToExistingPeak:
                        print("Peak overlaps existing peak: resample")
                        overlapsExistingPeak = True
                        break
            
            resampled_peaks_list.append(peak)
                
        resampled_group_df = pd.DataFrame(resampled_peaks_list)
        resampled_group_df['gene'] = gene
        resampled_group_df['condition'] = condition
        all_resampled_peaks_df_list.append(resampled_group_df)

    resampled_peaks_df = pd.concat(all_resampled_peaks_df_list,ignore_index=True)

    return resampled_peaks_df

def main(args):
    ### EXP peaks
    # Load experimental data and filter by benchmark genes
    path = args.exp_peaks_csv
    exp_df = pd.read_csv(path,index_col=0)
    # exp_df['gene'] = exp_df['gene'].str.replace('rpIL-coding-EcoliBL21DE3','rplL-coding-EcoliBL21DE3') # fix misnamed gene
    # filt_exp_df = filterLengthsByGene(exp_df[(exp_df['gene'].isin(benchmark_genes))]).copy(deep=True)
    filt_exp_df = exp_df[(exp_df['gene'].isin(benchmark_genes))]
    filt_exp_df = filterLengthsByGene(filt_exp_df).copy(deep=True)
    if (args.batch_id == 0):
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
    orig_filt_pred_df = filterLengthsByGene(pred_df[(pred_df['gene'].isin(benchmark_genes))&(pred_df['condition'].isin(benchmark_conditions))]).copy(deep=True)

    # Verify that we are considering a single set of peak prediction parameters
    params_list = ['n_contacts_cutoff','n_weighted_contacts_cutoff','iptm_cutoff','contact_distance_cutoff']
    assert orig_filt_pred_df.groupby(params_list).ngroups == 1

    # For each sample in this batch, find a new position for each predicted peak
    batch_size = math.ceil(args.n_samples/args.n_batches)
    k_start = args.batch_id*batch_size
    k_stop = min(args.n_samples,(args.batch_id+1)*batch_size)
    for k in range(k_start,k_stop):
        random.seed(args.random_seed+k)

        if (args.cluster_peaks_frac_overlap > 0 and args.cluster_peaks_frac_overlap < 1):
            print(f"Randomize with a check to ensure no peaks have > {args.cluster_peaks_frac_overlap} overlap")
            filt_pred_df = randomizePredictedPeakPositionsNoOverlap(orig_filt_pred_df,
                                                                    avg_peak_widths=args.avg_peak_widths,
                                                                    maxFracOverlapToExistingPeak=args.cluster_peaks_frac_overlap)
        else:
            print(f"Randomize the position of each peak independently")
            filt_pred_df = randomizePredictedPeakPositions(orig_filt_pred_df,args.avg_peak_widths)
        if (args.store_intermediate):
            filt_pred_df.to_csv(f"filtered_predicted_peaks_sample{k}.csv")
        print(f"There are {len(pred_df)} peaks in the original predicted peaks data and {len(filt_pred_df)} after filtering for predicted peaks in the benchmark")

        # Find overlap
        grouped_pred_df = filt_pred_df.groupby(params_list)
        print(grouped_pred_df.ngroups)
        res_overlap_frac_cutoffs = np.arange(15,31)/30 # convert to fractions to apply to different fragment lengths
        overlap_df_list = []
        stat_df_list = []
        for (n_contacts,nweighted_contacts,iptm,contact_cutoff),group_df in grouped_pred_df:
            print(n_contacts,nweighted_contacts,iptm,contact_cutoff)
            overlap_df = calculateOverlapBetweenPredandExp(group_df,filt_exp_df,res_overlap_frac_cutoffs)
            stat_df = calculateBenchmarkStatistics(overlap_df,filt_exp_df,30,args.by_gene)

            overlap_df['n_contacts_cutoff'] = n_contacts
            overlap_df['n_weighted_contacts_cutoff'] = nweighted_contacts
            overlap_df['iptm_cutoff'] = iptm
            overlap_df['contact_distance_cutoff'] = contact_cutoff
            overlap_df['batch'] = args.batch_id
            overlap_df['k_sample'] = k
            overlap_df_list.append(overlap_df)

            stat_df['n_contacts_cutoff'] = n_contacts
            stat_df['n_weighted_contacts_cutoff'] = nweighted_contacts
            stat_df['iptm_cutoff'] = iptm
            stat_df['contact_distance_cutoff'] = contact_cutoff
            stat_df['batch'] = args.batch_id
            stat_df['k_sample'] = k
            stat_df_list.append(stat_df)
        comb_overlap_df = pd.concat(overlap_df_list,ignore_index=True)
        comb_stat_df = pd.concat(stat_df_list,ignore_index=True)

        if (args.store_intermediate):
            comb_overlap_df.to_csv(f"overlap_sample{k}.csv")
        comb_stat_df.to_csv(f"benchmark_statistics_sample{k}.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'calculateBenchmarkStatistics',
        description = 'Script for calculating benchmark statistics given clustering results for a single set of parameters with randomized predicted peak positions')
    parser.add_argument('--batch_id',type=int,required=True)
    parser.add_argument('--n_batches',type=int,required=True)
    parser.add_argument('--n_samples',type=int,default=1000)
    parser.add_argument('--random_seed',type=int,default=42)
    parser.add_argument('--avg_peak_widths',action='store_true')
    parser.add_argument('--pred_peaks_csv',type=str,required=True)
    parser.add_argument('--exp_peaks_csv',type=str,required=True)
    parser.add_argument('--exp_peaks_known_csv',type=str,required=True)
    parser.add_argument('--store_intermediate',action='store_true')
    parser.add_argument('--by_gene',action='store_true')
    parser.add_argument('--cluster_peaks_frac_overlap',type=float,default=-1.0,
                        help='If provided (i.e. value is in the range (0,1)), will cluster predicted peaks before calculating overlap statistics')

    args = parser.parse_args()
    print(args)
    main(args)
    print('Done!')