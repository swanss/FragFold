import glob, json, os, shutil
from dataclasses import dataclass

import numpy as np
import pandas as pd
import seaborn as sns
import scipy
import matplotlib.pyplot as plt

from Bio.PDB import PDBParser

from src.analyze_predictions import *

# Code inspired by example from https://observablehq.com/@yurivish/peak-detection
# Also referenced http://constans.pbworks.com/w/file/fetch/120908295/Simple_Algorithms_for_Peak_Detection_in_Time-Series.pdf

# NOTE: code assumes that there are no missing positions in the experimental data (e.g. complete coverage of fragmentss)

def peakFindingFromExpData(exp_df: pd.DataFrame,
                           exp_val_name: str,
                           exp_val_threshold: float = -3.0,
                           lookaround_distance: int = 30,
                           grouping_distance: int = 25,
                           max_gap_distance: int = 5) -> pd.DataFrame:
    
    # Order fragments by their position in the gene and create "dummy fragments" if any are missing
    min_fragment_start = exp_df['fragment start (aa)'].min()
    max_fragment_start = exp_df['fragment start (aa)'].max()
    fragstart_exp_df = exp_df.set_index('fragment start (aa)')
    
    fragment_data = []
    n_missing_fragments = 0
    idx2fragmentcenter = dict()
    for i,fragment_start in enumerate(range(min_fragment_start,max_fragment_start)):
        if fragment_start in fragstart_exp_df.index:
            row = fragstart_exp_df.loc[fragment_start]
            fragment_data.append(row[exp_val_name])
            idx2fragmentcenter[i] = row['fragment center (aa)']
        else:
            # missing fragment
            fragment_data.append(np.nan)
            n_missing_fragments+=1
            
    print(f"{n_missing_fragments} missing fragments in experimental data")
    print(f"{len(fragment_data)} fragments in range {min_fragment_start}-{max_fragment_start}")
                
    peak_df,peak_score_list = peakFindingAlgorithm(data=fragment_data,peak_score_type='absolute',lookaround_distance=lookaround_distance,grouping_distance=grouping_distance,max_gap_distance=max_gap_distance,val_cutoff=exp_val_threshold)
    
    peak_df['fragment center (aa)'] = peak_df['peak_fragment_idx'].apply(lambda x: idx2fragmentcenter[x])
    peak_df['peak region first fragment center (aa)'] = peak_df['peak_region_first_idx'].apply(lambda x: idx2fragmentcenter[x])
    peak_df['peak region last fragment center (aa)'] = peak_df['peak_region_last_idx'].apply(lambda x: idx2fragmentcenter[x])

    peak_df = peak_df[['fragment center (aa)','peak region first fragment center (aa)','peak region last fragment center (aa)','peak_fragment_score']].merge(exp_df,on='fragment center (aa)')
    
    return peak_df,peak_score_list

def peakFindingAlgorithm(data: list,
                         peak_score_type: str,
                         lookaround_distance: int,
                         grouping_distance: int,
                         max_gap_distance: int,
                         val_cutoff: float):
    '''Finds peaks given experimentally measured inhibitory effect values
        
    Parameters
    ----------
    data : list
        The inhibitory effect values for each fragment of a gene. NOTE: algorithm assumes that fragments 
        are not missing, if experimental data is incomplete make sure to include missing fragments with
        null values (e.g. np.nan)
    peak_score_type : str
        The method for determining the score of each fragment. If 'absolute', then the score is 
        just the inhibitory effect value. If 'local' then the score is determined by considering
        a local context window with size determined by lookaround_distance
    lookaround_distance : int
        The number of fragments to consider in each direction when defining a peak score given a local
        context window.
    grouping_distance : int
        The maximum number of fragments to include in a single peak. For example, if 25 that will mean
        that all fragments share at least 5 residues
    max_gap_distance : int
        The maximum number of fragments that can be missing in a peak before it is considered to end.
    val_cutoff: int
        The cutoff that is used when determining whether a fragment is inhibitory
    '''
    
    print(peak_score_type,lookaround_distance,grouping_distance,max_gap_distance,val_cutoff)
    
    if peak_score_type not in {'absolute','local'}:
        raise ValueError(f"peak_score_type = {peak_score_type} is not recognized. Please use either 'absolute' or 'local'")
    
    # 1) Get "sharpness" score for each potential peak
    peak_scores = []
    peak_scores_idx = []
    for i,val in enumerate(data):
        if peak_score_type == 'local':
            left_vals = data[max(i-lookaround_distance,0):i]
            right_vals = data[i:min(i+lookaround_distance,len(data))]
            score = peakScoreLocal(val,left_vals,right_vals)
        else:
            score = val

        peak_scores.append(score)
        
        # 2) Ignore all peak centers that have a max value below the threshold
        if val < val_cutoff: 
            peak_scores_idx.append(i)
    
    # 3) Group peaks within distance (unless the gap size is too large)
    peak_groups = []
    gap_size = 0
    for i in peak_scores_idx:
        if len(peak_groups) == 0:
            peak_groups.append([i])
        else:
            peak_width = i - peak_groups[-1][0] + 1
            gap_size = i - peak_groups[-1][-1] - 1
            if (peak_width < grouping_distance) and (gap_size <= max_gap_distance):
                peak_groups[-1].append(i)
            else:
                peak_groups.append([i])
    print(f"Found {len(peak_groups)} peak groups in total")

    # 4) Select the peak representative by score
    representative_peaks = []
    representative_peak_scores = []
    start = []
    end = []
    for peak_group in peak_groups:
        peak_score = peak_scores[peak_group[0]]
        peak_idx = peak_group[0]
        print(peak_idx)
        for idx in peak_group:
            if peak_scores[idx] < peak_score:
                peak_score = peak_scores[idx]
                peak_idx = idx

        representative_peaks.append(peak_idx)
        representative_peak_scores.append(peak_score)
        start.append(peak_group[0])
        end.append(peak_group[-1])
    
    peak_df = pd.DataFrame({
        'peak_fragment_idx':representative_peaks,
        'peak_fragment_score':representative_peak_scores,
        'peak_region_first_idx':start,
        'peak_region_last_idx':end
    })
    return peak_df,peak_scores

def peakScoreLocal(val,left_vals,right_vals):
    if len(left_vals) == 0 and len(right_vals) == 0:
        raise ValueError("No neighboring values to define a peak score")
    window_vals = np.array(left_vals + right_vals)
#     print(window_vals)
    scores = val - window_vals
#     print(scores)
    return np.mean(scores)

def generatePredictionsFromAlphaFold(pred_df: pd.DataFrame,
                        pred_metric_name: str,
                        pred_metric_cutoff: float,
                        cluster_d_cutoff: float = 0.5,
                        only_output_df: bool = False) -> pd.DataFrame:

    # filter dataframe by selected metric
    filt_pred_df = pred_df[pred_df[pred_metric_name]>=pred_metric_cutoff]
    print(f"Original dataframe had {len(pred_df)} predicted fragments, while new dataframe has {len(filt_pred_df)} passing the cutoff")

    # group by condition and cluster by contacts
    parser = PDBParser(QUIET=True)
    df_list = []
    for (gene,condition),group_df in filt_pred_df.groupby(['gene','condition']):
        print(f"clustering {len(group_df)} predictions from {gene}_{condition} ...")

        # define contacts
        contacts_list = []
        for i,row in group_df.iterrows():
            s_extract = parser.get_structure("s", row['path'])
            fixResidueNumbers(s_extract[0]['B'],row['fragment start (aa)'])
            contacts_residues = getInterfaceContactsFromStructure(s_extract,{'A'},{'B'},4.0)
            # print(f"pred fragment {row['fragment_name']} has {len(contacts_residues)} contacts")
            contacts_list.append(contacts_residues)
        group_df['contact_set'] = contacts_list

        # get similarity matrix
        sim_list = []
        for frag_contacts_a in group_df['contact_set']:
            for frag_contacts_b in group_df['contact_set']:
                sim_list.append(contactOverlap(frag_contacts_a,frag_contacts_b))
        sim_matrix = np.array(sim_list).reshape(len(group_df['contact_set']),len(group_df['contact_set']))
        if not only_output_df:
            plt.figure(figsize = (40,40))
            ax = sns.heatmap(sim_matrix,
                        xticklabels=group_df['fragment start (aa)'],
                        yticklabels=group_df['fragment start (aa)'])
            ax.set_xlabel('fragment start (aa)')
            ax.set_ylabel('fragment start (aa)')
            plt.savefig(f"{gene}-{condition}_contactsim_fragments{pred_metric_cutoff}{pred_metric_name}.png",dpi=300)

        # hierarchically cluster
        dist_matrix = 1 - sim_matrix
        condensed_dist_matrix = scipy.spatial.distance.squareform(dist_matrix,force='tovector')
        linkage_matrix = scipy.cluster.hierarchy.linkage(condensed_dist_matrix, method='single', optimal_ordering=False)
        clusters = scipy.cluster.hierarchy.fcluster(linkage_matrix, t=cluster_d_cutoff, criterion='distance')
        print(f"{len(clusters)} fragments in {clusters.max()} clusters")
        print(clusters)
        group_df['cluster'] = clusters
        # Get the cluster representative for each cluster
        pred_clusrep_df = group_df.loc[group_df.groupby('cluster')[pred_metric_name].idxmax()].sort_values(by='fragment start (aa)')
        # get the first and last member for each cluster (to define the range of fragments that fall into it)
        cluster_size = []
        first_fragment_center = []
        last_fragment_center = []
        for i,row in pred_clusrep_df.iterrows():
            cluster_n = row['cluster']
            cluster_n_df = group_df[group_df['cluster']==cluster_n]
            cluster_size.append(len(cluster_n_df))
            first_fragment_center.append(cluster_n_df['fragment center (aa)'].min())
            last_fragment_center.append(cluster_n_df['fragment center (aa)'].max())
        pred_clusrep_df['cluster_size'] = cluster_size
        group_df['cluster first fragment center (aa)'] = first_fragment_center
        group_df['cluster last fragment center (aa)'] = first_fragment_center

        df_list.append(pred_clusrep_df)

        # create directory for storing cluster representatives
        if not only_output_df:
            dir_name = f"{gene}_{condition}"
            try:
                os.mkdir(dir_name)
            except FileExistsError:
                pass
            for i,row in pred_clusrep_df.iterrows():
                file = row['path'].split('/')[-1]
                shutil.copyfile(row['path'],os.path.join(dir_name,file))

        # write out cluster info to file
        if not only_output_df:
            plt.figure(figsize = (10,5))
            exp_df = pred_df[(pred_df['gene']==gene)&(pred_df['condition']==condition)]
            ax = sns.lineplot(data=exp_df,x='fragment center (aa)',y='E = inhibitory effect (enrichment)',zorder=1)
            ax.axhline(-3,c='r',ls='--')
            fcenter = pred_clusrep_df['fragment center (aa)']
            wcontacts = pred_clusrep_df['weighted_contacts']/pred_clusrep_df['weighted_contacts'].max()*exp_df['E = inhibitory effect (enrichment)'].min()
            sns.scatterplot(x=fcenter,y=wcontacts,color='black',ax=ax,zorder=2)
            plt.show()
            plt.savefig(f"{gene}-{condition}_{pred_metric_cutoff}{pred_metric_name}_inhibeffect_vs_clusrepfragments.png",dpi=300)

    return pd.concat(df_list)
