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

### Functions for predicting peaks using AlphaFold

def filterAlphaFoldPredictions(pred_df: pd.DataFrame,
                            n_contacts: int = 0,
                            n_wcontacts: int = 0,
                            iptm: float = 0,
                            verbose = True):

    filt_pred_df = pred_df[(pred_df['n_contacts']>=n_contacts)
                            &(pred_df['weighted_contacts']>=n_wcontacts)
                            &(pred_df['iptm']>=iptm)].copy(deep=True)
    if verbose:
        print(f"Original dataframe had {len(pred_df)} predicted fragments. After filtering by (n_contacts >= {n_contacts}) and (n_wcontacts >= {n_wcontacts}) and (iptm >= {iptm}) \n new dataframe has {len(filt_pred_df)} passing the cutoff")
    return filt_pred_df

def pairwise(iterable):
    it = iter(iterable)
    a = next(it, None)

    for b in it:
        yield (a, b)
        a = b

def splitDiscontinuousFragmentSets(pred_df: pd.DataFrame, verbose = True):
    df_list = []
    assert pred_df.gene.nunique() == pred_df.condition.nunique() == pred_df['fragment length (aa)'].nunique() == 1,f"{pred_df.gene.nunique()},{pred_df.condition.nunique(),pred_df['fragment length (aa)'].nunique()}"
    fragment_len = pred_df['fragment length (aa)'].unique()[0]

    # find boundaries between fragments (e.g. positions that share no residues) and group fragments between these boundaries
    pred_df.sort_values(by='fragment start (aa)',inplace=True)
    
    # group so that there is a single fragment per row, get gap, and merge back
    frag_df = pred_df.loc[pred_df.groupby('fragment_name',sort=False)['rank'].idxmin()].copy(True)
    frag_df['gap_to_next'] = -frag_df['fragment start (aa)'].diff(periods=-1)
    pred_df = pred_df.merge(frag_df[['fragment_name','gap_to_next']],on='fragment_name',how='left',validate='many_to_one')
    
    # get boundaries to find the contigs
    boundary_list = [0] + list(frag_df[frag_df['gap_to_next'] >= fragment_len]['fragment start (aa)']) + [frag_df['fragment start (aa)'].max()]
    pred_df['fragment contig'] = str(0)
    for fragset_n,(start,stop) in enumerate(pairwise(boundary_list)): # (start,stop]
        mask = (pred_df['fragment start (aa)']>start)&(pred_df['fragment start (aa)']<=stop)
        pred_df.loc[mask,'fragment contig'] = fragset_n
    if verbose:
        print(f"Found {pred_df.groupby('fragment contig').ngroups} non-overlapping groups of fragments")
    return pred_df

def clusterOverlappingFragments(pred_df: pd.DataFrame, 
                            cluster_dist_cutoff: float = 0.5, 
                            output_prefix: str = '',
                            verbose = True):
    # print(f"clustering {len(pred_df)} predictions from {gene}_{condition} ...")
    pred_df_copy = pred_df.copy(deep=True)

    # define contacts
    parser = PDBParser(QUIET=True)
    contacts_list = []
    for i,row in pred_df_copy.iterrows():
        s_extract = parser.get_structure("s", row['path'])
        fixResidueNumbers(s_extract[0]['B'],row['fragment start (aa)'])
        protein_chains = set(row['protein_chains'].split(','))
        fragment_chain = set(row['fragment_chain'])
        contacts_residues = getInterfaceContactsFromStructure(s_extract,protein_chains,fragment_chain,4.0)
        # print(f"pred fragment {row['fragment_name']} has {len(contacts_residues)} contacts")
        contacts_list.append(contacts_residues)
    pred_df_copy['contact_set'] = contacts_list

    if len(pred_df_copy) == 1:
        pred_df_copy['cluster'] = 1
    else:
        # get similarity matrix
        sim_list = []
        for frag_contacts_a in pred_df_copy['contact_set']:
            for frag_contacts_b in pred_df_copy['contact_set']:
                sim_list.append(contactOverlap(frag_contacts_a,frag_contacts_b))
        sim_matrix = np.array(sim_list).reshape(len(pred_df_copy['contact_set']),len(pred_df_copy['contact_set']))
        if output_prefix != '':
            plt.figure(figsize = (40,40))
            ax = sns.heatmap(sim_matrix,
                        xticklabels=pred_df_copy['fragment start (aa)'],
                        yticklabels=pred_df_copy['fragment start (aa)'])
            ax.set_xlabel('fragment start (aa)')
            ax.set_ylabel('fragment start (aa)')
            plt.savefig(f"{output_prefix}_contactsim_fragments.png",dpi=300)

        # hierarchically cluster
        dist_matrix = 1 - sim_matrix
        condensed_dist_matrix = scipy.spatial.distance.squareform(dist_matrix,force='tovector')
        linkage_matrix = scipy.cluster.hierarchy.linkage(condensed_dist_matrix, method='complete', optimal_ordering=False)
        clusters = scipy.cluster.hierarchy.fcluster(linkage_matrix, t=cluster_dist_cutoff, criterion='distance')
        if verbose:
            print(f"{len(clusters)} fragments in {clusters.max()} clusters")
            print(clusters)
        pred_df_copy['cluster'] = clusters

    # Get the cluster representative for each cluster
    pred_clusrep_df = pred_df_copy.loc[pred_df_copy.groupby('cluster')['iptm'].idxmax()].sort_values(by='fragment start (aa)')

    # get the first and last member for each cluster (to define the range of fragments that fall into it)
    cluster_size = []
    cluster_n_fragments = []
    first_fragment_center = []
    last_fragment_center = []
    cluster_member_idx_list = []
    for i,row in pred_clusrep_df.iterrows():
        cluster_n = row['cluster']
        cluster_n_df = pred_df_copy[pred_df_copy['cluster']==cluster_n]
        cluster_size.append(len(cluster_n_df))
        cluster_n_fragments.append(cluster_n_df['fragment_name'].nunique())
        first_fragment_center.append(cluster_n_df['fragment center (aa)'].min())
        last_fragment_center.append(cluster_n_df['fragment center (aa)'].max())
        cluster_member_idx_list.append(list(cluster_n_df.index))
    pred_clusrep_df['cluster_size'] = cluster_size
    pred_clusrep_df['cluster_n_fragments'] = cluster_n_fragments
    pred_clusrep_df['cluster first fragment center (aa)'] = first_fragment_center
    pred_clusrep_df['cluster last fragment center (aa)'] = last_fragment_center
    pred_clusrep_df['cluster_member_idx'] = cluster_member_idx_list
    pred_clusrep_df.drop(columns=['contact_set'],inplace=True)

    return pred_clusrep_df

# Functions for plotting clusters

def organizeClusters(cluster_df,baseline,scale,overlap_width=10):
    
    def getClusterWidth(cluster_tuple):
        return cluster_tuple[2] - cluster_tuple[1]

    def clustersOverlapping(cluster_tuple_A,cluster_tuple_B,overlap_width):
        if not ((cluster_tuple_A[2] < cluster_tuple_B[1] - overlap_width) or (cluster_tuple_B[2] < cluster_tuple_A[1] - overlap_width)):
            return True
        return False
        
    cluster_df = cluster_df.copy(deep=True)
    clusters = [(cluster,cluster_start,cluster_end) 
                 for i,(cluster,cluster_start,cluster_end) in 
                     cluster_df[['full_cluster_name','cluster first fragment center (aa)','cluster last fragment center (aa)']].iterrows()]
    clusters.sort(key=getClusterWidth,reverse=True)
    
    placed_clusters = []
    height2clusterlist = dict()
    for cluster in clusters:
        # print(cluster)
        # print(placed_clusters)
        height = 0
        overlaps = None
        if len(placed_clusters) == 0:
            placed_clusters.append((cluster,height))
            height2clusterlist[height] = [cluster]
        else:
            while overlaps is None or overlaps == True:
                overlaps = False
                
                # Try placing cluster, if it overlaps another at the same level, increment height and try again
                if height in height2clusterlist:
                    for placed_cluster in height2clusterlist[height]:
                        if clustersOverlapping(cluster,placed_cluster,overlap_width):
                            overlaps = True
                            break
                    if overlaps:
                        height += 1
                    else:
                        placed_clusters.append((cluster,height))
                        height2clusterlist[height].append(cluster)
                else:
                    overlaps = False
                    placed_clusters.append((cluster,height))
                    height2clusterlist[height] = [cluster]
                        
    print(f"{len(clusters)} before and {len(placed_clusters)} placed clusters")

    cluster_y_values = [(cluster[0],height) for cluster,height in placed_clusters]
    cluster_height_df = pd.DataFrame(cluster_y_values,columns=['full_cluster_name','height'])
    cluster_height_df['height'] = cluster_height_df['height']*scale 
    cluster_height_df['height'] = cluster_height_df['height']+baseline

    return cluster_df.merge(cluster_height_df,on='full_cluster_name',how='left')    

def plotClusters(pred_df,cluster_df,contig_df,name):
    xscale=int(10*pred_df['fragment center (aa)'].max()/300)
    fig = plt.figure(figsize=(xscale,5))
    # grouped_pred_df = pred_df.loc[pred_df.groupby('fragment_name')['n_contacts'].idxmax()]
    pred_df.groupby('fragment_name')['weighted_contacts'].mean().max()
    max_contacts_height = pred_df.groupby('fragment_name')['weighted_contacts'].mean().max()
    cluster_df = organizeClusters(cluster_df,max_contacts_height+1,.75,5)
    ax = sns.lineplot(data=pred_df,x='fragment center (aa)',y='weighted_contacts')
    plt.hlines(y=[max_contacts_height]*contig_df['fragment contig'].nunique(),
               xmin=list(contig_df.groupby('fragment contig')['fragment center (aa)'].min()),
               xmax=[x+1 for x in contig_df.groupby('fragment contig')['fragment center (aa)'].max()],
               color='red',linewidth=5,alpha=.5)
    plt.hlines(y=list(cluster_df['height']),
               xmin=list(cluster_df['cluster first fragment center (aa)']),
               xmax=list(cluster_df['cluster last fragment center (aa)']+1),
               color='blue',linewidth=5,alpha=.5)
    plt.savefig(f"{name}.png",dpi=300)
    plt.show()
    plt.close('all')
    return


# Functions for calculating the benchmark statistics

def rangeOverlap(startA,endA,startB,endB,fragLen=30,fragOverlapReq=15):
    '''Quick check for overlap between the range of peak A and peak B
    
    Parameters
    - startA: center aa of first fragment in peak A
    - endA:   center aa of last fragment in peak A
    - startB: center aa of first fragment in peak B
    - endB:   center aa of last fragment in peak B
    - fragLen number of residues in each fragment
    - fragOverlapReq:  number of residues that must be shared between two fragments to be considered overlapping
    
    There are four scenarios that can occur when comparing two ranges
    1) endA < startB 
    2) startB <= endA <= endB
    3) startB <= startA <= endB
    4) endB < startA
    
    1 and 4 correspond to no overlap, so we can check for those and negate to find if there is overlap
    
    We are also interested in the case where the ranges don't overlap, but the fragments within them share
    some residues. To handle this, we artificially extend the range according to the fragLen and 
    fragOverlapReq parameters. This is only performed for one of the ranges, to avoid double-counting
    '''
    assert startA <= endA,(f"{startA},{endA}")
    assert startB <= endB,(f"{startB},{endB}")
    rangeExtension = fragLen - fragOverlapReq
    # check if there is no overlap, then negate
    return not ((endA + rangeExtension < startB) or (endB + rangeExtension < startA))

# print(rangeOverlap(0,10,10,20,30,30))
# print(rangeOverlap(0,10,11,20,30,30))
# print(rangeOverlap(10,20,1,20,30,30))
# print(rangeOverlap(10,20,1,2,30,30))

# print(rangeOverlap(0,10,25,30,30,15))
# print(rangeOverlap(0,10,26,30,30,15))

# print(rangeOverlap(100,130,130,160,30,30))
# print(rangeOverlap(100,130,131,161,30,30))

# print(rangeOverlap(100,130,150,180,30,20))
# print(rangeOverlap(100,130,140,170,30,20))

def getResidueOverlapReq(fragmentLength,residueOverlapFrac):
    return round(fragmentLength * residueOverlapFrac)

# fragment_length = 30
# residue_overlap_frac = float(2/3)
# residue_overlap_req = getResidueOverlapReq(fragment_length,residue_overlap_frac)
# print(residue_overlap_req)

def calculateOverlapBetweenPredandExp(pred_df,exp_df,residue_overlap_frac_list):
    # lists for storing data
    gene_list = []
    ## experimental peak columns
    exp_cluster_id_list = []
    exp_start_list = []
    exp_end_list = []
    exp_inhibeffect_list = []
    exp_type_list = []
    ## pred peak columns
    pred_cluster_id_list = []
    pred_start_list = []
    pred_end_list = []
    pred_clustersize_list = []
    pred_clusternfragments_list = []
    ## overlap
    fragmentlength_list = []
    overlap_frac_list = []
    overlap_res_list = []
    overlap_list = []
    
    print(f"Residue overlap fraction cutoffs = {residue_overlap_frac_list}")
    for residue_overlap_frac in residue_overlap_frac_list:
        for i,(_,exp_peak) in enumerate(exp_df.iterrows()):
            res_overlap = getResidueOverlapReq(exp_peak['fragment length (aa)'],residue_overlap_frac)
            for j,(_,pred_peak) in enumerate(pred_df.iterrows()):
                # print(exp_peak['peak region first fragment center (aa)'],exp_peak['peak region last fragment center (aa)'])

                # if the peaks come from different genes, or different lengths, they can't overlap
                if (exp_peak['gene'] != pred_peak['gene']) or (exp_peak['fragment length (aa)'] != pred_peak['fragment length (aa)']):
                    continue
                overlap = rangeOverlap(exp_peak['peak region first fragment center (aa)'],exp_peak['peak region last fragment center (aa)'],
                                           pred_peak['cluster first fragment center (aa)'],pred_peak['cluster last fragment center (aa)'],
                                           exp_peak['fragment length (aa)'],res_overlap)

                # append data to create new row
                gene_list.append(exp_peak['gene'])
                ## experimental peak columns
                exp_cluster_id_list.append(i)
                exp_start_list.append(exp_peak['peak region first fragment center (aa)'])
                exp_end_list.append(exp_peak['peak region last fragment center (aa)'])
                exp_inhibeffect_list.append(exp_peak['E = inhibitory effect (enrichment)'])
                exp_type_list.append(exp_peak['peak type'])
                ## pred peak columns
                pred_cluster_id_list.append(j)
                pred_start_list.append(pred_peak['cluster first fragment center (aa)'])
                pred_end_list.append(pred_peak['cluster last fragment center (aa)'])
                pred_clustersize_list.append(pred_peak['cluster_size'])
                pred_clusternfragments_list.append(pred_peak['cluster_n_fragments'])
                ## overlap
                fragmentlength_list.append(exp_peak['fragment length (aa)'])
                overlap_frac_list.append(residue_overlap_frac)
                overlap_res_list.append(res_overlap)
                overlap_list.append(overlap)


    overlap_df = pd.DataFrame({
        'gene':gene_list,
        'exp cluster id':exp_cluster_id_list,
        'exp cluster first':exp_start_list,
        'exp cluster last' :exp_end_list,
        'E = inhibitory effect (enrichment)':exp_inhibeffect_list,
        'peak type':exp_type_list,
        'pred cluster id':pred_cluster_id_list,
        'pred cluster first':pred_start_list,
        'pred cluster last':pred_end_list,
        'pred cluster size':pred_clustersize_list,
        'pred cluster n fragments':pred_clusternfragments_list,
        'fragment length (aa)':fragmentlength_list,
        'overlap fraction req':overlap_frac_list,
        'overlap req':overlap_res_list,
        'overlap':overlap_list
    })
    return overlap_df

def calculateBenchmarkStatistics(overlap_df,exp_df,maxFragmentLength,byGene=False):

    nExpPeaks = len(exp_df)
    nExpKnownPeaks = len(exp_df[exp_df['peak type']=='known'])

    if byGene:
        nExpPeaks_byGene_df = exp_df.groupby('gene').size().reset_index()
        nExpPeaks_byGene = dict(zip(nExpPeaks_byGene_df['gene'],nExpPeaks_byGene_df[0]))
        nExpKnownPeaks_byGene_df = exp_df[exp_df['peak type']=='known'].groupby('gene').size().reset_index()
        nExpKnownPeaks_byGene = dict(zip(nExpKnownPeaks_byGene_df['gene'],nExpKnownPeaks_byGene_df[0]))

        nExpPeaks_list = []
        nExpKnownPeaks_list = []
        genename_list = []

    minClusterSize_list = []
    resOverlapReq_list = []
    resOverlapFrac_list = []

    ntotalpredpeaks_list = []
    noverlapexppeaks_list = []
    noverlappredpeaks_list = []

    noverlapexppeaks_known_list = []
    noverlappredpeaks_known_list = []

    geneName = ""
    genesList = list(overlap_df.gene.unique()) if byGene else [""]
    print(f"Genes: {genesList}")
    for minClusterSize in range(1,15):
        for resOverlapReq in range(15,maxFragmentLength+1):
            for geneName in genesList:
                resOverlapFrac = resOverlapReq/maxFragmentLength
                if geneName == "":
                    filt_overlap_df = overlap_df[(overlap_df['pred cluster n fragments']>=minClusterSize)&
                                                (overlap_df['overlap fraction req']==resOverlapFrac)]
                else:
                    filt_overlap_df = overlap_df[(overlap_df['pred cluster n fragments']>=minClusterSize)&
                                                (overlap_df['overlap fraction req']==resOverlapFrac)&
                                                (overlap_df['gene']==geneName)]
                ntotalpredpeaks = filt_overlap_df.groupby('pred cluster id').ngroups
                noverlapexppeaks = filt_overlap_df[filt_overlap_df['overlap']==True].groupby('exp cluster id').ngroups
                noverlappredpeaks = filt_overlap_df[filt_overlap_df['overlap']==True].groupby('pred cluster id').ngroups
                
                minClusterSize_list.append(minClusterSize)
                resOverlapReq_list.append(resOverlapReq)
                resOverlapFrac_list.append(resOverlapFrac)
                
                ntotalpredpeaks_list.append(ntotalpredpeaks)
                noverlapexppeaks_list.append(noverlapexppeaks)
                noverlappredpeaks_list.append(noverlappredpeaks)
                if geneName != "":
                    genename_list.append(geneName)
                    nExpPeaks_list.append(nExpPeaks_byGene[geneName])
                    nExpKnownPeaks_list.append(nExpKnownPeaks_byGene[geneName])

                # Filter for "known" peaks and repeat the calculations
                filt_overlap_known_df = filt_overlap_df[filt_overlap_df['peak type']=='known']
                noverlapexpknownpeaks = filt_overlap_known_df[filt_overlap_known_df['overlap']==True].groupby('exp cluster id').ngroups
                noverlappredknownpeaks = filt_overlap_known_df[filt_overlap_known_df['overlap']==True].groupby('pred cluster id').ngroups

                noverlapexppeaks_known_list.append(noverlapexpknownpeaks)
                noverlappredpeaks_known_list.append(noverlappredknownpeaks)
            
    if len(genename_list) == 0:
        stat_df = pd.DataFrame({'min cluster size':minClusterSize_list,
                                'overlap frac req':resOverlapFrac_list,
                                'overlap req':resOverlapReq_list,
                                '# exp peaks':nExpPeaks,
                                '# pred peaks':ntotalpredpeaks_list,
                                '# overlap exp peaks':noverlapexppeaks_list,
                                '# overlap pred peaks':noverlappredpeaks_list,
                                '# known exp peaks':nExpKnownPeaks,
                                '# overlap known exp peaks':noverlapexppeaks_known_list,
                                '# overlap known pred peaks':noverlappredpeaks_known_list
                            })
    else:
        stat_df = pd.DataFrame({'gene':genename_list,
                                'min cluster size':minClusterSize_list,
                                'overlap frac req':resOverlapFrac_list,
                                'overlap req':resOverlapReq_list,
                                '# exp peaks':nExpPeaks_list,
                                '# pred peaks':ntotalpredpeaks_list,
                                '# overlap exp peaks':noverlapexppeaks_list,
                                '# overlap pred peaks':noverlappredpeaks_list,
                                '# known exp peaks':nExpKnownPeaks_list,
                                '# overlap known exp peaks':noverlapexppeaks_known_list,
                                '# overlap known pred peaks':noverlappredpeaks_known_list
                            })
    stat_df['% inhibitory peaks predicted'] = 100 * stat_df['# overlap exp peaks'] / stat_df['# exp peaks']
    stat_df['% known inhibitory peaks predicted'] = 100 * stat_df['# overlap known exp peaks'] / stat_df['# known exp peaks']
    stat_df['% predicted peaks inhibitory'] = 100 * stat_df['# overlap pred peaks'] / stat_df['# pred peaks']
    return stat_df
