import glob, json, math, pathlib
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

from fragfold.src.analyze_predictions import *
from fragfold.src.peak_prediction import filterAlphaFoldPredictions,splitDiscontinuousFragmentSets,clusterOverlappingFragments,organizeClusters,plotClusters,clusterPeaksByOverlap

def loadParamsSetFromJSON(path,verbose=True):
    # Define parameter ranges
    # Load params from JSON file
    with open(path,"r") as file:
        params_dict = json.load(file)

    ncontacts_scan = np.arange(params_dict['n_contacts']['start'],
                               params_dict['n_contacts']['stop'],
                               params_dict['n_contacts']['step'])
    nwcontacts_scan = np.arange(params_dict['weighted_contacts']['start'],
                                params_dict['weighted_contacts']['stop'],
                                params_dict['weighted_contacts']['step'])
    iptm_scan = np.linspace(params_dict['iptm']['start'],
                            params_dict['iptm']['stop'],
                            params_dict['iptm']['n'])
    contactsimcutoff_scan = np.linspace(params_dict['cluster_dist_cutoff']['start'],
                                        params_dict['cluster_dist_cutoff']['stop'],
                                        params_dict['cluster_dist_cutoff']['n'])

    print(f"ncontacts: {ncontacts_scan}")
    print(f"weighted_contacts: {nwcontacts_scan}")
    print(f"iptm: {iptm_scan}")
    print(f"cluster_dist_cutoff: {contactsimcutoff_scan}")
    allparams_scan = list(itertools.product(ncontacts_scan,nwcontacts_scan,iptm_scan,contactsimcutoff_scan))
    print(f"{len(allparams_scan)} total parameter conditions")
    return allparams_scan

def createName(n_contacts,n_wcontacts,iptm,contactsimcutoff):
    return f"ncontacts{n_contacts}_nwcontacts{n_wcontacts}_iptm{iptm:.2f}_contactdist{contactsimcutoff:.2f}"

def clusterWithParams(params:tuple,pred_df:pd.DataFrame,outdir=''):
    n_contacts,n_wcontacts,iptm,contactsimcutoff = params[0],params[1],params[2],params[3]
    print(f"Cluster with parameters: n_contacts = {n_contacts}, n_weighted_contacts = {n_wcontacts}, iptm = {iptm}, cluster contact cutoff = {contactsimcutoff}")
    # print(n_contacts,n_wcontacts,iptm,contactsimcutoff)
    # Filter all fragments
    filt_pred_df = filterAlphaFoldPredictions(pred_df,n_contacts,n_wcontacts,iptm,True)

    if len(filt_pred_df) < 0:
        return

    paramsname = createName(n_contacts,n_wcontacts,iptm,contactsimcutoff)
    dirname = os.path.join(outdir,paramsname)
    pathlib.Path(dirname).mkdir(exist_ok=True)

    # Group by FragFold job
    df_list = []
    grouped_df = filt_pred_df.groupby(['fragment_parent_name','protein_name','fragment_length_aa','description'])
    print(f"After filtering, there are {grouped_df.ngroups} groups in the dataframe")
    for (fragment_parent_name,protein_name,fragmentlen,description),group_df in grouped_df:
        if args.verbose:
            print(f"Clustering {len(group_df)} predictions from {fragment_parent_name}_{protein_name}_{fragmentlen}_{description} ...")

        # Find discontinuous sets of fragments (i.e. fragments that don't share residues)
        contigs_df = splitDiscontinuousFragmentSets(group_df,args.verbose) # modify df in place

        # Cluster each set
        grouped_clusters_list = []
        for contign,single_contig_df in contigs_df.groupby('fragment contig'):
            cluster_df = clusterOverlappingFragments(single_contig_df,contactsimcutoff,'',False)
            cluster_df['fragment contig'] = contign
            cluster_df['full_cluster_name'] = cluster_df['fragment contig'].astype(str) + '_' + cluster_df['cluster'].astype(str)

            df_list.append(cluster_df)
            grouped_clusters_list.append(cluster_df.copy(deep=True))

        grouped_comb_df = pd.concat(grouped_clusters_list,ignore_index=True)
        # Plot the clusters
        cluster_plot_name =  os.path.join(dirname,f"{paramsname}_{fragment_parent_name}_{protein_name}_{fragmentlen}")
        all_fragments_df = pred_df[(pred_df['fragment_parent_name']==fragment_parent_name)
                                   &(pred_df['protein_name']==protein_name)
                                   &(pred_df['fragment_length_aa']==fragmentlen)
                                   &(pred_df['description']==description)]
        plotClusters(all_fragments_df,grouped_comb_df,contigs_df,cluster_plot_name)
        
    if len(df_list) > 0:
        comb_df = pd.concat(df_list,ignore_index=True)
        comb_df['n_contacts_cutoff'] = n_contacts
        comb_df['n_weighted_contacts_cutoff'] = n_wcontacts
        comb_df['iptm_cutoff'] = iptm
        comb_df['contact_distance_cutoff'] = contactsimcutoff
        print(f"Found {len(comb_df)} clusters")
        return comb_df
    else:
        print("Found 0 clusters")
        return pd.DataFrame()

def singleParamSet(args):
    print("Single set of parameters mode")

    if (args.n_contacts is None or args.n_weighted_contacts is None 
        or args.iptm is None or args.contact_distance is None):
        raise ValueError("If one parameter is set, all parameters must be set.")
    
    # Load DF
    pred_df = pd.read_csv(args.colabfold_data_csv,index_col=0)
    pred_df['fragment_length_aa'] = pred_df['fragment_end_aa'] - pred_df['fragment_start_aa'] + 1
    name = createName(args.n_contacts,args.n_weighted_contacts,args.iptm,args.contact_distance)

    # Define parameters
    outdir="cluster_info"
    pathlib.Path(outdir).mkdir(exist_ok=True)
    comb_df = clusterWithParams((args.n_contacts,args.n_weighted_contacts,args.iptm,args.contact_distance),pred_df,outdir)
    comb_df.to_csv(f"predictalphafoldpeaks_{name}.csv")

    if args.cluster_peaks_frac_overlap > 0.0 and args.cluster_peaks_frac_overlap < 1.0:
        print("Merging overlapping peaks...")
        clus_filt_pred_df = clusterPeaksByOverlap(comb_df,frac_overlap=args.cluster_peaks_frac_overlap,verbose=False)
        print(f"After merging, {len(clus_filt_pred_df)} peaks remain...")
        clus_filt_pred_df.to_csv(f"predictalphafoldpeaks_mergeoverlapping{args.cluster_peaks_frac_overlap:.2f}_{name}.csv")

def paramScan(args):
    print("Parameter scan")

    # Load DF
    pred_df = pd.read_csv(args.colabfold_data_csv,index_col=0)
    pred_df['fragment_length_aa'] = pred_df['fragment_end_aa'] - pred_df['fragment_start_aa'] + 1

    # Define parameter ranges
    allparams_scan = loadParamsSetFromJSON(args.paramscan_json)
    
    if args.n_batches is None or args.batch_id == 1:
        with open("paramscaninfo.txt","w") as file:
            file.write(args.paramscan_json+"\n")

    # If there are batches, split up the param settings list
    if args.n_batches is not None and args.batch_id is not None:
        print(f"{args.n_batches} total batches and ID = {args.batch_id}. Original list length = {len(allparams_scan)}")
        # shuffle list so try to even out the runtime for each batch
        print("Shuffling the parameter sets to try to even out the runtime per batch")
        random.Random(42).shuffle(allparams_scan)
        batch_size = math.ceil(len(allparams_scan) / args.n_batches)
        allparams_scan = allparams_scan[args.batch_id*batch_size:(args.batch_id+1)*batch_size]
        print(f"New list length = {len(allparams_scan)}")

    if allparams_scan:
        print("No parameters given, terminating...")

    # Create all possible combinations of the values
    df_list = []
    # n_workers = min(len(os.sched_getaffinity(0)),48)
    n_workers = args.n_processes
    print(f"N workers: {n_workers}")
    outdir="cluster_info"
    pathlib.Path(outdir).mkdir(exist_ok=True)
    if n_workers > 1:
        clusterWithParams_mapper = functools.partial(clusterWithParams,pred_df=pred_df,outdir=outdir)
        with mp.Pool(n_workers) as pool:
            df_list = pool.map(func=clusterWithParams_mapper,iterable=allparams_scan,chunksize=1)
    else:
        df_list = list()
        for params in allparams_scan:
            df_list.append(clusterWithParams(params,pred_df,outdir=outdir))

    if (len(df_list)>0):
        comb_df = pd.concat(df_list)
        comb_df.to_csv(f"predictalphafoldpeaks_paramscan_batch{args.batch_id}.csv")

    if args.cluster_peaks_frac_overlap > 0.0 and args.cluster_peaks_frac_overlap < 1.0:
        print("Merging overlapping peaks...")
        clus_filt_pred_df = clusterPeaksByOverlap(comb_df,frac_overlap=args.cluster_peaks_frac_overlap,verbose=False)
        print(f"After merging, {len(clus_filt_pred_df)} peaks remain...")
        clus_filt_pred_df.to_csv(f"predictalphafoldpeaks_mergeoverlapping{args.cluster_peaks_frac_overlap:.2f}_paramscan_batch{args.batch_id}.csv")

    n_peak_list = list()
    for df in df_list:
        n_peak_list.append(len(df))
    ncontsl,nweightl,iptml,contactdistl = zip(*allparams_scan)
    npeaks_df = pd.DataFrame({'n_contacts_cutoff':ncontsl,
                              'n_weighted_contacts_cutoff':nweightl,
                              'iptm_cutoff':iptml,
                              'contactdist':contactdistl,
                              'n_peaks':n_peak_list})
    npeaks_df.to_csv(f"npeaks_paramscan_batch{args.batch_id}.csv")

def cleanUp(args):
    # Load all csv files and concatenate
    csv_files = glob.glob("npeaks_paramscan_batch*.csv")
    assert len(csv_files) > 0
    npeaks_df = pd.concat([pd.read_csv(path,index_col=0) for path in csv_files])

    # Load the params file and recreate the parameter set
    with open("paramscaninfo.txt") as file:
        paramsjsonpath = file.readline().rstrip()

    allparams_scan = loadParamsSetFromJSON(args.paramscan_json)

    n_param_conditions = len(allparams_scan)
    n_conditions_in_df = npeaks_df.groupby(['n_contacts_cutoff','n_weighted_contacts_cutoff','iptm_cutoff','contact_distance_cutoff']).ngroups
    if len(allparams_scan) != n_conditions_in_df:
        print(f"The number of parameter conditions in the n_peaks dataframe {n_conditions_in_df} is less than the total number of conditions considered ({n_param_conditions})")
    
    missing_count = 0
    for i,(nconts,nweight,iptm,clustercontact) in allparams_scan:
        n = len(npeaks_df[(npeaks_df['n_contacts_cutoff']==nconts)&
                     (npeaks_df['n_weighted_contacts_cutoff']==nweight)&
                     (npeaks_df['iptm_cutoff']==iptm)&
                     (npeaks_df['contact_distance_cutoff']==clustercontact)])
        if n == 1:
            pass
        elif n < 1:
            missing_count+=1
            print(f"Missing parameters: n_contacts = {nconts}, n_weighted_contacts = {nweight}, iptm = {iptm}, cluster contact cutoff = {clustercontact}")
        else:
            raise ValueError("Shouldn't be possible: are batches overlapping?")
    
    print(f"In total {missing_count} conditions are missing")

def main(args):
    if (args.colabfold_data_csv == ""):
        assert args.n_batches == 1
        print("No input data was passed, will instead cleanup the directory by merging the output into a single file")
        cleanUp(args)
    else:
        print("Predicting peaks...")
        if (args.paramscan_json is not None):
            paramScan(args)
        else:
            singleParamSet(args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'predictAlphaFoldPeaks',
        description = 'Script for predicting inhibitory fragments using AlphaFold. The parameters used for filtering fragments and clustering are defined in the JSON file. Includes options for parallelizing, which are useful if performing a parameter scan.')
    parser.add_argument('--n_batches',type=int,default=1)
    parser.add_argument('--batch_id',type=int,default=0)
    parser.add_argument('--n_processes',type=int,default=1)
    parser.add_argument('--n_contacts',type=int,)
    parser.add_argument('--n_weighted_contacts',type=int)
    parser.add_argument('--iptm',type=float)
    parser.add_argument('--contact_distance',type=float)
    parser.add_argument('--paramscan_json',type=str)
    parser.add_argument('--colabfold_data_csv',type=str)
    parser.add_argument('--cluster_peaks_frac_overlap',type=float,default=1.0)
    parser.add_argument('--verbose',action='store_true')

    args = parser.parse_args()
    print(args)
    main(args)
    print('Done!')