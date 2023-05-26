import glob, json, argparse, os, sys, re, functools
import string
import multiprocessing as mp

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio.PDB import *

sys.path.insert(0,'/data1/groups/keatinglab/swans/savinovCollaboration/inhibitory_fragments_structure_prediction')
from src.colabfold_process_output_utils import *
from src.analyze_predictions import *

def getRankFromPath(path):
    pat = r"rank_(00)?(\d)"
    match = re.search(pattern=pat,string=path)
    if match is None:
        raise ValueError
    else:
        return int(match.group(2))
    

def load_confidence_data(path):
    name = path.split('/')[-3]
    pred_rank = getRankFromPath(path)
    start,end = name.split('_')[-1].split('-')
    start,end = int(start),int(end)
    center = (start+end)/2
    fragment_len = end - start + 1

    # Use the .a3m info line to determine the number of chains/residue lengths
    a3m_path_list = glob.glob(os.path.join(os.path.dirname(path),'*.a3m'))
    assert len(a3m_path_list) == 1, (print(a3m_path_list))
    n_copies_list,n_residue_list = getChainInfoFromA3M(a3m_path_list[0])
    # the following section creates a list with len == total number of residues, where the value is incremented by 1 between each chain
    chain_identifier = 0
    asym_id = []
    for i,n_residues_in_chain_x in enumerate(n_residue_list):
        for _ in range(n_copies_list[i]):
            for _ in range(n_residues_in_chain_x):
                asym_id.append(chain_identifier)
            chain_identifier+=1
    # asym_id = np.array([0]*306+[1]*30) previously set this manually
    asym_id = np.array(asym_id)
        
    with open(path,'r') as file:
        data = json.load(file)

    plddt = float(np.mean(np.array(data['plddt'][-fragment_len:])))        
    iptm = float(predicted_tm_score(np.array(data['pae']),asym_id,True))

    return [name,pred_rank,start,center,end,plddt,iptm]


def get_confidence_dataframe(colab_path,load_all=False,n_workers=1):

    if load_all:
        all_paths = glob.glob(os.path.join(colab_path,"data/*/output/*_rank_*.json"))
    else:
        all_paths = glob.glob(os.path.join(colab_path,"data/*/output/*_unrelaxed_rank_1_*.json"))
        if len(all_paths) == 0:
            all_paths = glob.glob(os.path.join(colab_path,"data/*/output/*scores_rank_001_*.json"))
        if len(all_paths) == 0:
            raise ValueError("Did not find .json files")
    
    if n_workers > 1:
        with mp.Pool(n_workers) as pool:
            data = pool.map(func=load_confidence_data,iterable=all_paths,chunksize=1)
    else:
        data = list()
        for path in all_paths:
            data.append(load_confidence_data(path))

    confidence_df = pd.DataFrame(data,columns=['fragment_name','rank','fragment start (aa)','fragment center (aa)','fragment end (aa)','iptm','plddt'])
    confidence_df = confidence_df.sort_values(by='fragment start (aa)')
    print(f"{len(all_paths)} total paths and {len(confidence_df)} lines in the dataframe")
    return confidence_df

# pool.map only accepts a function and an iterable
# since requires other fixed parameters, we bake those into the function using a decorator

# def load_contact_data_decorator(protein_chains,fragment_chains,distance_cutoff):
#     def load_contact_data_wrapped(path):
#         return load_contact_data(path,protein_chains,fragment_chains,distance_cutoff)
#     return load_contact_data_wrapped

def load_contact_data(path,protein_chains,fragment_chains,distance_cutoff):
    name = path.split('/')[-3]
    pred_rank = getRankFromPath(path)
    n_conts = countInterfaceContacts(path,protein_chains,fragment_chains,distance_cutoff)
    start,end = name.split('_')[-1].split('-')
    start,end = int(start),int(end)
    center = (start + end) / 2
    return name,pred_rank,start,center,end,n_conts,path

def get_contact_dataframe(colab_path,contact_distance_cutoff,load_all=False,n_workers=1):

    if load_all:
        all_paths = glob.glob(os.path.join(colab_path,"data/*/output/*_rank_*.pdb"))
    else:
        all_paths = glob.glob(os.path.join(colab_path,'data/*/output/*_unrelaxed_rank_1_*.pdb'))
        if len(all_paths)==0:
            all_paths = glob.glob(os.path.join(colab_path,'data/*/output/*_unrelaxed_rank_001_*.pdb'))
        if len(all_paths)==0:
            raise ValueError("Did not find .pdb files")
        
    # Use the .a3m info line to determine the number of chains/residue lengths (will be the same for each individual structure)
    a3m_path_list = glob.glob(os.path.join(os.path.dirname(all_paths[0]),'*.a3m'))
    assert len(a3m_path_list) == 1
    n_copies_list,_ = getChainInfoFromA3M(a3m_path_list[0])
    n_total_chains = np.sum(np.array([x for x in n_copies_list]))
    assert n_total_chains < len(string.ascii_uppercase), "ran out of chain names, are you sure there's not a bug?"
    protein_chains = set([chain_id for i,chain_id in enumerate(string.ascii_uppercase) if i < n_total_chains - 1])
    fragment_chains = set(string.ascii_uppercase[n_total_chains-1])
    
    if n_workers > 1:
        load_contact_data_mapper = functools.partial(load_contact_data,protein_chains=protein_chains,fragment_chains=fragment_chains,distance_cutoff=contact_distance_cutoff)
        with mp.Pool(n_workers) as pool:
            data = pool.map(func=load_contact_data_mapper,iterable=all_paths,chunksize=1)
    else:
        data = list()
        for path in all_paths:
            data.append(load_contact_data(path,protein_chains,fragment_chains,distance_cutoff=contact_distance_cutoff))

    conts_df = pd.DataFrame(data,columns=['fragment_name','rank','fragment start (aa)','fragment center (aa)','fragment end (aa)','n_contacts','path'])
    conts_df = conts_df.sort_values(by='fragment start (aa)')
    print(f"{len(all_paths)} total paths and {len(conts_df)} lines in the dataframe")
    return conts_df


def main(args):
    print('Loading JSON file specifying where colabfold results are located')
    # Load JSON file specifying where to import colabfold results from
    json_path = args.import_json
    load_all = args.load_all
    # n_workers = mp.cpu_count() 
    n_workers = len(os.sched_getaffinity(0))
    contact_distance_cutoff = args.contact_distance_cutoff
    experimental_data_path = args.experimental_data

    print(f"Loading data with {n_workers} workers")

    with open(json_path,"r") as file:
        colab_results = json.loads(file.read())

    # For each gene and condition, import available data and create individual dataframes
    print('Processing results...')
    df_list = []
    merge_on_list = ['fragment_name','rank','fragment start (aa)','fragment center (aa)','fragment end (aa)']
    for gene_name in colab_results:
        print('gene_name:',gene_name)
        for condition in colab_results[gene_name]:
            print('condition:',condition)

            # Load confidence (pLDDT/iPTM)
            confidence_df = get_confidence_dataframe(colab_results[gene_name][condition]['colabfold'],load_all,n_workers)

            # Count contacts
            contacts_df = get_contact_dataframe(colab_results[gene_name][condition]['colabfold'],contact_distance_cutoff,load_all,n_workers)

            # TODO calculate RMSD

            # Merge dataframes into a single one containing all types of data
            comb_df = confidence_df.merge(contacts_df,on=merge_on_list,indicator=True,validate="one_to_one")

            # Calculate the weighted contacts
            comb_df['weighted_contacts'] = comb_df['n_contacts'] * comb_df['iptm']/100

            comb_df['gene'] = gene_name
            comb_df['condition'] = condition

            print(f"Combined dataframe with {len(comb_df)} rows")

            df_list.append(comb_df)

    # Concatenate all dataframes into a single dataframe that stores all the data. Merge with experimental data df 
    concat_df = pd.concat(df_list,ignore_index=True)
    concat_df.to_csv("colabfold_predictions.csv")
    print(f"Dataframe with {len(concat_df)} entries total")

    # If experimental dataframe is provided, merge now
    if experimental_data_path != "":
        exp_df = pd.read_csv(experimental_data_path)
        print(f"Loaded experimental dataframe")
        merge_on_list = ['gene','fragment start (aa)','fragment center (aa)','fragment end (aa)']
        merge_df = concat_df.merge(exp_df,how='left',on=merge_on_list,validate='many_to_one')
        merge_df.to_csv("colabfold_predictions_expmerge.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'ColabFoldProcessOutput',
        description = 'Script for processing the output of colabfold jobs and combining into a single dataframe')
    parser.add_argument('--import_json',required=True)
    parser.add_argument('--experimental_data',required=False,default="")
    parser.add_argument('--load_all',required=False,action='store_true',
                        help='Load all predictions (not just top-ranked)')
    # parser.add_argument('--n_workers',required=False,default=1,type=int,
    #                     help='The number of processes that can be run concurrently (set this to the number of available CPU cores)')
    parser.add_argument('--contact_distance_cutoff',required=False,default=3.5,type=float,
                        help='The distance cutoff between heavy atoms of interface residues that defines a contact')

    args = parser.parse_args()
    main(args)

    print('Done!')