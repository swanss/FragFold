import glob, json, argparse, os, sys, re, functools
import string
import multiprocessing as mp

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio.PDB import *

from fragfold.src.colabfold_process_output_utils import *
from fragfold.src.analyze_predictions import *

def getRankFromPath(path):
    pat = r"rank_(00)?(\d)"
    match = re.search(pattern=pat,string=path)
    if match is None:
        raise ValueError
    else:
        return int(match.group(2))
    

def load_confidence_data(path):
    name = path.split('/')[-3]
    start,end = name.split('_')[-1].split('-')
    start,end = int(start),int(end)
    center = (start+end)/2

    # e.g. will match 2023-05-26 20:03:23,691 rank_001_alphafold2_ptm_model_1_seed_000 pLDDT=93 pTM=0.801 ipTM=0.257 
    pat = r"rank_00(\d)_alphafold2_ptm_model_\d_seed_\d{3} pLDDT=([+-]?[0-9]*[.]?[0-9]+) pTM=([+-]?[0-9]*[.]?[0-9]+) ipTM=([+-]?[0-9]*[.]?[0-9]+)"
    
    confidence_dict = {'fragment_name':[],'rank':[],
                       'fragment start (aa)':[],'fragment center (aa)':[],'fragment end (aa)':[],
                       'plddt':[],'ptm':[],'iptm':[]}

    match_count = 0
    with open(path,'r') as file:
        for line in file:
            match = re.search(pattern=pat,string=line)
            if match is not None:
                match_count+=1
                pred_rank = int(match.group(1))
                plddt = float(match.group(2))
                ptm = float(match.group(3))
                iptm = float(match.group(4))

                confidence_dict['fragment_name'].append(name)
                confidence_dict['rank'].append(pred_rank)
                confidence_dict['fragment start (aa)'].append(start)
                confidence_dict['fragment center (aa)'].append(center)
                confidence_dict['fragment end (aa)'].append(end)
                confidence_dict['plddt'].append(plddt)
                confidence_dict['ptm'].append(ptm)
                confidence_dict['iptm'].append(iptm)
    # if match_count < 5:
    #     raise ValueError(f"Expected to find 5 matching lines, instead found {match_count}")

    return pd.DataFrame(confidence_dict)


def get_confidence_dataframe(colab_path,n_workers=1):

    all_paths = glob.glob(os.path.join(colab_path,"data/*/output/log.txt"))

    if n_workers > 1:
        with mp.Pool(n_workers) as pool:
            data = pool.map(func=load_confidence_data,iterable=all_paths,chunksize=1)
    else:
        data = list()
        for path in all_paths:
            data.append(load_confidence_data(path))

    confidence_df = pd.concat(data,ignore_index=True)
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

def get_contact_dataframe(colab_path,contact_distance_cutoff,n_workers=1):

    all_paths = glob.glob(os.path.join(colab_path,'data/*/output/*_unrelaxed_rank_00?_*.pdb'))
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
    conts_df['protein_chains'] = ','.join(protein_chains)
    conts_df['fragment_chain'] = ','.join(fragment_chains)
    conts_df = conts_df.sort_values(by='fragment start (aa)')
    print(f"{len(all_paths)} total paths and {len(conts_df)} lines in the dataframe")
    return conts_df


def main(args):
    # Load JSON file specifying where to import colabfold results from
    json_path = args.import_json
    print(f"Loading JSON file specifying where colabfold results are located: {json_path}")
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
            confidence_df = get_confidence_dataframe(colab_results[gene_name][condition]['colabfold'],n_workers)

            # Count contacts
            contacts_df = get_contact_dataframe(colab_results[gene_name][condition]['colabfold'],contact_distance_cutoff,n_workers)

            # TODO calculate RMSD

            # Merge dataframes into a single one containing all types of data
            comb_df = confidence_df.merge(contacts_df,on=merge_on_list,indicator=True,validate="one_to_one")

            # Calculate the weighted contacts
            comb_df['weighted_contacts'] = comb_df['n_contacts'] * comb_df['iptm']

            comb_df['gene'] = gene_name
            comb_df['condition'] = condition

            print(f"Combined dataframe with {len(comb_df)} rows")

            df_list.append(comb_df)

    # Concatenate all dataframes into a single dataframe that stores all the data. Merge with experimental data df 
    concat_df = pd.concat(df_list,ignore_index=True)
    concat_df['fragment length (aa)'] = concat_df['fragment end (aa)'] - concat_df['fragment start (aa)'] + 1
    concat_df.to_csv("colabfold_predictions.csv")
    print(f"Dataframe with {len(concat_df)} entries total")

    # If experimental dataframe is provided, merge now
    if experimental_data_path != "":
        exp_df = pd.read_csv(experimental_data_path)
        print(f"Loaded experimental dataframe")
        merge_on_list = ['gene','fragment start (aa)','fragment center (aa)','fragment end (aa)','fragment length (aa)']
        # dataframe contains some replicated measurements (they are generally similar, so arbitrarily take the first)
        exp_df = exp_df.drop_duplicates(subset=merge_on_list)
        merge_df = concat_df.merge(exp_df,how='left',on=merge_on_list,validate='many_to_one')
        merge_df.to_csv("colabfold_predictions_expmerge.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'ColabFoldProcessOutput',
        description = 'Script for processing the output of colabfold jobs and combining into a single dataframe')
    parser.add_argument('--import_json',required=True)
    parser.add_argument('--experimental_data',required=False,default="")
    # parser.add_argument('--n_workers',required=False,default=1,type=int,
    #                     help='The number of processes that can be run concurrently (set this to the number of available CPU cores)')
    parser.add_argument('--contact_distance_cutoff',required=False,default=3.5,type=float,
                        help='The distance cutoff between heavy atoms of interface residues that defines a contact')

    args = parser.parse_args()
    main(args)

    print('Done!')