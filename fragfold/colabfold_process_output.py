from pathlib import Path
from datetime import datetime
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
from fragfold.src.plot_utils import *

def getRankFromPath(path):
    pat = r"rank_(00)?(\d)"
    match = re.search(pattern=pat,string=path)
    if match is None:
        raise ValueError
    else:
        return int(match.group(2))
    

def load_confidence_data(path):
    # e.g. will match 2024-05-04 12:32:21,588 Query 1/1: ftsZ1copies_10-316_ftsZ_166-195 (length 337)
    fragment_name_pat = r"Query 1\/1: ((\D+)(\d)copies_(\d+)-(\d+)_(\D+)_(\d+)-(\d+)) \(length (\d+)\)"

    # e.g. will match 2023-05-26 20:03:23,691 rank_001_alphafold2_ptm_model_1_seed_000 pLDDT=93 pTM=0.801 ipTM=0.257 
    confidence_pat = r"rank_00(\d)_alphafold2_ptm_model_\d_seed_\d{3} pLDDT=([+-]?[0-9]*[.]?[0-9]+) pTM=([+-]?[0-9]*[.]?[0-9]+) ipTM=([+-]?[0-9]*[.]?[0-9]+)"
    
    confidence_dict = {'fragment name':[],'rank':[],
                       'fragment start (aa)':[],'fragment center (aa)':[],'fragment end (aa)':[],
                       'plddt':[],'ptm':[],'iptm':[]}

    fragment_name_match_count,conf_match_count = 0,0
    with open(path,'r') as file:
        for i,line in enumerate(file):
            # search for fragment name
            match = re.search(pattern=fragment_name_pat,string=line)
            if match is not None:
                fragment_name_match_count += 1
                name = match[1]
                start,end = int(match[7]),int(match[8])
                center = (start+end)/2
                continue

            # search for confidence data
            match = re.search(pattern=confidence_pat,string=line)
            if match is not None:
                if fragment_name_match_count == 0:
                    raise ValueError('Fragment name should come before the confidence metrics')
                conf_match_count+=1
                pred_rank = int(match.group(1))
                plddt = float(match.group(2))
                ptm = float(match.group(3))
                iptm = float(match.group(4))

                confidence_dict['fragment name'].append(name)
                confidence_dict['rank'].append(pred_rank)
                confidence_dict['fragment start (aa)'].append(start)
                confidence_dict['fragment center (aa)'].append(center)
                confidence_dict['fragment end (aa)'].append(end)
                confidence_dict['plddt'].append(plddt)
                confidence_dict['ptm'].append(ptm)
                confidence_dict['iptm'].append(iptm)
    if conf_match_count != 5 or fragment_name_match_count != 1:
        raise ValueError(f"Expected to find 5 lines with confidence metrics and 1 with the fragment name, \
                         instead found ({conf_match_count},{fragment_name_match_count}) when searching {path}")

    return pd.DataFrame(confidence_dict)


def get_confidence_dataframe(all_paths,n_workers=1):

    if len(all_paths) <= 0:
        raise ValueError(f"Unable to find colabfold output in the directory: {glob_path}")

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
    name = path.split('_unrelaxed')[0]
    pred_rank = getRankFromPath(path)
    n_conts = countInterfaceContacts(path,protein_chains,fragment_chains,distance_cutoff)
    start,end = name.split('_')[-1].split('-')
    start,end = int(start),int(end)
    center = (start + end) / 2
    return name,pred_rank,start,center,end,n_conts,path

def get_contact_dataframe(all_paths,contact_distance_cutoff,n_workers=1):

    # if len(all_paths)==0:
    #     raise ValueError("Did not find .pdb files")
        
    # # Use the .a3m info line to determine the number of chains/residue lengths (will be the same for each individual structure)
    # a3m_path_list = glob.glob(os.path.join(os.path.dirname(all_paths[0]),'*.a3m'))
    # assert len(a3m_path_list) > 0
    # n_copies_list,_ = getChainInfoFromA3M(a3m_path_list[0])
    # n_total_chains = np.sum(np.array([x for x in n_copies_list]))
    # assert n_total_chains < len(string.ascii_uppercase), "ran out of chain names, are you sure there's not a bug?"
    # protein_chains = set([chain_id for i,chain_id in enumerate(string.ascii_uppercase) if i < n_total_chains - 1])
    # fragment_chains = set(string.ascii_uppercase[n_total_chains-1])
    protein_chains = set('A')
    fragment_chains = set('B')
    
    if n_workers > 1:
        load_contact_data_mapper = functools.partial(load_contact_data,protein_chains=protein_chains,fragment_chains=fragment_chains,distance_cutoff=contact_distance_cutoff)
        with mp.Pool(n_workers) as pool:
            data = pool.map(func=load_contact_data_mapper,iterable=all_paths,chunksize=1)
    else:
        data = list()
        for path in all_paths:
            data.append(load_contact_data(path,protein_chains,fragment_chains,distance_cutoff=contact_distance_cutoff))

    conts_df = pd.DataFrame(data,columns=['fragment name','rank','fragment start (aa)','fragment center (aa)','fragment end (aa)','n_contacts','path'])
    conts_df['protein_chains'] = ','.join(protein_chains)
    conts_df['fragment_chain'] = ','.join(fragment_chains)
    conts_df = conts_df.sort_values(by='fragment start (aa)')
    print(f"{len(all_paths)} total paths and {len(conts_df)} lines in the dataframe")
    return conts_df


def main(args):
    # # Load JSON file specifying where to import colabfold results from
    # json_path = Path(args.import_json)
    # assert json_path.is_file()
    # print(f"Loading JSON file specifying where colabfold results are located: {json_path}")
    # n_workers = mp.cpu_count() 
    n_workers = len(os.sched_getaffinity(0))
    contact_distance_cutoff = args.contact_distance_cutoff
    experimental_data_path = args.experimental_data

    print(f"Loading data with {n_workers} workers")

    # with open(json_path,"r") as file:
    #     colab_results = json.loads(file.read())

    # For each gene and condition, import available data and create individual dataframes
    print('Processing results...')
    df_list = []
    merge_on_list = ['fragment name','rank','fragment start (aa)','fragment center (aa)','fragment end (aa)']
    fragment_protein_name = args.fragment_protein
    full_protein_name = args.full_protein

    # Load confidence (pLDDT/iPTM)
    confidence_log_paths = args.confidence_logs
    print(confidence_log_paths)
    confidence_df = get_confidence_dataframe(confidence_log_paths,n_workers)

    # Count contacts
    predicted_pdb_paths = args.predicted_pdbs
    print(predicted_pdb_paths)
    contacts_df = get_contact_dataframe(predicted_pdb_paths,contact_distance_cutoff,n_workers)

    # TODO calculate RMSD

    # Merge dataframes into a single one containing all types of data
    comb_df = confidence_df.merge(contacts_df,on=merge_on_list,indicator=True,validate="one_to_one")

    # Calculate the weighted contacts
    comb_df['weighted_contacts'] = comb_df['n_contacts'] * comb_df['iptm']

    comb_df['fragment_parent_name'] = fragment_protein_name
    comb_df['protein_name'] = full_protein_name
    comb_df['fragment length (aa)'] = comb_df['fragment end (aa)'] - comb_df['fragment start (aa)'] + 1
    print(f"Combined dataframe with {len(comb_df)} rows")

    # Merge with experimental data df 
    comb_df.to_csv("colabfold_predictions.csv")
    print(f"Dataframe with {len(comb_df)} entries total")
    # confidence_df.to_csv("colabfold_predictions.csv")

    # If experimental dataframe is provided, merge now
    if experimental_data_path != "":
        exp_df = pd.read_csv(experimental_data_path)
        print(f"Loaded experimental dataframe")
        merge_on_list = ['gene','fragment start (aa)','fragment center (aa)','fragment end (aa)','fragment length (aa)']
        # dataframe contains some replicated measurements (they are generally similar, so arbitrarily take the first)
        exp_df = exp_df.drop_duplicates(subset=merge_on_list)
        merge_df = comb_df.merge(exp_df,how='left',on=merge_on_list,validate='many_to_one')
        merge_df.to_csv(f"results_expmerge.csv")

    if args.generate_plots:
        for gene,group_df in comb_df.groupby('gene'):
            ax = plotRawValuesOnSingle(group_df)
            plt.savefig(f"{datetime.today().strftime('%y%m%d')}_fragmentcenter_vs_weightedcontacts_combined_fragments-{gene}.png",
                        dpi=300,bbox_inches='tight')
            g = plotRawValuesOnFacetGrid(group_df)
            plt.savefig(f"{datetime.today().strftime('%y%m%d')}_fragmentcenter_vs_weightedcontacts_facetgrid_fragments-{gene}.png",
                        dpi=300,bbox_inches='tight')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'ColabFoldProcessOutput',
        description = 'Script for processing the output of colabfold jobs and combining into a single dataframe')
    parser.add_argument('--predicted_pdbs',nargs='+',required=True,
                        help='List of paths to all of the PDBs predicted by ColabFold')
    parser.add_argument('--confidence_logs',nargs='+',required=True,
                        help='List of log.txt files from the ColabFold Predictions')
    parser.add_argument('--full_protein',required=True)
    parser.add_argument('--fragment_protein',required=True)
    parser.add_argument('--experimental_data',required=False,default="")
    parser.add_argument('--contact_distance_cutoff',required=False,default=3.5,type=float,
                        help='The distance cutoff between heavy atoms of interface residues that defines a contact')
    parser.add_argument('--generate_plots',action='store_true',
                        help='If provided, will generate plots for each fragmented gene + condition')
    args = parser.parse_args()
    main(args)

    print('Done!')