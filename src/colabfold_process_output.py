import glob, json, argparse, os, sys
import string

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio.PDB import *

sys.path.insert(0,'/data1/groups/keatinglab/swans/savinovCollaboration/inhibitory_fragments_structure_prediction')
from src.colabfold_process_output_utils import *
from src.colabfold_process_output import *
from src.analyze_predictions import *

def get_confidence_dataframe(colab_path):
    iptm_list = []
    plddt_list = []
    fragment_name = []
    fragment_start = []
    fragment_center = []
    fragment_end = []

    all_paths = glob.glob(os.path.join(colab_path,"data/*/output/*_unrelaxed_rank_1_*.json"))
    if len(all_paths) == 0:
        all_paths = glob.glob(os.path.join(colab_path,"data/*/output/*scores_rank_001_*.json"))
    if len(all_paths) == 0:
        raise ValueError("Did not find .json files")
    for path in all_paths:
        name = path.split('/')[-3]
        (start,end) = name.split('_')[-1].split('-')
        fragment_len = int(end) - int(start) + 1

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
        plddt_list.append(plddt)
            
        iptm = float(predicted_tm_score(np.array(data['pae']),asym_id,True))
        iptm_list.append(iptm)
        
        fragment_name.append(name)
        fragment_start.append(int(start))
        fragment_end.append(int(end))
        fragment_center.append((fragment_start[-1]+fragment_end[-1])/2)

    confidence_df = pd.DataFrame({
        'fragment_name':fragment_name,
        'fragment start (aa)':fragment_start,
        'fragment center (aa)':fragment_center,
        'fragment end (aa)':fragment_end,
        'iptm':iptm_list,
        'plddt':plddt_list
    })
    confidence_df = confidence_df.sort_values(by='fragment start (aa)')
    return confidence_df

def get_contact_dataframe(colab_path):
    fragment_name = []
    fragment_start = []
    fragment_end = []
    fragment_center = []
    n_contacts = []

    parser = PDBParser()

    all_paths = glob.glob(os.path.join(colab_path,'data/*/output/*_unrelaxed_rank_1_*.pdb'))
    if (len(all_paths)==0):
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

    for path in all_paths:
        n_conts = countInterfaceContacts(path,protein_chains,fragment_chains,3.5)
        n_contacts.append(n_conts)
            
        name = path.split('/')[-3]
        fragment_name.append(name)
        (start,end) = name.split('_')[-1].split('-')
        fragment_start.append(int(start))
        fragment_end.append(int(end))
        fragment_center.append((fragment_start[-1]+fragment_end[-1])/2)

    conts_df = pd.DataFrame({
        'fragment_name':fragment_name,
        'fragment start (aa)':fragment_start,
        'fragment center (aa)':fragment_center,
        'fragment end (aa)':fragment_end,
        'n_contacts':n_contacts,
        'path':all_paths
    })
    conts_df = conts_df.sort_values(by='fragment start (aa)')
    return conts_df

def main(args):
    print('Loading JSON file specifying where colabfold results are located')
    # Load JSON file specifying where to import colabfold results from
    json_path = args.import_json
    with open(json_path,"r") as file:
        colab_results = json.loads(file.read())

    # For each gene and condition, import available data and create individual dataframes
    print('Processing results...')
    df_list = []
    merge_on_list = ['fragment_name','fragment start (aa)','fragment center (aa)','fragment end (aa)']
    for gene_name in colab_results:
        print('gene_name:',gene_name)
        for condition in colab_results[gene_name]:
            print('condition:',condition)

            # Load confidence (pLDDT/iPTM)
            confidence_df = get_confidence_dataframe(colab_results[gene_name][condition]['colabfold'])

            # Count contacts
            contacts_df = get_contact_dataframe(colab_results[gene_name][condition]['colabfold'])

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
    concat_df.to_csv("colabfold_predictions.csv")
    print(f"Dataframe with {len(concat_df)} entries total")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'ColabFoldProcessOutput',
        description = 'Script for processing the output of colabfold jobs and combining into a single dataframe')
    parser.add_argument('--import_json',required=True)
    parser.add_argument('--experimental_data',required=False)

    args = parser.parse_args()
    main(args)

    print('Done!')