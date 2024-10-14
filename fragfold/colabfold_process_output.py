from pathlib import Path
from datetime import datetime
import glob, json, argparse, os, sys, re, functools
import string
import multiprocessing as mp

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio.PDB import PDBParser

from fragfold.src.colabfold_process_output_utils import *
from fragfold.src.analyze_predictions import *
from fragfold.src.plot_utils import *

def get_confidence_paths(dir_path):
    dir_path = Path(dir_path)
    assert dir_path.is_dir()
    all_paths = glob.glob(str(dir_path / "*/output/log.txt"))
    if len(all_paths) == 0:
        raise ValueError(f"Did not find any output in directory: {dir_path}")
    return all_paths

def load_confidence_data(path):
    # e.g. will match 2024-05-04 12:32:21,588 Query 1/1: ftsZ1copies_10-316_ftsZ_166-195 (length 337)
    fragment_name_pat = r"Query 1\/1: ((\D+)(\d)copies_(\d+)-(\d+)_(\D+)_(\d+)-(\d+)) \(length (\d+)\)"

    # e.g. will match 2023-05-26 20:03:23,691 rank_001_alphafold2_ptm_model_1_seed_000 pLDDT=93 pTM=0.801 ipTM=0.257 
    confidence_pat = r"rank_00(\d)_alphafold2_ptm_model_\d_seed_\d{3} pLDDT=([+-]?[0-9]*[.]?[0-9]+) pTM=([+-]?[0-9]*[.]?[0-9]+) ipTM=([+-]?[0-9]*[.]?[0-9]+)"
    
    confidence_dict = {'fragment_name':[],'rank':[],
                       'fragment_start_aa':[],'fragment_center_aa':[],'fragment_end_aa':[],
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

                confidence_dict['fragment_name'].append(name)
                confidence_dict['rank'].append(pred_rank)
                confidence_dict['fragment_start_aa'].append(start)
                confidence_dict['fragment_center_aa'].append(center)
                confidence_dict['fragment_end_aa'].append(end)
                confidence_dict['plddt'].append(plddt)
                confidence_dict['ptm'].append(ptm)
                confidence_dict['iptm'].append(iptm)
    if conf_match_count != 5:# or fragment_name_match_count != 1:
        raise ValueError(f"Expected to find 5 lines with confidence metrics and 1 with the fragment name, \
                         instead found ({conf_match_count},{fragment_name_match_count}) when searching {path}")

    return pd.DataFrame(confidence_dict)

def get_confidence_dataframe(all_paths,n_workers=1):

    if len(all_paths) == 0:
        raise ValueError(f"Unable to find colabfold output in the directory: {all_paths}")

    if n_workers > 1:
        with mp.Pool(n_workers) as pool:
            data = pool.map(func=load_confidence_data,iterable=all_paths,chunksize=1)
    else:
        data = list()
        for path in all_paths:
            data.append(load_confidence_data(path))

    confidence_df = pd.concat(data,ignore_index=True)
    confidence_df = confidence_df.sort_values(by='fragment_start_aa')
    print(f"{len(all_paths)} total paths and {len(confidence_df)} lines in the dataframe")
    return confidence_df

def get_pdb_paths(dir_path):
    dir_path = Path(dir_path)
    assert dir_path.is_dir()
    all_paths = glob.glob(str(dir_path / '*/output/*_unrelaxed_rank_00?_*.pdb'))
    if len(all_paths) == 0:
        raise ValueError(f"Did not find any output in directory: {dir_path}")
    return all_paths

def get_chains_from_structure(path):
    parser = PDBParser()
    structure = parser.get_structure(Path(path).stem, path)
    chain_id_list = [c.id for c in structure.get_chains()]
    return chain_id_list

def getRankFromPath(path):
    pat = r"rank_(00)?(\d)"
    match = re.search(pattern=pat,string=path)
    if match is None:
        raise ValueError
    else:
        return int(match.group(2))

def load_contact_data(path,protein_chains,fragment_chains,distance_cutoff):
    name = Path(path).stem.split('_unrelaxed')[0]
    pred_rank = getRankFromPath(path)
    
    start,end = name.split('_')[-1].split('-')
    start,end = int(start),int(end)
    center = (start + end) / 2

    parser = PDBParser(QUIET=True)
    s = parser.get_structure("s", path)
    fragment_chain = s[0][list(fragment_chains)[0]]
    fixResidueNumbers(fragment_chain,start)
    contacts = convertToString(getInterfaceContactsFromStructure(s,protein_chains,fragment_chains,distance_cutoff))
    n_conts = len(contacts)

    path = Path(path).resolve()
    return name,pred_rank,start,center,end,contacts,n_conts,path

def get_contact_dataframe(all_paths,contact_distance_cutoff,n_workers=1):

    if len(all_paths)==0:
        raise ValueError("No pdb files were provided")
        
    # Assumptions: last chain is the fragment, preceding chain(s) are the protein
    # Load the first structure to get names and exact lengths

    chain_list = get_chains_from_structure(all_paths[0])
    protein_chains = set(chain_list[:-1])
    fragment_chain = set(chain_list[-1])
    
    if n_workers > 1:
        # An alternative to starmap
        load_contact_data_mapper = functools.partial(load_contact_data,protein_chains=protein_chains,fragment_chains=fragment_chain,distance_cutoff=contact_distance_cutoff)
        with mp.Pool(n_workers) as pool:
            data = pool.map(func=load_contact_data_mapper,iterable=all_paths,chunksize=1)
    else:
        data = list()
        for path in all_paths:
            data.append(load_contact_data(path,protein_chains,fragment_chain,distance_cutoff=contact_distance_cutoff))

    conts_df = pd.DataFrame(data,columns=['fragment_name','rank','fragment_start_aa','fragment_center_aa','fragment_end_aa','contacts','n_contacts','path'])
    conts_df['protein_chains'] = ','.join(protein_chains)
    conts_df['fragment_chain'] = ','.join(fragment_chain)
    conts_df = conts_df.sort_values(by='fragment_start_aa')
    print(f"{len(all_paths)} total paths and {len(conts_df)} lines in the dataframe")
    return conts_df


def main(args):

    parse_v2_output = (args.predicted_pdbs is not None and args.confidence_logs is not None and 
                       args.full_protein is not None and args.fragment_protein is not None)
    parse_v1_output = (args.import_json is not None)
    if not parse_v2_output and not parse_v1_output:
        raise ValueError("Must provide --predicted_pdbs and --confidence_logs OR --import_json")

    n_workers = len(os.sched_getaffinity(0))
    contact_distance_cutoff = args.contact_distance_cutoff
    experimental_data_path = args.experimental_data

    print(f"Loading data with {n_workers} workers")

    print('Processing results...')

    if parse_v2_output:
        # Load confidence (pLDDT/iPTM)
        confidence_log_paths = args.confidence_logs
        # print(confidence_log_paths)
        confidence_df = get_confidence_dataframe(confidence_log_paths,n_workers)

        # Count contacts
        predicted_pdb_paths = args.predicted_pdbs
        # print(predicted_pdb_paths)
        contacts_df = get_contact_dataframe(predicted_pdb_paths,contact_distance_cutoff,n_workers)

        # Merge dataframes into a single one containing all types of data
        merge_on_list = ['fragment_name','rank','fragment_start_aa','fragment_center_aa','fragment_end_aa']
        comb_df = confidence_df.merge(contacts_df,on=merge_on_list,indicator=True,validate="one_to_one")
        print(f"{len(comb_df)} lines after merging")

        comb_df['fragment_parent_name'] =  args.fragment_protein
        comb_df['protein_name'] = args.full_protein
        comb_df['description'] = args.description
    else:
        json_path = Path(args.import_json)
        assert json_path.is_file()
        print(f"Loading JSON file specifying where colabfold results are located: {json_path}")
        with open(json_path,"r") as file:
            colab_results = json.loads(file.read())

        comb_df_list = []
        for fragfold_job_info in colab_results:
            dir_path = fragfold_job_info[0]

            # Load confidence (pLDDT/iPTM)
            confidence_log_paths = get_confidence_paths(dir_path)
            # print(confidence_log_paths)
            confidence_df = get_confidence_dataframe(confidence_log_paths,n_workers)
            # print(confidence_df.head())

            # Count contacts
            predicted_pdb_paths = get_pdb_paths(dir_path)
            # print(predicted_pdb_paths)
            contacts_df = get_contact_dataframe(predicted_pdb_paths,contact_distance_cutoff,n_workers)
            # print(contacts_df.head())

            # Merge dataframes into a single one containing all types of data
            merge_on_list = ['fragment_name','rank','fragment_start_aa','fragment_center_aa','fragment_end_aa']
            comb_df = confidence_df.merge(contacts_df,on=merge_on_list,indicator=True,validate="one_to_one")
            print(f"{len(comb_df)} lines after merging")

            comb_df['fragment_parent_name'] = fragfold_job_info[1]
            comb_df['protein_name'] = fragfold_job_info[2]
            comb_df['description'] = fragfold_job_info[3]
            comb_df_list.append(comb_df)
        
        comb_df = pd.concat(comb_df_list,ignore_index=True)

    # Calculate the weighted contacts
    comb_df['weighted_contacts'] = comb_df['n_contacts'] * comb_df['iptm']
    comb_df['fragment_length_aa'] = comb_df['fragment_end_aa'] - comb_df['fragment_start_aa'] + 1
    print(f"Combined dataframe with {len(comb_df)} rows")

    # Merge with experimental data df 


    # If experimental dataframe is provided, merge now
    if experimental_data_path != "":
        exp_df = pd.read_csv(experimental_data_path)
        print(f"Loaded experimental dataframe")
        merge_on_list = ['fragment_parent_name','fragment_start_aa','fragment_center_aa','fragment_end_aa','fragment_length_aa']
        # dataframe contains some replicated measurements (they are generally similar, so arbitrarily take the first)
        exp_df = exp_df.drop_duplicates(subset=merge_on_list)
        merge_df = comb_df.merge(exp_df,how='left',on=merge_on_list,validate='many_to_one')
        merge_df.to_csv("colabfold_predictions.csv")
    else:
        comb_df.to_csv("colabfold_predictions.csv")
        print(f"Dataframe with {len(comb_df)} entries total")

    if args.generate_plots:
        for groupers,group_df in comb_df.groupby(['fragment_parent_name','protein_name','fragment_length_aa','description']):
            ax = plotRawValuesOnSingle(group_df)
            plt.savefig(f"{datetime.today().strftime('%y%m%d')}_fragmentcenter_vs_weightedcontacts_combined_fragments-{'_'.join([str(x) for x in groupers])}.png",
                        dpi=300,bbox_inches='tight')
            g = plotRawValuesOnFacetGrid(group_df)
            plt.savefig(f"{datetime.today().strftime('%y%m%d')}_fragmentcenter_vs_weightedcontacts_facetgrid_fragments-{'_'.join([str(x) for x in groupers])}.png",
                        dpi=300,bbox_inches='tight')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'ColabFoldProcessOutput',
        description = 'Script for processing the output of colabfold jobs and combining into a single dataframe')
    parser.add_argument('--predicted_pdbs',nargs='+',required=False,
                        help='The paths to the PDBs predicted by ColabFold')
    parser.add_argument('--confidence_logs',nargs='+',required=False,
                        help='The paths to the log.txt files from the ColabFold predictions')
    parser.add_argument('--import_json',required=False,
                        help='A JSON file specifying the location of output from v1 of FragFold. If given, will ignore --predicted_pdbs and --confidence_logs')
    parser.add_argument('--full_protein',required=False)
    parser.add_argument('--fragment_protein',required=False)
    parser.add_argument('--description',required=False,default="")
    parser.add_argument('--experimental_data',required=False,default="")
    parser.add_argument('--contact_distance_cutoff',required=False,default=3.5,type=float,
                        help='The distance cutoff between heavy atoms of interface residues that defines a contact')
    parser.add_argument('--generate_plots',action='store_true',
                        help='If provided, will generate plots for each colabfold job')
    args = parser.parse_args()
    main(args)

    print('Done!')