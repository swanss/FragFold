import glob, json
from dataclasses import dataclass

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio.PDB import *

import glob, json, argparse, os, sys, re, functools
import string
import multiprocessing as mp

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from src.analyze_predictions import *

def plotContactRecovery(pred_df,name,dir_name,pos='',int_rec_ls='-'):
    
    plt.rcParams["figure.figsize"] = (7.5,3.75)
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    
    # weighted contacts gives a sense of how many peaks were predicted overall
    sns.lineplot(data=pred_df,x='fragment center (aa)',y='weighted contacts',ax=ax1)
    ax1.set_ylabel('ipTM-weighted contacts')
    ax1.set_xlabel('Fragment center (aa)')

    # recovery highlights how many of the fragments recover binding site residues/specific fragment-protein contacts
    sns.lineplot(data=pred_df,x='fragment center (aa)',y='binding site residue recovery',ax=ax2,color=sns.color_palette()[1])
    sns.lineplot(data=pred_df,x='fragment center (aa)',y='contact recovery',ax=ax2,color=sns.color_palette()[2],ls=int_rec_ls)
    ax2.set_xlabel('Fragment center (aa)')
    ax2.set_ylabel('Recovery (%)')
    ax2.set_ylim(0,1)
    
    legend_elements = [Patch(facecolor=sns.color_palette()[1], label='Binding site residue recovery'),
                       Patch(facecolor=sns.color_palette()[2], label='Interface contact recovery')]

    loc = 'upper right' if pos == '' else pos
    ax2.legend(handles=legend_elements, loc=loc)
    
    plt.savefig(f"{dir_name}/{name}.png",dpi=300)

def calculateContactRecoveryFromDFRow(native_contacts,contact_distance_cutoff,row):
    '''
    For each structure predicted by colabfold we will calculate the contact recovery
    '''
    # Load colabfold prediction, fix residue numbers, get contacts
    parser = PDBParser(QUIET=True)
    colabfold_structure = parser.get_structure("", row['path'])
    fixResidueNumbers(colabfold_structure[0][row['colabfold_fragment_chain']],row['fragment start (aa)'])
    prot_chain_set = set(row['colabfold_protein_chains'])
    frag_chain_set = set(row['colabfold_fragment_chain'])
    pred_contacts = getInterfaceContactsFromStructure(colabfold_structure,prot_chain_set,frag_chain_set,contact_distance_cutoff)
    cont_recovered = contactsRecovered(native_contacts,pred_contacts)
    bindingres_recovered = bindingSiteResiduesRecovered(native_contacts,pred_contacts)
    
    row_dict = row.to_dict()
    row_dict['n_contacts_recovered'] = cont_recovered[0]
    row_dict['frac_contacts_recovered'] = cont_recovered[0]/cont_recovered[1]
    row_dict['n_bindingres_recovered'] = bindingres_recovered[0]
    row_dict['frac_bindingres_recovered'] = bindingres_recovered[0]/bindingres_recovered[1]
    row_dict['all_contacts'] = [residueContactName(x,y) for x,y in pred_contacts]
    return row_dict


def main(args):
    print('Loading JSON file specifying where colabfold results are located')
    # Load JSON file specifying where to import colabfold results from
    json_path = args.import_json
    native_pdbs_dir_path = args.native_pdbs
    colabfold_data_csv_path  = args.colabfold_data_csv
    n_workers = len(os.sched_getaffinity(0))
    contact_distance_cutoff = args.contact_distance_cutoff

    print(f"Loading data with {n_workers} workers")

    with open(json_path,"r") as file:
        all_contact_comparison_data = json.loads(file.read())
    print(f"Loaded JSON file with values for {len(all_contact_comparison_data)} contact recovery comparisons")

    colabfold_data_df = pd.read_csv(colabfold_data_csv_path)
    print(f"Loaded colabfold data with {len(colabfold_data_df)} rows")

    dir_name = "contact_recovery_analysis"
    os.mkdir("dir_name")

    # For each gene and condition, import available data and create individual dataframes
    print('Calculating contact recovery...')
    df_list = []
    merge_on_list = ['fragment_name','rank','fragment start (aa)','fragment center (aa)','fragment end (aa)']
    for cont_comp_dict in all_contact_comparison_data:
        print(cont_comp_dict)

        # Load the native pdb, extract fragment + full-length protein chain(s), and find contacts
        start = cont_comp_dict['fragment_res_start']
        end = cont_comp_dict['fragment_res_start'] + cont_comp_dict['fragment_res_length'] - 1
        nat_path = os.path.join(native_pdbs_dir_path,cont_comp_dict['native_structure'])
        res_range_dict = {x:resRange() for x in cont_comp_dict['native_protein_chains']} #add the protein chain(s)
        res_range_dict[cont_comp_dict['native_fragment_chain']] = resRange(start,end) # add the fragment chain
        s_extract = extractFragmentFromNativeStructure(nat_path,res_range_dict)
        protein_chain_set = set(cont_comp_dict['native_protein_chains'])
        fragment_chain_set = set(cont_comp_dict['native_fragment_chain'])
        native_contacts_residues = getInterfaceContactsFromStructure(s_extract,protein_chain_set,fragment_chain_set,contact_distance_cutoff)

        # Get the rows corresponding to the colabfold gene/condition
        filt_df = colabfold_data_df[(colabfold_data_df['gene']==cont_comp_dict['colabfold_gene'])&
                                    (colabfold_data_df['condition']==cont_comp_dict['colabfold_condition'])].copy(deep=True)
       # Add all of the information from the json
        for key,val in cont_comp_dict:
            filt_df[key] = val

        # Calculate the contact recovery for each prediction
        with mp.Pool(n_workers) as pool:
            calculateContactRecoveryFromDFRow_mapper = functools.partial(calculateContactRecoveryFromDFRow,native_contacts_residues,contact_distance_cutoff)
            data = pool.map(func=calculateContactRecoveryFromDFRow_mapper,iterable=filt_df.iterrows(),chunksize=1)
        df = pd.DataFrame(data)
        print(f"{len(filt_df)} lines in the filtered dataframe and {len(df)} lines after calculating contact recovery")
        name = f"{cont_comp_dict['colabfold_gene']}_{cont_comp_dict['colabfold_condition']}"
        plotContactRecovery(df,name,dir_name)

        df_list.append(df)

    # Concatenate all dataframes into a single dataframe that stores all the data.
    concat_df = pd.concat(df_list,ignore_index=True)
    concat_df.to_csv("colabfold_contact_recovery.csv")
    print(f"Final dataframe with {len(concat_df)} entries total")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'ColabFoldContactRecovery',
        description = 'Script for processing the predicted structures from colabfold and calculating contact recovery. The JSON file controls the details of the analysis')
    parser.add_argument('--import_json',required=True)
    parser.add_argument('--colabfold_data_csv',required=True)
    parser.add_argument('--native_pdbs',required=True,
                        help='The path to the directory containing the native PDB structures')
    parser.add_argument('--contact_distance_cutoff',required=False,default=3.5,type=float,
                        help='The distance cutoff between heavy atoms of interface residues that defines a contact')

    args = parser.parse_args()
    main(args)

    print('Done!')