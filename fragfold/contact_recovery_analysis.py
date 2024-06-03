import glob, json, pathlib, copy
from dataclasses import dataclass
from shutil import rmtree

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import time

from Bio.PDB import *

import glob, json, argparse, os, sys, re, functools
import string
import multiprocessing as mp

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from fragfold.src.analyze_predictions import *

# Create global variables for PDB parsing/writing
io = PDBIO()
pdbparser = PDBParser(QUIET=True)

def getJobName(contact_comp_dict):
    return f"{contact_comp_dict['colabfold_gene']}_{contact_comp_dict['colabfold_condition']}_{contact_comp_dict['fragment_res_start']}_{contact_comp_dict['fragment_res_length']}"

def plotContactRecovery(pred_df,name,dir_name,pos='',int_rec_ls='-'):
    
    plt.rcParams["figure.figsize"] = (7.5,3.75)
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    
    # weighted contacts gives a sense of how many peaks were predicted overall
    sns.lineplot(data=pred_df,x='fragment_center_aa',y='weighted_contacts',ax=ax1)
    ax1.set_ylabel('ipTM-weighted contacts')
    ax1.set_xlabel('fragment_center_aa')

    # recovery highlights how many of the fragments recover binding site residues/specific fragment-protein contacts
    sns.lineplot(data=pred_df,x='fragment_center_aa',y='frac_bindingres_recovered',ax=ax2,color=sns.color_palette()[1])
    sns.lineplot(data=pred_df,x='fragment_center_aa',y='frac_contacts_recovered',ax=ax2,color=sns.color_palette()[2],ls=int_rec_ls)
    ax2.set_xlabel('fragment_center_aa')
    ax2.set_ylabel('Recovery (%)')
    ax2.set_ylim(0,1)
    
    legend_elements = [Patch(facecolor=sns.color_palette()[1], label='Binding site residue recovery'),
                       Patch(facecolor=sns.color_palette()[2], label='Interface contact recovery')]

    loc = 'upper right' if pos == '' else pos
    ax2.legend(handles=legend_elements, loc=loc)
    
    plt.savefig(f"{dir_name}/{name}.png",dpi=300)

def calculateStructureRecoveryFromDFRow(native_contacts,contact_distance_cutoff,lowres_interface_residues,contact_comparison_data,index_value_tuple):
    '''
    For each structure predicted by colabfold we will calculate the contact recovery, ligand RMSD, and interface RMSD
    '''
    print(f"row: {index_value_tuple[0],index_value_tuple[1]}")
    row = index_value_tuple[1]
    native_structure = lowres_interface_residues[0].get_parent().get_parent().get_parent()
    job_name = "native_"+getJobName(contact_comparison_data)

    ### Load colabfold prediction, fix residue numbers, get contacts
    colabfold_structure = pdbparser.get_structure("", row['path'])
    
    colab2native_chainid = {colab_chain_id:native_chain_id for colab_chain_id,native_chain_id in zip(contact_comparison_data['colabfold_protein_chains'],contact_comparison_data['native_protein_chains'])}
    colab_chainid2resstart = {chain_id:res_start for chain_id,res_start in zip(contact_comparison_data['colabfold_protein_chains'],contact_comparison_data['colabfold_protein_chains_start'])}
    
    colabfold_chains = []
    for chain in list(colabfold_structure[0].get_chains()):
        colabfold_structure[0].detach_child(chain.id)
        if (chain.id == contact_comparison_data['colabfold_fragment_chain']):
            # print(f"Adjust fragment chain, start: {row['fragment_id']}")
            fixResidueNumbers(chain,row['fragment_id'])
            chain.id = contact_comparison_data['native_fragment_chain']
        else:
            # print(f"Adjust protein chain, start: {colab_chainid2resstart[chain.id]}")
            fixResidueNumbers(chain,colab_chainid2resstart[chain.id])
            chain.id = colab2native_chainid[chain.id]
        colabfold_chains.append(chain)
    for chain in colabfold_chains:
        colabfold_structure[0].add(chain)

    ### Calculate contact recovery
    prot_chain_set = set(list(contact_comparison_data['native_protein_chains']))
    frag_chain_set = set(list(contact_comparison_data['native_fragment_chain']))
    pred_contacts = getInterfaceContactsFromStructure(colabfold_structure,prot_chain_set,frag_chain_set,contact_distance_cutoff)
    cont_recovered = contactsRecovered(native_contacts,pred_contacts)
    bindingres_recovered = bindingSiteResiduesRecovered(native_contacts,pred_contacts)
    frac_bindingres_recovered = bindingres_recovered[0]/bindingres_recovered[1] if bindingres_recovered[1] > 0 else 0

    ### Calculate ligand RMSD
    fixed_chains_res_range_map = createResRangeMap(list(native_structure[0].get_chains()))
    target_rmsd = alignStructureByChains(native_structure,contact_comparison_data['native_protein_chains'],
                           fixed_chains_res_range_map,colabfold_structure)
    native_peptide_res = list(native_structure[0][contact_comparison_data['native_fragment_chain']].get_residues())
    colab_peptide_res = list(colabfold_structure[0][contact_comparison_data['native_fragment_chain']].get_residues())
    if (len(native_peptide_res) == len(colab_peptide_res)):
        peptide_rmsd = calcRMSDFromRes(native_peptide_res,colab_peptide_res)
    else:
        if (hasChainBreak(native_peptide_res)):
            filt_colab_peptide_res = [colab_peptide_res[idx] for idx in findGapsFromResNum(native_peptide_res)]
            peptide_rmsd = calcRMSDFromRes(native_peptide_res,filt_colab_peptide_res)
        else:
            peptide_rmsd = calcRMSDFromRes(native_peptide_res,colab_peptide_res[:len(native_peptide_res)])
    print(f"peptide RMSD: {peptide_rmsd}")
    
    colabfold_name = f"colabfold_{contact_comparison_data['colabfold_gene']}_{row['fragment_id']}_{row['rank']}"
    peptide_rmsd_path = os.path.join("contact_recovery_analysis",job_name,colabfold_name+"_ligandrmsd.pdb")
    if peptide_rmsd < 10 or frac_bindingres_recovered > 0.3:
        io.set_structure(colabfold_structure)
        io.save(peptide_rmsd_path)

    ## Calculate interface RMSD
    native_structure_copy = copy.deepcopy(native_structure)
    interface_rmsd,realigned_colabfold_structure,filt_colabfold_res = calculateInterfaceRMSDFromSelectedRes(lowres_interface_residues,native_structure_copy,colabfold_structure)
    print(f"interface RMSD: {interface_rmsd}")
    
    interface_rmsd_path = ""
    if interface_rmsd < 5.0:
        resSelector = SelectFromSet(filt_colabfold_res)
        io.set_structure(realigned_colabfold_structure)
        interface_rmsd_path = os.path.join("contact_recovery_analysis",job_name,colabfold_name+"_interfacermsd.pdb")
        io.save(interface_rmsd_path,resSelector)

    row_dict = row.to_dict()
    row_dict['n_contacts_recovered'] = cont_recovered[0]
    row_dict['frac_contacts_recovered'] = cont_recovered[0]/cont_recovered[1] if cont_recovered[1] > 0 else 0
    row_dict['n_bindingres_recovered'] = bindingres_recovered[0]
    row_dict['frac_bindingres_recovered'] = frac_bindingres_recovered
    row_dict['all_contacts'] = ','.join([residueContactName(x,y) for x,y in pred_contacts])
    row_dict['n_contacts'] = len(pred_contacts)
    row_dict['target_rmsd'] = target_rmsd
    row_dict['peptide_rmsd'] = peptide_rmsd
    row_dict['peptide_rmsd_structure_path'] = os.path.join(os.getcwd(),peptide_rmsd_path)
    row_dict['n_peptide_res'] = len(native_peptide_res)
    row_dict['interface_rmsd'] = interface_rmsd
    row_dict['interface_rmsd_structure_path'] = os.path.join(os.getcwd(),interface_rmsd_path)
    row_dict['n_interface_res'] = len(filt_colabfold_res) if filt_colabfold_res is not None else 0
    return row_dict

def main(args):
    print('Loading JSON file specifying where colabfold results are located')
    # Load JSON file specifying where to import colabfold results from
    json_path = args.import_json
    native_pdbs_dir_path = args.native_pdbs
    colabfold_data_csv_path  = args.colabfold_data_csv
    n_workers = min(48,len(os.sched_getaffinity(0)))
    # n_workers = 1
    contact_distance_cutoff = args.contact_distance_cutoff

    print(f"Loading data with {n_workers} workers")

    with open(json_path,"r") as file:
        all_contact_comparison_data = json.loads(file.read())
    print(f"Loaded JSON file containing {len(all_contact_comparison_data)} native interface structures to compare to...")

    colabfold_data_df = pd.read_csv(colabfold_data_csv_path,index_col=0)
    print(f"Loaded colabfold data with {len(colabfold_data_df)} total predicted structures and {colabfold_data_df.fragment_name.nunique()} unique fragments")
    print("Columns:")
    print(colabfold_data_df.columns)

    intermediate_df_dir = "intermediate_data"
    pathlib.Path(intermediate_df_dir).mkdir(exist_ok=True)

    dir_name = "contact_recovery_analysis"
    pathlib.Path(dir_name).mkdir(exist_ok=True)

    # For each gene and condition, import available data and create individual dataframes
    print('Calculating contact recovery...')
    df_list = []
    # merge_on_list = ['fragment_name','rank','fragment_id','fragment_center_aa','fragment_end_aa']
    for cont_comp_dict in all_contact_comparison_data:
        print(cont_comp_dict)
        start = time.time()

        job_name = "native_"+getJobName(cont_comp_dict)
        subdir_name = os.path.join(dir_name,job_name)
        done_file = os.path.join(subdir_name,"done.txt")
        if os.path.isdir(subdir_name) and os.path.isfile(done_file):
            print(f"Skipping {subdir_name}, as the jobs are already complete")
            df = pd.read_csv(os.path.join(intermediate_df_dir,f"{job_name}.csv"),index_col=0)
        else:
            try:
                pathlib.Path(subdir_name).mkdir(exist_ok=False)
            except OSError:
                # If the directory exists, but the jobs didn't complete we'll just delete and rerun
                rmtree(subdir_name)
                pathlib.Path(subdir_name).mkdir(exist_ok=False)

            # Load the native pdb, extract fragment + full-length protein chain(s), and find atomic contacts
            start = cont_comp_dict['fragment_res_start']
            end = cont_comp_dict['fragment_res_start'] + cont_comp_dict['fragment_res_length'] - 1
            nat_path = os.path.join(native_pdbs_dir_path,cont_comp_dict['native_structure'])
            res_range_dict = {x:resRange() for x in list(cont_comp_dict['native_protein_chains'])} #add the protein chain(s)
            res_range_dict[cont_comp_dict['native_fragment_chain']] = resRange(start,end) # add the fragment chain
            s_extract = extractFragmentFromNativeStructure(nat_path,res_range_dict)
            io.set_structure(s_extract)
            io.save(job_name+".pdb")

            protein_chain_set = set(list(cont_comp_dict['native_protein_chains']))
            fragment_chain_set = set(list(cont_comp_dict['native_fragment_chain']))
            hires_native_contacts_residues = getInterfaceContactsFromStructure(s_extract,protein_chain_set,fragment_chain_set,contact_distance_cutoff)
            hires_native_contacts_str = ','.join([residueContactName(x,y) for x,y in hires_native_contacts_residues])
            print(f"{len(hires_native_contacts_residues)} high-resolution contacts ({contact_distance_cutoff} Å): {hires_native_contacts_str}")

            # Define a second of set of interface contacts with a higher cutoff
            lores_interface_residues = getInterfaceResidues(s_extract,protein_chain_set, fragment_chain_set, contact_distance=8.0)
            print(f"{len(lores_interface_residues)} interface residues defined from low-resolution (8 Å)")

            if len(hires_native_contacts_residues) == 0 or len(lores_interface_residues) == 0:
                raise ValueError("At least some residues must be defined at the interface")
                
            resSelector = SelectFromSet(lores_interface_residues)
            io.save(job_name+"_8Ainterfaceresidues.pdb",resSelector)

            # Get the rows corresponding to the colabfold gene/condition
            filt_df = colabfold_data_df[(colabfold_data_df['gene']==cont_comp_dict['colabfold_gene'])&
                                        (colabfold_data_df['condition']==cont_comp_dict['colabfold_condition'])].copy(deep=True).reset_index(drop=True)
            if len(filt_df) <= 0:
                raise ValueError(f"Dataframe does not contain values matching gene: {cont_comp_dict['colabfold_gene']} and condition: {cont_comp_dict['colabfold_condition']}")

            # Calculate the contact recovery for each prediction
            if n_workers == 1:
                data = []
                for i,row in filt_df.iterrows():
                    data.append(calculateStructureRecoveryFromDFRow(hires_native_contacts_residues,contact_distance_cutoff,lores_interface_residues,cont_comp_dict,(i,row)))
            else:
                with mp.Pool(n_workers) as pool:
                    calculateContactRecoveryFromDFRow_mapper = functools.partial(calculateStructureRecoveryFromDFRow,hires_native_contacts_residues,contact_distance_cutoff,lores_interface_residues,cont_comp_dict)
                    data = pool.map(func=calculateContactRecoveryFromDFRow_mapper,iterable=filt_df.iterrows(),chunksize=1)

            df = pd.DataFrame(data)

            print(f"{len(filt_df)} lines in the filtered dataframe and {len(df)} lines after calculating contact recovery")
            name = f"{cont_comp_dict['colabfold_gene']}_{cont_comp_dict['colabfold_condition']}"
            df['gene'] = cont_comp_dict['colabfold_gene']
            df['condition'] = cont_comp_dict['colabfold_condition']
            df['native_fragment'] = job_name
            df['native_path'] = os.path.join(os.getcwd(),job_name+".pdb")
            df['all_native_contacts'] = hires_native_contacts_str
            df['fragment_chain_id'] = cont_comp_dict['native_fragment_chain']
            plotContactRecovery(df,name,dir_name)

            df.to_csv(os.path.join(intermediate_df_dir,f"{job_name}.csv"))

        df_list.append(df)

        with open(done_file,"w") as file:
            # If all the jobs completed without error, make a record of it
            pass

        stop = time.time()
        print(f"Elapsed: {stop - start} s")

    # Concatenate all dataframes into a single dataframe that stores all the data.
    concat_df = pd.concat(df_list,ignore_index=True)
    concat_df.to_csv("colabfold_contact_recovery.csv")
    print(f"Final dataframe with {len(concat_df)} entries total")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'ColabFoldContactRecovery',
        description = 'Script for processing the predicted structures from colabfold and calculating contact recovery. The JSON file controls the details of the analysis')
    parser.add_argument('--import_json',required=True,type=str)
    # parser.add_argument('--overwrite',default=True,action='store_true')
    parser.add_argument('--colabfold_data_csv',required=True,type=str)
    parser.add_argument('--native_pdbs',required=True,type=str,
                        help='The path to the directory containing the native PDB structures')
    parser.add_argument('--contact_distance_cutoff',required=False,default=3.5,type=float,
                        help='The distance cutoff between heavy atoms of interface residues that defines a contact')

    args = parser.parse_args()
    print(args)
    main(args)
    print('Done!')