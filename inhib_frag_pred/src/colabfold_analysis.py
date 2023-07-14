import glob, json

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio.PDB import *

def createiPTMDF(all_pae_paths,prot_len,frag_len=30):
    iptm_list = []
    plddt_list = []
    fragment_name = []
    fragment_start = []
    fragment_center = []
    fragment_end = []

    # asym_id = np.array([0]*307+[1]*30)
    asym_id = np.array([0]*prot_len+[1]*frag_len)
    for path in all_pae_paths:
        with open(path,'r') as file:
            data = json.load(file)
            
        plddt = float(np.mean(np.array(data['plddt'][-frag_len:])))
        plddt_list.append(plddt)
            
        iptm = float(predicted_tm_score(np.array(data['pae']),asym_id,True))
        iptm_list.append(iptm)
        
        name = path.split('/')[-3]
        fragment_name.append(name)
        (start,end) = name.split('_')[-1].split('-')
        fragment_start.append(int(start))
        fragment_end.append(int(end))
        fragment_center.append((fragment_start[-1]+fragment_end[-1])/2)

    iptm_df = pd.DataFrame({
        'fragment_name':fragment_name,
        'start':fragment_start,
        'aa.fragmentCenter':fragment_center,
        'end':fragment_end,
        'iptm':iptm_list,
        'plddt':plddt_list
    })
    iptm_df = iptm_df.sort_values(by='start')
    return iptm_df

def getContDF(all_pdb_paths,contact_distance=3.5):
    fragment_name = []
    fragment_start = []
    fragment_end = []
    fragment_center = []
    n_contacts = []

    parser = PDBParser()
    for path in all_pdb_paths:
        s = parser.get_structure("", path)
        n_conts = countInterfaceContacts(s,contact_distance)
        n_contacts.append(n_conts)
            
        name = path.split('/')[-3]
        fragment_name.append(name)
        (start,end) = name.split('_')[-1].split('-')
        fragment_start.append(int(start))
        fragment_end.append(int(end))
        fragment_center.append((fragment_start[-1]+fragment_end[-1])/2)

    conts_df = pd.DataFrame({
        'fragment_name':fragment_name,
        'start':fragment_start,
        'aa.fragmentCenter':fragment_center,
        'end':fragment_end,
        'n_contacts':n_contacts,
        'path':all_pdb_paths
    })
    conts_df = conts_df.sort_values(by='start')
    return conts_df

# Function originally from AlphaFold, modified to work with different input
def predicted_tm_score(
    predicted_aligned_error: np.ndarray,
    asym_id:[np.ndarray] = None,
    interface: bool = False) -> np.ndarray:
    """Computes predicted TM alignment or predicted interface TM alignment score.

    Args:
    predicted_aligned_error: [num_res, num_res] the expected aligned distance
      error for each pair of residues.
    asym_id: [num_res] the asymmetric unit ID - the chain ID. Only needed for
      ipTM calculation, i.e. when interface=True.
    interface: If True, interface predicted TM score is computed.

    Returns:
    ptm_score: The predicted TM alignment or the predicted iTM score.
    """
    
    num_res = int(predicted_aligned_error.shape[0])
    # Clip num_res to avoid negative/undefined d0.
    clipped_num_res = max(num_res, 19)

    # Compute d_0(num_res) as defined by TM-score, eqn. (5) in Yang & Skolnick
    # "Scoring function for automated assessment of protein structure template
    # quality", 2004: http://zhanglab.ccmb.med.umich.edu/papers/2004_3.pdf
    d0 = 1.24 * (clipped_num_res - 15) ** (1./3) - 1.8

    # TM-Score term for every bin.
    predicted_tm = 1. / (1 + np.square(predicted_aligned_error) / np.square(d0))

    pair_mask = np.ones(shape=(num_res, num_res), dtype=bool)
    if interface:
        pair_mask *= asym_id[:, None] != asym_id[None, :]

    predicted_tm *= pair_mask
    
    residue_weights = np.ones(predicted_aligned_error.shape[0])
    pair_residue_weights = pair_mask * (
      residue_weights[None, :] * residue_weights[:, None])
    normed_residue_mask = pair_residue_weights / (1e-8 + np.sum(
      pair_residue_weights, axis=-1, keepdims=True))
    per_alignment = np.sum(predicted_tm * normed_residue_mask, axis=-1)
    return np.asarray(per_alignment[(per_alignment * residue_weights).argmax()])


def countInterfaceContacts(s, contact_distance = 4):
    s = parser.get_structure("test", path)
    ns = NeighborSearch([x for x in s.get_atoms()])
    nearby_res = ns.search_all(contact_distance,'R')
#     print(len(nearby_res),nearby_res[0:5])
    interface_contacts = [x.get_parent()!=y.get_parent() for x,y in nearby_res]
    return np.array(interface_contacts).sum()

def rescaleSeriesAtoB(df_series_a,df_series_b,a_min=None,b_min=None):
    max_val_a = df_series_a.max()
    if a_min is not None:
        min_val_a = a_min
    else:
        min_val_a = df_series_a.min()
    max_val_b = df_series_b.max()
    if b_min is not None:
        min_val_b = b_min
    else:
        min_val_b = df_series_b.min()
    print(max_val_a,min_val_a)
    print(max_val_b,min_val_b)
    return ((max_val_b-min_val_b)*((df_series_a - min_val_a)/(max_val_a-min_val_a)))+min_val_b