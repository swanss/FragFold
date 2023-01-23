import glob, json

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio.PDB import *

def getChainInfoFromA3M(a3m_path):
    with open(a3m_path,'r') as file:
        chain_info_line = file.readline().rstrip()[1:] #strip trailing \n and preceding #
    n_residues,n_copies = chain_info_line.split('\t')
    n_residue_list = [int(x) for x in n_residues.split(',')]
    n_copies_list = [int(x) for x in n_copies.split(',')]
    assert len(n_residue_list) == len(n_copies_list)
    return n_copies_list,n_residue_list

## interface pTM

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

## Heavy-atom contacts

def isContact(chain_group_a,chain_group_b,res_1,res_2):
    if res_1.get_parent().id in chain_group_a and res_2.get_parent().id in chain_group_b:
            return True
    if res_1.get_parent().id in chain_group_b and res_2.get_parent().id in chain_group_a:
            return True
    return False

def countInterfaceContacts(structure_path, chain_group_a=set('A'), chain_group_b=set('B'), contact_distance = 4):
    parser = PDBParser()
    s = parser.get_structure("s", structure_path)
    ns = NeighborSearch([x for x in s.get_atoms()])
    nearby_res = ns.search_all(contact_distance,'R')
#     print(len(nearby_res),nearby_res[0:5])
    interface_contacts = [isContact(chain_group_a,chain_group_b,x,y) for x,y in nearby_res]
    return np.array(interface_contacts).sum()

## RMSD to native (or other) structures (after alignment by target)

def fixResidueNumbers(chain,start_res_num):
    # add an offset to residue numbers
    res_list = [res for res in chain.get_residues()]
    current_start_res_num = res_list[0].id[1]
    offset = start_res_num - current_start_res_num
    if offset == 0:
        return
    
    # set the direction of iteration depending on whether residue numbers will be 
    # shifted up or down
    if offset > 0:
        res_list.reverse()
    for res in res_list:
        res_info = list(res.id)
        res_info[1] = res_info[1] + offset
        res.id = tuple(res_info)

def alignStructureByChains(fixed,fixed_chain,mobile,mobile_chain,sup=None):
    if sup == None:
        sup = Superimposer()
    # find the optimal superposition between the selected chains
    #some of the residues may be missing, so walk through the structures manually
    numRes = len([x for x in fixed[0][fixed_chain].get_residues()])
    fixed_atoms,mobile_atoms = getAtomsInRange(fixed,fixed_chain,mobile,mobile_chain,11,numRes)
    
    sup.set_atoms(fixed_atoms,mobile_atoms)
    print('Aligned chains with RMSD:',sup.rms)
    
    # apply to the whole mobile structure
    sup.apply([a for a in mobile.get_atoms()])
    return sup.rms
        
def getAtomsInRange(s1,s1_chain_id,s2,s2_chain_id,start_res_num,n_res):
    s1_resmap = createResidueMap(s1[0][s1_chain_id])
    s2_resmap = createResidueMap(s2[0][s2_chain_id])
    # try to get the backbone atoms for all residues in the range, but if
    # either s1 or s2 doesn't have a residue, we skip it
    s1_atoms = []
    s2_atoms = []
    for res_num in range(start_res_num,start_res_num+n_res):
        r1 = s1_resmap[res_num] if res_num in s1_resmap else None
        r2 = s2_resmap[res_num] if res_num in s2_resmap else None
        if r1 == None or r2 == None:
#             print(f"Residues are missing at residue number {res_num}. R1: {r1}. R2: {r2}")
            continue
        r1_bb_atoms = getBBAtoms(r1)
        r2_bb_atoms = getBBAtoms(r2)
        s1_atoms.extend(r1_bb_atoms)
        s2_atoms.extend(r2_bb_atoms)
    return s1_atoms,s2_atoms

def createResidueMap(chain):
    resNum_to_res = dict()
    for res in chain.get_residues():
        res_number = res.id[1]
        assert res_number not in resNum_to_res
        resNum_to_res[res_number] = res
    return resNum_to_res

def getBBAtoms(residue):
    bb_atoms = []
    all_atoms = {a.name:a for a in residue.get_atoms()}
    for atom_name in ['N','CA','C','O']:
        if atom_name in all_atoms:
            bb_atoms.append(all_atoms[atom_name])
        else:
            raise ValueError(f"Backbone atom {atom_name} not found in {residue}")
    return bb_atoms

def calcRMSD(atoms1,atoms2):
    atoms1 = np.array(atoms1)
    atoms2 = np.array(atoms2)
    return np.sqrt(np.mean(np.square(atoms1-atoms2)))

## Interface RMSD

def getInterfaceContacts(res_list_1,res_list_2):
    res_list_1.extend(res_list_2)
    atom_list = [atom 
                 for res in res_list_1
                 for atom in res.get_atoms()]
    ns = NeighborSearch(atom_list)
    nearby_res = ns.search_all(contact_distance,'R')
    interface_contacts = [(x,y) for x,y in nearby_res if x.get_parent()!=y.get_parent()]
    return interface_contacts