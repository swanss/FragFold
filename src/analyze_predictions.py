import glob, json
from dataclasses import dataclass

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio.PDB import *

# The following functions are for analyzing generic structure predictions

## Heavy-atom contacts
def contactOverlap(contacts_a,contacts_b):
    if len(contacts_a) == 0 or len(contacts_b) == 0:
        return 0.0
    contacts_a,contacts_b = set(contacts_a),set(contacts_b)
    return len(contacts_a.intersection(contacts_b))/len(contacts_a.union(contacts_b))

def contactsRecovered(contacts_a,contacts_b):
    # How contacts in a are recovered in b
    if len(contacts_a) == 0 or len(contacts_b) == 0:
        return (0,0)
    contacts_a,contacts_b = set(contacts_a),set(contacts_b)
    return len(contacts_a.intersection(contacts_b)),len(contacts_a)

def bindingSiteResiduesRecovered(contacts_a,contacts_b):
    # the first residue is considered to come from the binding site
    res_a = set([contact[0] for contact in contacts_a])
    res_b = set([contact[0] for contact in contacts_b])
    return len(res_a.intersection(res_b)),len(res_a)

def residueContactName(resi,resj):
    if resi.get_parent().id <= resj.get_parent().id:
        return f"{resi.get_parent().id}_{resi.id[1]}-{resj.get_parent().id}_{resj.id[1]}"
    else:
        return f"{resj.get_parent().id}_{resj.id[1]}-{resi.get_parent().id}_{resi.id[1]}"

def isContact(chain_group_a,chain_group_b,res_1,res_2):
    if res_1.get_parent().id in chain_group_a and res_2.get_parent().id in chain_group_b:
            return True
    if res_1.get_parent().id in chain_group_b and res_2.get_parent().id in chain_group_a:
            return True
    return False

def countInterfaceContacts(structure_path, chain_group_a=set('A'), chain_group_b=set('B'), contact_distance = 4, structure = None):
    if structure == None:
        parser = PDBParser(QUIET=True)
        s = parser.get_structure("s", structure_path)
    else:
        s = structure
    contacts = getInterfaceContactsFromStructure(s,chain_group_a,chain_group_b,contact_distance)
    return len(contacts)

def getInterfaceContactsFromStructure(structure, chain_group_a=set('A'), chain_group_b=set('B'), contact_distance=4):
    ns = NeighborSearch([x for x in structure.get_atoms()])
    nearby_res = ns.search_all(contact_distance,'R')
    contacts = [(x,y) for x,y in nearby_res if isContact(chain_group_a,chain_group_b,x,y)]
    sorted_contacts = [(x,y) if x.get_parent().id < y.get_parent().id else (y,x) for x,y in contacts]
    return contacts

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

@dataclass
class resRange:
    first: float = -1
    last: float = -1

def extractFragmentFromNativeStructure(structure_path,chain_residue_selection):
    '''
    chain_residue_selection: dict
        key is chain ID, value is resRange
    '''

    parser = PDBParser(QUIET=True)
    s = parser.get_structure("s", structure_path)
    
    # we need to get the full list before iterating (as we may delete members during iteration)
    chain_list = [chain for chain in s.get_chains()]
    for chain in chain_list:
        # print(chain.id)
        if chain.id not in chain_residue_selection:
            # print(f"Deleting chain {chain.id}")
            s[0].detach_child(chain.id)
            continue
        res_sel = chain_residue_selection[chain.id]
        if res_sel.first == -1 or res_sel.last == -1:
            # take all the residues
            continue
        # print('range:',res_sel.first,res_sel.last)
        residue_list = [residue for residue in chain.get_residues()]
        for residue in residue_list:
            # print(residue.id)
            res_number = int(residue.id[1])
            if res_number >= res_sel.first and res_number <= res_sel.last:
                # keep residue
                pass
            else:
                # residue falls out of range, detach
                # print(f"Deleting residue {residue.id} from chain {chain.id}")
                s[0][chain.id].detach_child(residue.id)
    return s

## RMSD to native (or other) structures (after alignment by target)

def alignStructureByChains(fixed,fixed_chain_ids_list,fixed_chains_res_range_map,
                           mobile,
                           sup=None):
    if sup == None:
        sup = Superimposer()
    # find the optimal superposition between the selected chains
    #some of the residues may be missing, so walk through the structures manually
    numRes = len([x for fixed_chain_id in fixed_chain_ids_list for x in fixed[0][fixed_chain_id].get_residues()])
    fixed_atoms,mobile_atoms = getAtomsInRange(fixed,fixed_chain_ids_list,fixed_chains_res_range_map,mobile)
    
    sup.set_atoms(fixed_atoms,mobile_atoms)
    print('Aligned chains with RMSD:',sup.rms)
    
    # apply to the whole mobile structure
    sup.apply([a for a in mobile.get_atoms()])
    return sup.rms
        
def getAtomsInRange(s1,s1_chain_id_list,s1_chain_res_range_map,s2):
    s1_resmap = createResidueMap([s1[0][chain_id] for chain_id in s1_chain_id_list])
    s2_resmap = createResidueMap([s2[0][chain.id] for chain in s2.get_chains()])
    # try to get the backbone atoms for all residues in the range, but if
    # either s1 or s2 doesn't have a residue, we skip it
    s1_atoms = []
    s2_atoms = []
    for chain_id in s1_chain_id_list:
        chain  = s1[0][chain_id]
        for res_num in range(s1_chain_res_range_map[chain.id].first,s1_chain_res_range_map[chain.id].last):
            res_info = (chain_id,res_num)
            r1 = s1_resmap[res_info] if res_info in s1_resmap else None
            r2 = s2_resmap[res_info] if res_info in s2_resmap else None
            if r1 == None or r2 == None:
                print(f"Residues are missing at residue number {res_num}. R1: {r1}. R2: {r2}")
                continue
            r1_bb_atoms = getBBAtoms(r1)
            r2_bb_atoms = getBBAtoms(r2)
            s1_atoms.extend(r1_bb_atoms)
            s2_atoms.extend(r2_bb_atoms)
    return s1_atoms,s2_atoms

def createResidueMap(chain_list):
    resInfo_to_res = dict()
    for chain in chain_list:
        for res in chain.get_residues():
            res_info = (res.parent.id,res.id[1]) #(chainid,resnum)
            assert res_info not in resInfo_to_res
            resInfo_to_res[res_info] = res
    return resInfo_to_res

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
