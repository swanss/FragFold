import glob, json, copy, os, shutil
from dataclasses import dataclass

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib

from Bio.PDB import *

# The following functions are for analyzing generic structure predictions

def convertToString(contacts):
    return [residueContactName(rA,rB) for rA,rB in contacts]

## Heavy-atom contacts
def contactOverlap(contacts_a,contacts_b,min_overlap=1):
    if len(contacts_a) == 0 or len(contacts_b) == 0:
        return 0.0
    contacts_a,contacts_b = set(contacts_a),set(contacts_b)
    intersection = len(contacts_a.intersection(contacts_b))
    if intersection < min_overlap:
        return 0
    return intersection/len(contacts_a.union(contacts_b))

def contactsRecovered(contacts_a,contacts_b):
    contacts_a = convertToString(contacts_a)
    contacts_b = convertToString(contacts_b)
    # How contacts in a are recovered in b
    if len(contacts_a) == 0 or len(contacts_b) == 0:
        return (0,0)
    contacts_a,contacts_b = set(contacts_a),set(contacts_b)
    return len(contacts_a.intersection(contacts_b)),len(contacts_a)

def bindingSiteResiduesRecovered(contacts_a,contacts_b):
    # the first residue is considered to come from the binding site
    res_a = set([residueName(contact[0]) for contact in contacts_a])
    res_b = set([residueName(contact[0]) for contact in contacts_b])
    return len(res_a.intersection(res_b)),len(res_a)

def residueName(res):
    return f"{res.get_parent().id}_{res.id[1]}"

def residueContactName(resi,resj):
    if resi.get_parent().id <= resj.get_parent().id:
        return f"{residueName(resi)}-{residueName(resj)}"
    else:
        return f"{residueName(resj)}-{residueName(resi)}"

def isContact(chain_group_a,chain_group_b,res_1,res_2):
    if res_1.get_parent().id in chain_group_a and res_2.get_parent().id in chain_group_b:
            return True
    if res_1.get_parent().id in chain_group_b and res_2.get_parent().id in chain_group_a:
            return True
    return False

def isIntraChainContact(res_1,res_2,chain_id,min_distance_in_chain=5):
    if res_1.get_parent().id != chain_id or res_2.get_parent().id != chain_id:
        return False
    if abs(int(res_2.id[1]) - int(res_1.id[1])) < min_distance_in_chain:
        return False
    return True

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
    # sorted_contacts = [(x,y) if x.get_parent().id < y.get_parent().id else (y,x) for x,y in contacts]
    return contacts

def getIntraChainContactsFromStructure(structure, chain_id='A', contact_distance=4.0, min_distance_in_chain=5):
    ns = NeighborSearch([x for x in structure.get_atoms()])
    nearby_res = ns.search_all(contact_distance,'R')
    intra_chain_contacts = [(x,y) for x,y in nearby_res if isIntraChainContact(x,y,chain_id,min_distance_in_chain)]
    return intra_chain_contacts

def fixResidueNumbers(chain,start_res_num):
    res_list = [res for res in chain.get_residues()]
    current_start_res_num = res_list[0].id[1]
    offset = start_res_num - current_start_res_num
    # print(f"Offset is {offset}")
    if offset == 0:
        return
    for residue in res_list:
        res_info = list(residue.id)
        residue.detach_parent()
        res_info[1] = res_info[1] + offset
        residue.id = tuple(res_info)
    for residue in res_list:
        residue.set_parent(chain)
    

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
                           mobile,sup=None):
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
        
def getAtomsInRange(s1,s1_chain_id_list,s1_chain_res_range_map,s2,verbose=False):
    s1_resmap = createResidueMap([s1[0][chain_id] for chain_id in s1_chain_id_list])
    s2_resmap = createResidueMap([s2[0][chain.id] for chain in s2.get_chains()])
    # try to get the backbone atoms for all residues in the range, but if
    # either s1 or s2 doesn't have a residue, we skip it
    s1_atoms = []
    s2_atoms = []
    skip_count = 0
    for chain_id in s1_chain_id_list:
        chain  = s1[0][chain_id]
        for res_num in range(s1_chain_res_range_map[chain.id].first,s1_chain_res_range_map[chain.id].last):
            res_info = (chain_id,res_num)
            r1 = s1_resmap[res_info] if res_info in s1_resmap else None
            r2 = s2_resmap[res_info] if res_info in s2_resmap else None
            if r1 == None or r2 == None:
                skip_count+=1
                if verbose:
                    print(f"Residues are missing at residue number {res_num}. R1: {r1}. R2: {r2}")
                continue
            r1_bb_atoms = getBBAtoms(r1)
            r2_bb_atoms = getBBAtoms(r2)
            s1_atoms.extend(r1_bb_atoms)
            s2_atoms.extend(r2_bb_atoms)
    print(f"In total, missing {skip_count} residues s1 and or s2")
    assert len(s1_atoms) > 0,"When searching for residues, at least one matching residue should be identified (and likely many more)"
    return s1_atoms,s2_atoms

def createResidueMap(chain_list):
    resInfo_to_res = dict()
    for chain in chain_list:
        for res in chain.get_residues():
            res_info = (res.parent.id,res.id[1]) #(chainid,resnum)
            assert res_info not in resInfo_to_res
            resInfo_to_res[res_info] = res
    return resInfo_to_res

def createResRangeMap(chain_list):
    res_range_map = dict()
    for chain in chain_list:
        chain_residues = list(chain.get_residues())
        assert len(chain_residues) > 0
        first_residue_num = chain_residues[0].get_id()[1]
        last_residue_num = chain_residues[-1].get_id()[1]
        res_range_map[chain.id] = resRange(first_residue_num,last_residue_num)
    return res_range_map

def getBBAtoms(residue):
    bb_atoms = []
    all_atoms = {a.name:a for a in residue.get_atoms()}
    for atom_name in ['N','CA','C','O']:
        if atom_name in all_atoms:
            bb_atoms.append(all_atoms[atom_name])
        else:
            raise ValueError(f"Backbone atom {atom_name} not found in {residue}")
    return bb_atoms

def getBBAtomsFromResidueList(residue_list):
    all_bb_atoms = []
    for residue in residue_list:
        all_bb_atoms.extend(getBBAtoms(residue))
    return all_bb_atoms

def getBBAtomsFromStructure(structure):
    residue_list = []
    for chain in structure[0].get_chains():
        for residue in chain.get_residues():
            residue_list.append(residue)
    return getBBAtomsFromResidueList(residue_list)

def calcRMSD(atoms1,atoms2):
    atoms1 = np.array(atoms1)
    atoms2 = np.array(atoms2)
    return np.sqrt(np.mean(np.square(atoms1-atoms2)))

def calcRMSDFromRes(residues1,residues2):
    if (len(residues1) != len(residues2)):
        return np.nan
    atoms1 = getBBAtomsFromResidueList(residues1)
    atoms2 = getBBAtomsFromResidueList(residues2)
    return calcRMSD(atoms1,atoms2)

## Interface RMSD

def getInterfaceResidues(structure, chain_group_a, chain_group_b, contact_distance):
    interface_contacts = getInterfaceContactsFromStructure(structure, chain_group_a, chain_group_b, contact_distance)

    # group residues by parent chain ID
    chainid2resset = dict()
    for (res_a,res_b) in interface_contacts:
        chain_a = res_a.get_parent()
        if chain_a not in chainid2resset:
            chainid2resset[chain_a] = {res_a}
        else:
            chainid2resset[chain_a].add(res_a)
        
        chain_b = res_b.get_parent()
        if chain_b not in chainid2resset:
            chainid2resset[chain_b] = {res_b}
        else:
            chainid2resset[chain_b].add(res_b)
    
    sorted_chain_ids = sorted(list(chainid2resset.keys()))
    residue_list = []
    for chain_id in sorted_chain_ids:
        chain_residues = sorted(list(chainid2resset[chain_id]))
        residue_list.extend(chain_residues)
    return residue_list
        
# def getInterfaceResidues(res_list_1,res_list_2,contact_distance):
#     # Get interface residues
#     res_list_1.extend(res_list_2)
#     atom_list = [atom 
#                  for res in res_list_1
#                  for atom in res.get_atoms()]
#     ns = NeighborSearch(atom_list)
#     nearby_res = ns.search_all(contact_distance,'R')
#     interface_contacts = [(x,y) for x,y in nearby_res if x.get_parent()!=y.get_parent()]
#     return interface_contacts

def getIntersectionOfResidues(native_structure,colabfold_structure,native_selected_residues=[]):
    # NOTE: the residues that are returned are not necessarily the selected residues
    if len(native_selected_residues) == 0:
        for chain in native_structure[0].get_chains():
            for residue in chain.get_residues():
                native_selected_residues.append(residue)
    colabfoldresmap = createResidueMap(list(colabfold_structure[0].get_chains()))
    nativeresmap = createResidueMap(list(native_structure[0].get_chains()))
    native_intersection_res = []
    colabfold_intersection_res = []
    for residue in native_selected_residues:
        res_info = (residue.parent.id,residue.id[1])
        if res_info in colabfoldresmap and res_info in nativeresmap:
            native_intersection_res.append(nativeresmap[res_info])
            colabfold_intersection_res.append(colabfoldresmap[res_info])
    print(f"Of {len(native_selected_residues)} selected residues in the native structure, there are {len(native_intersection_res)} matching residues in the colabfold structure")
    return native_intersection_res,colabfold_intersection_res

def calculateInterfaceRMSD(native_structure,
                           protein_chains:list,
                           fragment_chains:list,
                           colabfold_structure,
                           contact_distance=8.0,
                           min_interface_res=10):
    # define the residues at the interface 
    native_contacts_residues = getInterfaceResidues(native_structure,set(protein_chains),set(fragment_chains),contact_distance)
    
    return calculateInterfaceRMSDFromSelectedRes(native_contacts_residues,native_structure,colabfold_structure,min_interface_res)

def verifyResiduesForInterfaceRMSDCalc(native_res_list,colabfold_res_list,minimum_res_per_chain=5):
    chain2rescount = {chain_id:0 for chain_id in set([res.get_parent().id for res in native_res_list])}
    for res in set(colabfold_res_list):
        chain_id = res.get_parent().id
        chain2rescount[chain_id] += 1
    for chain_id,count in chain2rescount.items():
        if count < minimum_res_per_chain:
            return False
    return True

def calculateInterfaceRMSDFromSelectedRes(native_interface_residues,
                                          native_structure,
                                          colabfold_structure,
                                          min_res_per_chain=5):
    # find matching residues in the colabfold structure
    filt_native_res,filt_colabfold_res = getIntersectionOfResidues(native_structure,colabfold_structure,native_interface_residues)

    if not verifyResiduesForInterfaceRMSDCalc(native_interface_residues,filt_colabfold_res,min_res_per_chain):
        return np.nan,None,None

    # optimally superimpose the two sets and report RMSD
    sup = Superimposer()
    sup.set_atoms(getBBAtomsFromResidueList(filt_native_res),getBBAtomsFromResidueList(filt_colabfold_res))
    colabfold_structure_copy = copy.deepcopy(colabfold_structure)
    sup.apply(colabfold_structure_copy.get_atoms())
    rmsd = sup.rms
    return rmsd,colabfold_structure_copy,filt_colabfold_res

class SelectFromSet(Select):
    def __init__(self,residue_list):
        self.sel_set = set()
        for residue in residue_list:
            self.sel_set.add((residue.parent.id,residue.id[1]))
        
    def accept_residue(self, residue):
        """Overload this to reject residues that do not match chain ID and residue number."""
        key = (residue.parent.id,residue.id[1])
        if key in self.sel_set:
            return 1
        return 0
    
def hasChainBreak(residue_list):
    prev_res = None
    for res in residue_list:
        if prev_res is not None:
            if areResDiscontiguous(prev_res,res):
                return True
    return False

def areResDiscontiguous(r1,r2):
    assert "C" in r1
    prev_C = r1["C"]
    assert "N" in r2
    current_N = r2["N"]
    distance = np.linalg.norm(prev_C.coord - current_N.coord)
    if distance > 1.5:
        return True
    return False

def findGapsFromResNum(res_list):
    assert len(res_list) > 0
    res_index = []
    start_res_num = res_list[0].get_id()[1]
    for i,res in enumerate(res_list):
        res_index_in_list = res.get_id()[1] - start_res_num
        res_list.append(res_index_in_list)
    return res_index


def createContactRecoveryDirectory(best_df,path_to_pml="/data1/groups/keatinglab/swans/savinovCollaboration/inhibitory_fragments_structure_prediction/inhib_frag_pred/pymol_scripts/interface_recovery.pml"):
    # create main dir
    dir_name = "native_contact_recovery"
    pathlib.Path(dir_name).mkdir(exist_ok=False)

    for i,row in best_df.iterrows():
        # for each native fragment, create a new subdir
        subdir_name = os.path.join(dir_name,row["fragment_name"])
        pathlib.Path(subdir_name).mkdir(exist_ok=False)

        # copy the native/colabfold
        new_native_structure = os.path.join(subdir_name,row['native_path'].split('/')[-1])
        shutil.copyfile(row['native_path'],new_native_structure)
        new_colabfold_structure = os.path.join(subdir_name,row['peptide_rmsd_structure_path'].split('/')[-1])
        shutil.copyfile(row['peptide_rmsd_structure_path'],new_colabfold_structure)

        # copy the script and replace values
        new_script = os.path.join(subdir_name,path_to_pml.split('/')[-1])
        shutil.copyfile(path_to_pml,new_script)

        with open(new_script,'r') as file:
            filedata = file.read()
        filedata = filedata.replace('NATIVE_STRUCTURE', f"'{row['native_path'].split('/')[-1]}'")
        filedata = filedata.replace('COLABFOLD_STRUCTURE', f"'{row['peptide_rmsd_structure_path'].split('/')[-1]}'")
        filedata = filedata.replace('CHAIN_ID', f"'{row['fragment_chain_id']}'")
        filedata = filedata.replace('NATIVE_CONTACTS', f"'{row['all_native_contacts']}'")
        filedata = filedata.replace('COLABFOLD_CONTACTS', f"'{row['all_contacts']}'")
        with open(new_script, 'w') as file:
            file.write(filedata)