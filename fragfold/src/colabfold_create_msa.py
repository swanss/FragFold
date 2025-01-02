import math
import numpy as np
from pathlib import Path
import os
import random
from multiprocessing import Pool
from fragfold.src import a3m_tools
import copy


def calcHammingDistance(s1,s2):
    assert len(s1) == len(s2)
    return np.sum([a!=b for a,b in zip(s1,s2)])

def shuffleSeq(seq,maxPercentIdentity=0.3):
    assert len(seq) > 1
    minHammingDistance = math.floor((1-maxPercentIdentity)*len(seq))
    hammingDistance = None
    shuffled_seq = None
    while hammingDistance is None or hammingDistance < minHammingDistance:
        seq_list = list(seq)
        random.shuffle(seq_list)
        shuffled_seq = ''.join(seq_list)
        hammingDistance = calcHammingDistance(shuffled_seq,seq)
        # print(shuffled_seq,hammingDistance)
    return shuffled_seq


def subsampleMSA(msa: a3m_tools.MSAa3m, subsample: int):
    # subsample the MSA
    if subsample > 0:
        print("Subsampling MSA")
        msa.sequences = msa.sequences[:subsample]
    return msa


def createMSA(
    inputMSA: str, 
    protein_range: tuple, 
    fragment_range: tuple, 
    a3m_output_path: str,
    subsample: int = -1,
    protein_copies: int = 1, 
    fragment_single_sequence=False, 
    fragment_shuffle=False
):
    msa = a3m_tools.MSAa3m.from_a3m_file(inputMSA)
    # assuming that `protein_range` and `fragment_range` are 1-based
    protein_msa = msa[protein_range[0]-1:protein_range[1]]
    fragment_msa = msa[fragment_range[0]-1:fragment_range[1]]
    if fragment_shuffle:
        # get rid of all but the query sequence in the fragment MSA
        # shuffle the query sequence
        shuffled_query = a3m_tools.ProteinSequence(
            fragment_msa.query.header,
            shuffleSeq(fragment_msa.query.seq_str)
        )
        fragment_msa = a3m_tools.MSAa3m(
            fragment_msa.info_line, 
            shuffled_query, 
            []
        )
    if fragment_single_sequence:
        # get rid of all but the query sequence in the fragment MSA
        fragment_msa = a3m_tools.MSAa3m(
            fragment_msa.info_line, 
            fragment_msa.query, 
            []
        )
    protein_msa = subsampleMSA(protein_msa,subsample)
    if protein_copies > 1:
        protein_msa.set_cardinality(protein_copies)
    fragment_msa = subsampleMSA(fragment_msa,subsample)
    combined_msa = protein_msa + fragment_msa
    combined_msa.save(a3m_output_path)


def createMSAHeteromicInteraction(
    fullLengthProteinMSA: str, 
    protein_range: tuple, 
    fragmentProteinMSA: str, 
    fragment_range: tuple, 
    a3m_output_path: str, 
    subsample: int = -1,
    protein_copies: int = 1,
    fragmentSingleSequence=False, 
    fragmentShuffle=False
):
    protein_msa = a3m_tools.MSAa3m.from_a3m_file(fullLengthProteinMSA)
    # assuming that `protein_range` and `fragment_range` are 1-based
    protein_msa = protein_msa[protein_range[0]-1:protein_range[1]]

    fragment_msa = a3m_tools.MSAa3m.from_a3m_file(fragmentProteinMSA)
    fragment_msa = fragment_msa[fragment_range[0]-1:fragment_range[1]]
    if fragmentShuffle:
        # get rid of all but the query sequence in the fragment MSA
        # shuffle the query sequence
        shuffled_query = a3m_tools.ProteinSequence(
            fragment_msa.query.header,
            shuffleSeq(fragment_msa.query.seq_str)
        )
        fragment_msa = a3m_tools.MSAa3m(
            fragment_msa.info_line, 
            shuffled_query, 
            []
        )
    if fragmentSingleSequence:
        # get rid of all but the query sequence in the fragment MSA
        fragment_msa = a3m_tools.MSAa3m(
            fragment_msa.info_line, 
            fragment_msa.query, 
            []
        )
    protein_msa = subsampleMSA(protein_msa,subsample)
    if protein_copies > 1:
        protein_msa.set_cardinality(protein_copies)
    fragment_msa = subsampleMSA(fragment_msa,subsample)
    combined_msa = protein_msa + fragment_msa
    combined_msa.save(a3m_output_path)

def readA3MProteinLength(a3m_path):
    with open(a3m_path,'r') as file:
        first_line = file.readline()
    if first_line[0] != "#":
        raise ValueError(f"the first line of .a3m file: {a3m_path} should begin with '#'")
    nRes = int(first_line[1:].split('\t')[0])
    return nRes

def replaceNullProteinRange(protein_length,start,stop):
    # If necessary, define ranges from the provided sequences
    start = start if start != -1 else 1
    stop = stop if stop != -1 else protein_length
    return (start,stop)

def replaceNullFragmentRange(protein_length,fragment_length,start,stop):
    start = start if start != -1 else 1
    stop = stop if stop != -1 else protein_length - fragment_length + 1
    return (start,stop)

def verifyFragmentNterminalResRange(fragment_start_range,fragment_length,protein_n_res):
    if (len(fragment_start_range) != 2) or \
       (fragment_start_range[1] < fragment_start_range[0]) or \
       (fragment_start_range[0] < 1) or \
       (fragment_start_range[1] > protein_n_res - fragment_length + 1):
        raise ValueError(f"Provided fragment n-terminal residue range: ({fragment_start_range[0]},{fragment_start_range[1]}) with fragment length: {fragment_length} and a total protein length of {protein_n_res} is invalid")
    return True

def verifyProteinRange(protein_range,protein_n_res):
    # Verify that is is possible to extract the desired residues from the protein
    if (len(protein_range) != 2) or \
       (protein_range[1] < protein_range[0]) or \
       (protein_range[0] < 1) or \
       (protein_range[1] > protein_n_res):
        raise ValueError(f"Provided protein residue range: ({protein_range[0]},{protein_range[1]}) is invalid")
        
def createIndividualMSAsFullLengthFragment(a3m_path,name,protein_range,fragment_start_range,fragment_length,protein_copies=1,fragment_single_sequence=False,fragment_shuffle_sequence=False):
    '''Loads a monomer MSA and creates new MSAs that contain 1) a large section of the monomer and 2) a short fragment of the monomer
    
    Args
    ----
    a3m_path : str
        path to an a3m formatted msa file describing the monomeric protein (obtained from colab_mmseqs)
    name : str
        name of the monomeric protein
    protein_range : list
        an inclusive range defining which positions of the MSA to take from the monomeric protein (sometimes the ends need to be chopped off)
    fragment_start_range : list
        an inclusive range defining the range of starting residues for fragments of length `fragment_length`
    fragment_length : int
        the number of residues to take when defining a fragment
    protein_copies : int
        the cardinality of the protein chain
    fragment_single_sequence : bool
        option to use single sequence for the fragment
    fragment_shuffle_sequence : bool
        option to shuffle the sequence of the fragment 
    '''
    print("Generating MSAs for a monomeric interaction...")

    # verify residue selections
    protein_n_res = readA3MProteinLength(a3m_path)
    protein_range = replaceNullProteinRange(protein_n_res,protein_range[0],protein_range[1])
    fragment_start_range = replaceNullFragmentRange(protein_n_res,fragment_length,fragment_start_range[0],fragment_start_range[1])
    verifyFragmentNterminalResRange(fragment_start_range,fragment_length,protein_n_res)
    verifyProteinRange(protein_range,protein_n_res)

    dir_name = Path(f"{name}{protein_copies}copies_tile{str(fragment_length)}aa")
    try:
        dir_name.mkdir(parents=False,exist_ok=False)
    except FileExistsError:
        print('Directory already exists, possibly overwriting existing files')

    fragment_start_iter = range(fragment_start_range[0],min(fragment_start_range[1]+1,protein_n_res-fragment_length+2))
    with Pool() as p:
        a3m_out_path_list = p.starmap(createIndividualMSAsFullLengthFragment_starmap,[(a3m_path,fragment_start,fragment_length,dir_name,name,protein_copies,protein_range,fragment_single_sequence,fragment_shuffle_sequence) for fragment_start in fragment_start_iter])

    return a3m_out_path_list

def createIndividualMSAsFullLengthFragment_starmap(a3m_path,fragment_start,fragment_length,dir_name,name,protein_copies,protein_range,fragment_single_sequence,fragment_shuffle_sequence):
    fragment_range = (fragment_start,fragment_start+fragment_length-1) # range is inclusive
    a3m_out_path = dir_name.joinpath(f"{name}{protein_copies}copies_{protein_range[0]}-{protein_range[1]}_{name}_{fragment_range[0]}-{fragment_range[1]}.a3m")
    abs_a3m_out_path = a3m_out_path.absolute()
    print(f"Creating .a3m file: {abs_a3m_out_path}")
    createMSA(a3m_path, protein_range, fragment_range, abs_a3m_out_path, -1, protein_copies, fragment_single_sequence, fragment_shuffle_sequence)
    return abs_a3m_out_path
        
def createIndividualMSAsFullLengthFragmentHeteromeric(fulllength_a3m_path,fulllength_name,fragment_a3m_path,fragment_name,protein_range,fragment_start_range,fragment_length,protein_copies=1,fragment_single_sequence=False,fragment_shuffle_sequence=False):
    '''Loads a monomer MSA and creates new MSAs that contain 1) a large section of the monomer and 2) a short fragment of the monomer
    
    Args
    ----
    fulllength_a3m_path : str
        path to an a3m formatted msa file describing the monomeric protein (obtained from colab_mmseqs)    
    fragment_a3m_path : str
        path to an a3m formatted msa file describing the protein broken into fragments (obtained from colab_mmseqs)
    name : str
        name of the monomeric protein
    protein_range : list
        an inclusive range defining which positions of the MSA to take from the monomeric protein (sometimes the ends need to be chopped off)
    fragment_start_range : list
        an inclusive range defining the range of starting residues for fragments of length `fragment_length`
    fragment_length : int
        the number of residues to take when defining a fragment
    protein_copies : int
        the cardinality of the protein chain
    fragment_single_sequence : bool
        option to use single sequence for the fragment
    fragment_shuffle_sequence : bool
        option to shuffle the sequence of the fragment 
    '''
    print("Generating MSAs for a heteromeric interaction...")

    # verify residue selections
    fulllengthprotein_n_res = readA3MProteinLength(fulllength_a3m_path)
    fragmentprotein_n_res = readA3MProteinLength(fragment_a3m_path)
    protein_range = replaceNullProteinRange(fulllengthprotein_n_res,protein_range[0],protein_range[1])
    fragment_start_range = replaceNullFragmentRange(fragmentprotein_n_res,fragment_length,fragment_start_range[0],fragment_start_range[1])
    verifyFragmentNterminalResRange(fragment_start_range,fragment_length,fragmentprotein_n_res)
    verifyProteinRange(protein_range,fulllengthprotein_n_res)

    dir_name = Path(f"{fulllength_name}{protein_copies}copies_{fragment_name}tile{str(fragment_length)}aa")
    try:
        dir_name.mkdir(parents=False,exist_ok=False)
    except FileExistsError:
        print('Directory already exists, possibly overwriting existing files')

    fragment_start_iter = range(fragment_start_range[0],min(fragment_start_range[1]+1,fragmentprotein_n_res-fragment_length+2))
    with Pool() as p:
        a3m_out_path_list = p.starmap(createIndividualMSAsFullLengthFragmentHeteromeric_starmap,[(fulllength_a3m_path,fragment_a3m_path,fragment_start,fragment_length,dir_name,fulllength_name,fragment_name,protein_copies,protein_range,fragment_single_sequence,fragment_shuffle_sequence) for fragment_start in fragment_start_iter])

    return a3m_out_path_list

def createIndividualMSAsFullLengthFragmentHeteromeric_starmap(fulllength_a3m_path,fragment_a3m_path,fragment_start,fragment_length,dir_name,fulllength_name,fragment_name,protein_copies,protein_range,fragment_single_sequence,fragment_shuffle_sequence):
    fragment_range = (fragment_start,fragment_start+fragment_length-1) # range is inclusive
    a3m_out_path = dir_name.joinpath(f"{fulllength_name}{protein_copies}copies_{protein_range[0]}-{protein_range[1]}_{fragment_name}_{fragment_range[0]}-{fragment_range[1]}.a3m")
    abs_a3m_out_path = a3m_out_path.absolute()
    print(f"Creating .a3m file: {abs_a3m_out_path}")
    createMSAHeteromicInteraction(fulllength_a3m_path, protein_range, fragment_a3m_path, fragment_range, abs_a3m_out_path, -1, protein_copies, fragment_single_sequence, fragment_shuffle_sequence)
    return abs_a3m_out_path
