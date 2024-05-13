import glob, json

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio.PDB import *

# The following functions are specific to processing colabfold output
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

