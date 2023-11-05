import argparse
import os,sys
import pickle
import pandas as pd
from pathlib import Path

from fragfold.src.colabfold_create_msa import createIndividualMSAsFullLengthFragment,createIndividualMSAsFullLengthFragmentHeteromeric

def main(args):
    # In either mode, there will be information detailing the fragments
    fragment_a3m = args.fragment_a3m_input
    fragment_name = Path(fragment_a3m).stem
    fragment_start_range = (args.fragment_ntermres_start,args.fragment_ntermres_final)
    fragment_length = args.fragment_length
    protein_copies = args.protein_copies

    protein_range = (args.protein_ntermres,args.protein_ctermres)

    if not args.protein_a3m_input or args.protein_a3m_input == args.fragment_a3m_input:
        # Create monomeric interaction MSAs: fragments of protein A + full-length protein A
        createIndividualMSAsFullLengthFragment(fragment_a3m,
                                               fragment_name,
                                               protein_range,
                                               fragment_start_range,
                                               fragment_length,
                                               protein_copies)
    else:
        # Create heteromeric interaction MSAs: fragments of protein A + full-length protein B

        # In heteromeric interaction mode, we parse more arguments
        fulllengthprotein_a3m_path = args.protein_a3m_input
        fulllengthprotein_name = Path(fulllengthprotein_a3m_path).stem
        createIndividualMSAsFullLengthFragmentHeteromeric(fulllengthprotein_a3m_path,
                                                          fulllengthprotein_name,
                                                          fragment_a3m,
                                                          fragment_name,
                                                          protein_range,
                                                          fragment_start_range,
                                                          fragment_length,
                                                          protein_copies)

if __name__ == "__main__":
    parser = argparse()
    parser.add_argument(
        "--fragment_a3m_input",
        type=Path,
        help="",
        required=True
    )
    parser.add_argument(
        "--fragment_ntermres_start",
        type=int,
        help="",
    )
    parser.add_argument(
        "--fragment_ntermres_final",
        type=int,
        help="",
    )
    parser.add_argument(
        "--fragment_length",
        type=int,
        help="",
        default=30
    )
    parser.add_argument(
        "--protein_a3m_input",
        type=Path,
        help="",
        default=None
    )
    parser.add_argument(
        "--protein_ntermres",
        type=int,
        help="",
        default=None
    )
    parser.add_argument(
        "--protein_ctermres",
        type=int,
        help="",
        default=None
    )
    parser.add_argument(
        "--protein_copies",
        type=int,
        help="",
        default=1
    )
    args = parser.parse_args()
    main(args)