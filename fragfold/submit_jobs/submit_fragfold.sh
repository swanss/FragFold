#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name fragfold
#SBATCH --tasks-per-node=1
#SBATCH --time=3-00:00:00
#SBATCH --mem=8G
#SBATCH -o fragfold.%j.log

### -------------------------------------------------------------------- ###
# Set the following parameters

REPO_DIR=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold
COLABFOLD_ENV_DIR=~/localcolabfold/colabfold-conda
CONDA_ENV_NAME=fragfold

# MSA generation
query_seq=ftsZ.fasta # can be a single .fasta, or a directory containing many .fasta

# MSA processing
fragment_a3m_name="ftsZ" # should match the name of the fasta file from which fragments will be derived
fullprotein_a3m_name="ftsZ" # should match the name of fasta file from the full-length protein will be derived
fragment_ntermres_start=1
fragment_ntermres_final=383
fragment_length=30
protein_ntermres=10
protein_ctermres=316
protein_copies=1

### -------------------------------------------------------------------- ###

# The following script can be used to run FragFold
# 1. Generate MSA(s) with MMseqs2
# 2. Process MSAs to create fragment + full-length protein MSAs
# 3. Submit ColabFold jobs
# After all jobs are completed, the results can be analyzed

### Generate MSAs
# NOTE: this needs to be run in an environment with internet connection
if bash ${REPO_DIR}/fragfold/submit_jobs/submit_mmseqs2.sh $query_seq $COLABFOLD_ENV_DIR $REPO_DIR; then
    echo "Generated MSA(s) with MMseqs2"
else
    echo "Unable to generate MSA(s), terminating."
    exit 1
fi

### Create fragment + full-length MSAs
if conda run --no-capture-output -n $CONDA_ENV_NAME python ${REPO_DIR}/fragfold/create_fragment_msa.py \
    --fragment_a3m_input ${PWD}/mmseqs2_a3m/{$fragment_a3m_name}.a3m \
    --fragment_ntermres_start $fragment_ntermres_start \
    --fragment_ntermres_final $fragment_ntermres_final \
    --fragment_length $fragment_length \
    --protein_a3m_input ${PWD}/mmseqs2_a3m/{$fullprotein_a3m_name}.a3m \
    --protein_ntermres $protein_ntermres \
    --protein_ctermres $protein_ctermres \
    --protein_copies $protein_copies \
    ; then 
    echo "Created fragment + full-length processed MSAs"
else
    echo "Unable to process MSA(s), terminating."
    exit 1
fi

echo "Success: submitting ColabFold jobs and terminating..."