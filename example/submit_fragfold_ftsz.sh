#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name fragfold
#SBATCH --tasks-per-node=1
#SBATCH --time=3-00:00:00
#SBATCH --mem=8G
#SBATCH -o fragfold.%j.log

### -------------------------------------------------------------------- ###
# Set the following parameters

# Environments and installation paths
repo_dir=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold
colabfold_dir=~/localcolabfold
conda_env_name=fragfold
cuda_version=11.1

# FASTA files that will be used to generate MSA
query_seq=ftsZ.fasta 

# MSA processing parameters
fragment_a3m_name=ftsZ 
fullprotein_a3m_name=ftsZ 
fragment_ntermres_start=1
# fragment_ntermres_final=353
fragment_ntermres_final=10
fragment_length=30
protein_ntermres=10
protein_ctermres=316
protein_copies=1

# Job array mode:
array_mode=slurm_array #alternative: "llsub"

# LLsub job submission (these do not matter if running in slurm_array mode)
n_nodes=4
n_gpu=2
n_threads=1

### -------------------------------------------------------------------- ###

### parameters that are defined automatically from user-provided parameters###
colabfold_env_dir=$colabfold_dir"/colabfold-conda"
alphafold_params_dir=$colabfold_dir"/colabfold"
### ###

# The following script can be used to run FragFold
# 1. Generate MSA(s) with MMseqs2
# 2. Process MSAs to create fragment + full-length protein MSAs
# 3. Submit ColabFold jobs
# After all jobs are completed, the results can be analyzed

### Generate MSAs
# NOTE: this needs to be run in an environment with internet connection
if bash ${repo_dir}/fragfold/submit_jobs/submit_mmseqs2.sh $query_seq $colabfold_env_dir $repo_dir; then
    echo "Generated MSA(s) with MMseqs2"
else
    echo "Error: unable to generate MSA(s), terminating..."
    exit 1
fi

### Create fragment + full-length MSAs
if conda run --no-capture-output -n $conda_env_name python ${repo_dir}/fragfold/create_fragment_msa.py \
    --fragment_a3m_input ${PWD}/mmseqs2_a3m/${fragment_a3m_name}.a3m \
    --fragment_ntermres_start $fragment_ntermres_start \
    --fragment_ntermres_final $fragment_ntermres_final \
    --fragment_length $fragment_length \
    --protein_a3m_input ${PWD}/mmseqs2_a3m/${fullprotein_a3m_name}.a3m \
    --protein_ntermres $protein_ntermres \
    --protein_ctermres $protein_ctermres \
    --protein_copies $protein_copies \
    ; then 
    echo "Created fragment + full-length processed MSAs"
else
    echo "Error: unable to process MSA(s), terminating..."
    exit 1
fi

### Submit ColabFold jobss
if [[ $array_mode == "slurm_array" ]]; then
    # Create slurm array job submission script
    cp ${repo_dir}/fragfold/submit_jobs/submit_colabfold_slurmarray.sh .

    # replace variables
    a3m_list_file=$(ls $(pwd -P)/a3m_list.txt)
    
    readarray -t RUNLIST < $a3m_list_file
    ntotaljobs=${#RUNLIST[@]}
    ntotaljobs=$(($ntotaljobs-1))

    sed -i -e "s/NTOTALJOBS/$ntotaljobs/g" submit_colabfold_slurmarray.sh
    sed -i -e "s#A3M_LIST_FILE#$a3m_list_file#g" submit_colabfold_slurmarray.sh
    sed -i -e "s#COLABFOLD_ENV_DIR#$colabfold_env_dir#g" submit_colabfold_slurmarray.sh
    sed -i -e "s#ALPHAFOLD_PARAMS_DIR#$alphafold_params_dir#g" submit_colabfold_slurmarray.sh
    sed -i -e "s#REPO_DIR#$repo_dir#g" submit_colabfold_slurmarray.sh
    sed -i -e "s/CUDA_VERSION/$cuda_version/g" submit_colabfold_slurmarray.sh

    echo "Submitting ColabFold jobs in SLURM array mode..."
    if sbatch submit_colabfold_slurmarray.sh; then
        echo "Successfully submitted jobs"
    else
        echo "Error: unable to submit jobs, terminating..."
    fi
elif [[ $array_mode == "llsub" ]]; then
    # Create llsub job submission script
    cp ${repo_dir}/fragfold/submit_jobs/submit_colabfold_triples.sh .

    # replace variables
    a3m_list_file=$(ls $(pwd -P)/a3m_list.txt)
    sed -i -e "s#A3M_LIST_FILE#$a3m_list_file#g" submit_colabfold_triples.sh
    sed -i -e "s#COLABFOLD_ENV_DIR#$colabfold_env_dir#g" submit_colabfold_triples.sh
    sed -i -e "s#ALPHAFOLD_PARAMS_DIR#$alphafold_params_dir#g" submit_colabfold_triples.sh
    sed -i -e "s#REPO_DIR#$repo_dir#g" submit_colabfold_triples.sh
    sed -i -e "s/CUDA_VERSION/$cuda_version/g" submit_colabfold_triples.sh

    echo "Submitting ColabFold jobs in triples mode..."
    chmod u+x submit_colabfold_triples.sh
    echo $PWD
    echo $CWD
    echo "LLsub ./submit_colabfold_triples.sh "["$n_nodes","$n_gpu","$n_threads"]" -g volta:1"
    if LLsub ./submit_colabfold_triples.sh "["$n_nodes","$n_gpu","$n_threads"]" -g volta:1; then
        echo "Successfully submitted jobs"
    else
        echo "Error: unable to submit jobs, terminating..."
    fi
else
    echo "array mode = "$array_mode" is not recognized, terminating..."
fi
