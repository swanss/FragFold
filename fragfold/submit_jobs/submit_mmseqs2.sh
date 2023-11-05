#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name submit_mmseqs
#SBATCH --tasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=1G
#SBATCH -o submit_mmseqs.%j.log

### -------------------------------------------------------------------- ###
# Set the following parameters

query="ftsZ.fasta"
COLABFOLD_ENV_DIR=/home/gridsan/asavinov/localcolabfold/colabfold_batch/colabfold-conda
REPO_DIR=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold

### -------------------------------------------------------------------- ###

# If values were passed as arguments, override whatever was set in this file
if [ ! -z "$1" ]; then
    query=$1
fi

if [ ! -z "$2" ]; then
    COLABFOLD_ENV_DIR=$2
fi

if [ ! -z "$3" ]; then
    REPO_DIR=$3
fi

conda run --no-capture-output -p $COLABFOLD_ENV_DIR python $REPO_DIR/mmseqs2.py --query $query