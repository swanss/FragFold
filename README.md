# Fragment binding prediction with ColabFold

Scripts for predicting how short fragments of natural proteins bind to full-length proteins, as described in the [manuscript](https://www.biorxiv.org/content/10.1101/2023.12.19.572389v1). This program is built on top of MMseqs2 and ColabFold, extending them to efficiently predict interactions between a full-length protein and fragments derived from a protein.

This code is associated with the following article:
A. Savinov, S. Swanson, A. E. Keating, G.-W. Li. High-throughput computational discovery of inhibitory protein fragments with AlphaFold. bioRxiv (2023). doi: 10.1101/2023.12.19.572389. https://www.biorxiv.org/content/10.1101/2023.12.19.572389v1

Please cite this article if you make use of FragFold.

Associated Source Data can also be found here: 
https://figshare.com/articles/dataset/Source_Data_for_Savinov_and_Swanson_et_al_2023/24841269  
doi: 10.6084/m9.figshare.24841269

# Installing FragFold

## Overview of the steps

### 1. Install local colabfold

Follow the instructions at the local colabfold [repo](https://github.com/YoshitakaMo/localcolabfold) to install ColabFold on your system. We developed and tested FragFold on this specific [commit](https://github.com/YoshitakaMo/localcolabfold/tree/88d174ffa7a7bc76a644db14ba0099ceb0606aed).

### 2. Install FragFold

Use the following commands to clone the repo and install FragFold and its dependencies on your computer.

```bash
git clone https://github.com/swanss/FragFold.git
cd FragFold
conda create -n fragfold python=3.7
conda activate fragfold
pip install .
```

# Running FragFold

To run the FragFold pipeline, you can submit a single bash script (`fragfold/submit_jobs/submit_fragfold.sh`) to run all the steps together. 

## FtsZ example

### Running submit_fragfold.sh

Navigate to the example directory (`example`) and set the variables in `submit_fragfold_ftsz.sh`

```bash
# Environments and installation paths
repo_dir=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold # path to the cloned repo
colabfold_env_dir=~/localcolabfold/ # path to where localcolabfold was installed
conda_env_name=fragfold # name of the conda environment
cuda_version=11.1 

# FASTA files that will be used to generate MSA
query_seq=ftsZ.fasta # can be a single .fasta, or a directory containing many .fasta

# MSA processing parameters
fragment_a3m_name=ftsZ # should match the name of the fasta file from which fragments will be derived
fullprotein_a3m_name=ftsZ # should match the name of fasta file from the full-length protein will be derived
fragment_length=30 # the number of residues in each fragment
fragment_ntermres_start=1 # the first residue of the first fragment
fragment_ntermres_final=354 # the first residue of the last fragment (set this to ~10 for a quick test)
protein_ntermres=10 # the first residue of the protein (generally 1)
protein_ctermres=316 # the last residue of the protein 
protein_copies=1

# Job array mode:
array_mode=slurm_array #alternative for supercloud: "llsub"

# LLsub job submission (these do not matter if running in slurm_array mode)
n_nodes=4
n_gpu=2
n_threads=1
```

If you're working on a cluster where the compute nodes have internet access, run `submit_fragfold.sh` after copying/editing the file with your arguments.

```bash
cd /FragFold/example
sbatch submit_fragfold_ftsz.sh 
```

Otherwise, you will need to run the MSA generation stage separately in an environment with internet access (e.g. on the login node). After the MSAs are ready, the mmseqs step will be bypassed when running the bash script. Note that if you change the fasta files you will need to rerun the mmseqs script.

```bash
cd /FragFold/example
cp ../fragfold/submit_jobs/submit_mmseqs2.sh .
bash submit_mmseqs2.sh
sbatch submit_fragfold_ftsz.sh 
```

### Process output

First, create a `.json` file with paths to the output of the FragFold jobs.

```json
{
    "ftsZ-coding-EcoliBL21DE3": # the full-length protein name
        {
        "30aa_monomer_ftsZ": # the name of the fragment + full-length protein screen
            {
            "colabfold":"REPODIR/example" #replace with absolute path to directory where colabfold jobs were submitted
            }
        }
}
```

Copy the script to the directory, set the input variables, and then run `colabfold_process_output.py`

```bash
cd example/processoutput
sbatch run_colabfold_process_output.sh
```

### Downstream analysis

Given the colabfold predictions, you can predict peaks by running `predict_alphafold_peaks.py`.

```bash
cd example/ftsZ_predictpeaks
sbatch run_predict_alphafold_peaks.sh
```