# Fragment binding prediction with ColabFold

Scripts for predicting how short fragments of natural proteins bind to full-length proteins, as described in the [manuscript](link/to/paper). This program is built on top of MMseqs2 and ColabFold, extending them to efficiently predict interactions between a full-length protein and fragments derived from a protein.

# Installation

## Overview of the steps

### 1. Install local colabfold

Follow the instructions at the local colabfold [repo](https://github.com/YoshitakaMo/localcolabfold) to install ColabFold on your system. We developed and tested FragFold on this specific [commit](https://github.com/YoshitakaMo/localcolabfold/tree/88d174ffa7a7bc76a644db14ba0099ceb0606aed).

### 2. Install fragfold

Use the following commands to clone the repo and install FragFold and its dependencies on your computer.

```
git clone https://github.com/swanss/FragFold.git
cd FragFold
conda create -n fragfold python=3.7
conda activate fragfold
pip install .
```

### NOTE

As an alternative to the above steps, you can also use the provided docker container: `fragfold_img`

# Examples

The following are examples derived from the manuscript.

## Predict homomeric interactions between full-length FtsZ and fragments

If you're working on a cluster where the compute nodes have internet access

```
cd /FragFold/example
cp ../fragfold/submit_jobs/submit_fragfold.sh .
sbatch submit_fragfold.sh 
```

Otherwise, you will need to run the MSA generation stage separately in an environment with internet access (e.g. on the login node)

```
cd /FragFold/example
cp ../fragfold/submit_jobs/submit_fragfold.sh .
cp ../fragfold/submit_jobs/submit_mmseqs2.sh .
bash submit_mmseqs2.sh
sbatch submit_fragfold.sh 
```

