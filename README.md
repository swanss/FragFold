# Fragment binding prediction with ColabFold

Scripts for predicting how short fragments of natural proteins bind to full-length proteins, as described in the [manuscript](link/to/paper). This program is built on top of MMseqs2 and ColabFold, extending them to efficiently predict interactions between a full-length protein and fragments derived from a protein.

# Installation

## Overview of the steps

### 1. Install local colabfold

Follow the instructions at [https://github.com/YoshitakaMo/localcolabfold]() to install ColabFold on your system. We developed and and tested FragFold on this specific [commit](https://github.com/YoshitakaMo/localcolabfold/tree/88d174ffa7a7bc76a644db14ba0099ceb0606aed).

### 2. Clone fragfold repo

Clone the repo to your computer via ssh or https:

```
git clone INSERTPATHHERE
cd fragfold
```

### 3. Set up environment

Using conda or the faster alternative [mamba](https://github.com/mamba-org/mamba), install the FragFold  environment.

`conda env create -f fragfold.yml`

### NOTE

As an alternative to the above steps, you can also use the provided docker container: `fragfold_img`

# Examples

The following are examples derived from the manuscript.

## Predict homomeric interactions

In this example we will predict interactions between a full-length protein and fragments derived from it.

### Generate MSA describing the full-length protein using MMseqs2

#### Option 1: Use a Google ColabFold Notebook

1. Select your protein sequence
2. Use google colabfold or a local install of mmseqs to generate a msa
3. copy the MSA to dir

#### Option 2: Use a local install of MMseqs2


### Process the MSA to create the fragment+protein MSAs

1. Process the MSA using `script.py`

### Run colabfold job

2. Run the jobs in parallel using 'script.py'

## Predict heterometic interactions

In this example we will predict interactions between a full-length protein and fragments derived from another protein.

### Generate MSA describing the full-length protein using MMseqs2

### Process the MSA to create the fragment+protein MSAs

### Run colabfold job

