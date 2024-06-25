#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh

conda create -n fragfold5 -c bioconda python==3.9 nextflow -y && \
    conda activate fragfold5 && \
    pip install . && \
    echo "Installation complete"