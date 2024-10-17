#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh

conda create -n fragfold -c bioconda python==3.9 nextflow -y && \
    conda activate fragfold && \
    pip install . && \
    echo "Installation complete"
