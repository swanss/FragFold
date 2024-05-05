#!/bin/bash

conda create -n fragfold -c bioconda python=3.7 nextflow -y 
    && conda activate fragfold 
    && pip install . 
    && echo "Installation complete"