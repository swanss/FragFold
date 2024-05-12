#!/bin/bash

conda create -n fragfold4 -c bioconda python==3.9 nextflow -y 
    && conda activate fragfold4
    && pip install . 
    && echo "Installation complete"