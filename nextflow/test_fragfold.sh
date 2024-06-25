#!/bin/bash

NEXTFLOWDIR=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold/nextflow #directory containing nextflow scripts
WORKDIR=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold/example #directory where results will be stored
NF_CFG=${NEXTFLOWDIR}/nextflow.config
PARAMS=${NEXTFLOWDIR}/params/ftsZ_monomeric_example.yaml
mkdir -p $WORKDIR
nextflow run ${NEXTFLOWDIR}/main.nf -w $WORKDIR -c $NF_CFG -params-file $PARAMS -resume
