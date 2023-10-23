#!/bin/bash

silentfile=../230803/batch_${LLSUB_RANK}.silent
args=""

rosettaDir=/data1/groups/keatinglab/rosetta/rosetta_src_2021.16.61629_bundle/main
intAnaBin=$rosettaDir/source/bin/InterfaceAnalyzer.linuxgccrelease

SECONDS=0
echo $intAnaBin -in:file:silent $silentfile -out:file:score_only -out:file:scorefile score.sc $args
$intAnaBin -in:file:silent $silentfile -out:file:score_only -out:file:scorefile score.sc $args
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED