#!/bin/bash

listfile=/data1/groups/keatinglab/swans/savinovCollaboration/analysis/rosetta_relax/relax_lists/batch_${LLSUB_RANK}_list.txt
script=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold/fragfold/rosetta/relaxscript_modified_InterfaceRelax2019.txt
silentfile_name="structures_${LLSUB_RANK}.bin"
args="--relax:script $script -constrain_relax_to_start_coords"

rosettaDir=/data1/groups/keatinglab/rosetta/rosetta_src_2021.16.61629_bundle/main
relaxBin=$rosettaDir/source/bin/relax.linuxgccrelease

SECONDS=0
echo $relaxBin -in:file:l $listfile -out:file:silent $silentfile_name -out:file:scorefile score.sc $args
$relaxBin -in:file:l $listfile -out:file:silent $silentfile_name -out:file:scorefile score.sc $args
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED