#!/bin/bash
#SBATCH -J interfaceAnalyzer
#SBATCH -n 1 #Request 1 task (core)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 0-6:00
#SBATCH --mem-per-cpu=4G #Request 4G of memory per CPU
#SBATCH -o interfaceAnalyzer.%J.out #redirect output to output_JOBID.txt
#SBATCH -e interfaceAnalyzer.%J.err #redirect errors to error_JOBID.txt

structureList=structure_list.txt
chain="0"
packstat_oversample=1 #increase to 100 for higher accuracy, but considerably longer runtime

rosettaDir=/data1/groups/keatinglab/rosetta/rosetta_src_2021.16.61629_bundle/main
intAnaBin=$rosettaDir/source/bin/InterfaceAnalyzer.linuxgccrelease

SECONDS=0

echo $intAnaBin -in:file:l $structureList -add_regular_scores_to_scorefile -compute_packstat -packstat::oversample $packstat_oversample -tracer_data_print
$intAnaBin -in:file:l $structureList -add_regular_scores_to_scorefile -compute_packstat -packstat::oversample $packstat_oversample -tracer_data_print

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
