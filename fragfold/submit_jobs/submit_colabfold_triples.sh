#!/bin/bash

. /etc/profile.d/modules.sh
list_file=/home/gridsan/sswanson/local_code_mirror/inhibitory_fragments_structure_prediction/groEL2copies_tile30aa_list.txt
readarray -t RUNLIST < $list_file
RUNLIST_LEN=${#RUNLIST[@]}
batch_size=$(( $(($RUNLIST_LEN/$LLSUB_SIZE)) +1 ))
runfile=$PWD/run_colabfold_triples.sh # the script that will be run for each element in the file

# get the batch boundaries
let start=$batch_size*$LLSUB_RANK
let next=$LLSUB_RANK+1
let next=$batch_size*$next
if [[ $next -gt $RUNLIST_LEN ]]
then
        let end=$RUNLIST_LEN
else
        let end=$next
fi

SECONDS=0
# run the batch
i=$start
while [[ $i -lt $end ]]
do
        element=${RUNLIST[$i]}
        echo bash $runfile $element
        bash $runfile $element
        i=$(($i + 1))
done

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED