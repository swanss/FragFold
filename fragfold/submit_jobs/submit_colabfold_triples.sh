#!/bin/bash

a3m_list_file=A3M_LIST_FILE
colabfold_env_dir=COLABFOLD_ENV_DIR
alphafold_params_dir=ALPHAFOLD_PARAMS_DIR
repo_dir=REPO_DIR
cuda_version=CUDA_VERSION

echo "a3m_list_file: "$a3m_list_file
echo "colabfold_env_dir: "$colabfold_env_dir
echo "alphafold_params_dir: "$alphafold_params_dir
echo "repo_dir: "$repo_dir
echo "cuda_version: "$cuda_version

echo "Submitting colabfold jobs..."

. /etc/profile.d/modules.sh
readarray -t RUNLIST < $a3m_list_file
RUNLIST_LEN=${#RUNLIST[@]}
batch_size=$(( $(($RUNLIST_LEN/$LLSUB_SIZE)) +1 ))
runfile=${repo_dir}/fragfold/submit_jobs/run_colabfold.sh # the script that will be run for each element in the file

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

echo "total number of jobs: "$RUNLIST_LEN
echo "batch size: "$batch_size
echo "start: "$start
echo "end: "$end

SECONDS=0
# run the batch
i=$start
while [[ $i -lt $end ]]
do
        a3m=${RUNLIST[$i]}
        echo bash $runfile $a3m $colabfold_env_dir $alphafold_params_dir $cuda_version
        bash $runfile $a3m $colabfold_env_dir $alphafold_params_dir $cuda_version
        i=$(($i + 1))
done

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED