#!/bin/bash

#SBATCH --job-name=colabfold_array
#SBATCH --time=2-00:00:00
#SBATCH -e colabfold_%A_%a_err.txt
#SBATCH -o colabfold_%A_%a_out.txt
#SBATCH --qos=low
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem=32000
#SBATCH --array=0-NTOTALJOBS

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

echo "Running colabfold, job: ${SLURM_ARRAY_TASK_ID}..."

. /etc/profile.d/modules.sh
readarray -t RUNLIST < $a3m_list_file
RUNLIST_LEN=${#RUNLIST[@]}
batch_size=$(( $(($RUNLIST_LEN/$LLSUB_SIZE)) +1 ))
runfile=${repo_dir}/fragfold/submit_jobs/run_colabfold.sh # the script that will be run for each element in the file

SECONDS=0

a3m=${RUNLIST[$SLURM_ARRAY_TASK_ID]}
echo bash $runfile $a3m $colabfold_env_dir $alphafold_params_dir $cuda_version
bash $runfile $a3m $colabfold_env_dir $alphafold_params_dir $cuda_version

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED