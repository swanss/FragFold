#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name submit_colabfold
#SBATCH --tasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=1G
#SBATCH -o submit_colabfold.%j.log

# Job submission script with the following features:
# - Waits between submissions to keep the number of pending jobs from going too high
# - Checks to make sure that AssocMaxSubmitJobLimit is not exceeded before submitting more jobs
# - Checks to make sure a job is completed, and if not, cleans the directory and resubmits

max_jobs=150
msa_list=/home/gridsan/sswanson/local_code_mirror/ColabFoldCustomMSA/lptG_tile14aa_list.txt
#sbatch=$PWD/dummy_test.sh
sbatch=$PWD/run_colabfold.sh

result_dir=$PWD/data
if [ ! -d $result_dir ]
then
	echo mkdir $result_dir
        mkdir $result_dir
fi
cd $result_dir

SECONDS=0
COUNTER=1
while IFS= read -r path; do
	echo $COUNTER
	echo $path
	## Decide whether it is time to submit the job
	val=$(squeue | wc -l)
	#this is total jobs: could be in PD/R state, or could be other types of jobs
	n_jobs=$(($val - 1))
	echo $n_jobs
	echo $max_jobs
	while [[ $n_jobs -ge $max_jobs ]]; do
		echo "sleep 5"
		sleep 5m
		val=$(squeue | wc -l)
		n_jobs=$(($val - 1))
	done

	## Check if the job has completed
	file=${path##*/}
	name=${file%.a3m}
	echo $name
	
	# Check if the directory exists
	if [ -d $name ]; then
		echo "directory exists"
		if test -f "$name/output/$name.done.txt"; then
			echo "job is completed, do not resubmit"
			continue
		else
			echo "job did not complete, delete directory and resubmit"
			rm -r $name
		fi
	fi
	
	# (re)submit job
	mkdir $name
	cd $name
	echo sbatch $sbatch $path
	sbatch $sbatch $path
	cd $result_dir #return to main using absolute path to avoid propagating errors
	COUNTER=$[COUNTER +1]
done < $msa_list

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED
