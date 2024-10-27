#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name nf-manager
#SBATCH --tasks-per-node=1
#SBATCH --time=4-04:00:00
#SBATCH --mem=1G
#SBATCH --partition=xeon-p8
#SBATCH --output=nf-fragfold.%j.log
#SBATCH --signal=B:USR1@600

ENV=fragfold
WORKFLOW=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold/nextflow/main.nf
NF_CFG=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold/nextflow/nextflow.config
PARAMS=/data1/groups/keatinglab/swans/savinovCollaboration/FragFold/nextflow/params/ftsZ_monomeric_example.yaml
WORK_DIR=$(pwd -P)/nextflow_files
LOGS_DIR=${WORK_DIR}/nextflow_logs

if [ ! -d "$WORK_DIR" ]; then
    mkdir -p $WORK_DIR
fi

# Go to a tmp dir
USER=$(whoami) && cd $TMPDIR
echo "tmpdir: "$TMPDIR

# If this has already been run, copy the logs back to the tmp dir so that we can resume
if [ -d "$LOGS_DIR" ]; then
    cp -r $LOGS_DIR/.nextflow .
    cp $LOGS_DIR/.nextflow.log* .
else
    mkdir -p ${LOGS_DIR}
fi

# Terminate early if job is about to run out of time to ensure that metadata is copied so that the job can be resumed
trap 'echo Signal USR1 received! terminating early and storing progress before walltime ends; kill ${PID}; wait ${PID}' USR1
conda run -n $ENV --no-capture-output nextflow run $WORKFLOW -w $WORK_DIR -c $NF_CFG -params-file $PARAMS -resume & # launch my_script as a background job
PID=$!; echo "waiting for PID: "$PID; wait ${PID}

echo 'Finished job'

cp -r .nextflow* ${LOGS_DIR} && \
    echo "Copied logs from $TMPDIR"

cp *.csv $WORK_DIR/.. && \
    echo "Copied csv from $TMPDIR"

cp -r --preserve=links output_colabfold $WORK_DIR/.. && \
    echo "Copied colabfold files from $TMPDIR"

cp -r --preserve=links output_peakprediction $WORK_DIR/.. && \
    echo "Copied peakprediction files from $TMPDIR"
