On supercloud, request an interactive node

```
LLsub -i --time=1-00:00 --partition=xeon-p8
```

Go to the local directory

```
USER=$(whoami) && mkdir -p /state/partition1/user/$USER && cd /state/partition1/user/$USER
```

Run the command

```
NEXTFLOWDIR=/home/gridsan/sswanson/keatinglab_shared/swans/savinovCollaboration/FragFold/nextflow
WORKDIR=/home/gridsan/sswanson/keatinglab_shared/swans/savinovCollaboration/FragFold/nextflow/practice
nextflow run ${NEXTFLOWDIR}/fragfold_dummy.nf -w $WORKDIR
```

```
FASTA=/home/gridsan/sswanson/keatinglab_shared/swans/savinovCollaboration/FragFold/example/ftsZ.fasta
NEXTFLOWDIR=/home/gridsan/sswanson/keatinglab_shared/swans/savinovCollaboration/FragFold/nextflow
WORKDIR=/home/gridsan/sswanson/keatinglab_shared/swans/savinovCollaboration/FragFold/nextflow/practice
cp $FASTA .
nextflow run ${NEXTFLOWDIR}/fragfold_dummy_2.nf -w $WORKDIR
```

Get information about the pipeline 

```
nextflow log $RUN_NAME -f name,status,workdir | tee nextflow_report.txt
```

final command 

```
nextflow run ${NEXTFLOWDIR}/fragfold_example.nf -w $WORKDIR -resume
```