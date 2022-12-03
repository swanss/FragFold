#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name dummy_test
#SBATCH --tasks-per-node=1
#SBATCH --time=0:10:00
#SBATCH -o test.%J.log

echo "test"
sleep 5m
