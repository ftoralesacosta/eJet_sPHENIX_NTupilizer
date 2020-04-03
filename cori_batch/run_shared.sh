#!/bin/bash
#SBATCH --image=docker:ddixit/fun4all:eicresearch
#SBATCH --qos=shared
#SBATCH --constraint=haswell
#SBATCH --time=00:10:00
#SBATCH --array=0-10
shifter ./sPHENIX_shifter.sh $SLURM_ARRAY_TASK_ID 2
