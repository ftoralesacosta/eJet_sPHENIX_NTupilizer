#!/bin/bash
#SBATCH --image=docker:ddixit/fun4all:eicresearch
#SBATCH --qos=shared
#SBATCH --constraint=haswell
#SBATCH --time=1:00:00
#SBATCH --array=0-2999
shifter ./slurm_EIC_shifter.sh 500 $SLURM_ARRAY_TASK_ID*$RANDOM+$RANDOM $SLURM_ARRAY_TASK_ID 
