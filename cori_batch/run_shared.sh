#!/bin/bash
#SBATCH --image=docker:ddixit/fun4all:eicresearch
#SBATCH --qos=shared
#SBATCH --constraint=haswell
#SBATCH --time=5:00:00
#SBATCH --array=0-40
shifter ./sPHENIX_shifter.sh $SLURM_ARRAY_TASK_ID 250
