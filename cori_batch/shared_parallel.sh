#!/bin/bash
#SBATCH --image=docker:ddixit/fun4all:eicresearch
#SBATCH --qos=regular
#SBATCH --constraint=haswell
#SBATKCH --Nodes=4
module load parallel
srun parallel --jobs 64 
