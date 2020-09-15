#!/bin/bash

let rando=$RANDOM*$2+$2
#random number seed uses built in $Random and unique gnu parallel slot number

shifter --image=docker:ddixit/fun4all:eicresearch /project/projectdirs/alice/ftorales/new_Singularity/Singularity/e_Jet_sPHENIX/cori_batch/gnu_parallel_EIC_shifter.sh $1 $2 $rando
# First argument = number of events
# Second argument = gnu parallel slot number for output file name
# Third argument = random number for MC generator SEED

#command to call this script:
#seq N_Jobs | parallel ./parallel_EIC_shifter.sh N_Events {%}

#Note that max N_Jobs should be the number of threads.
#and N_Events should be able to complete within the time you have the interactive node (max 4 hours)
