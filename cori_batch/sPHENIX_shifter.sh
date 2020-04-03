#!/bin/bash

dir="/global/project/projectdirs/alice/ftorales" #replace per user, often with ~
#dir="~"
if [[ ! -e $dir/eP_MC ]]; then
    mkdir $dir/eP_MC
fi

source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=$dir/Singularity/install
source $dir/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_local.sh $MYINSTALL

cd $dir/Singularity/macros/macros/g4simulations/
root -b -q "Fun4All_G4_sPHENIX.C($2, \"\", \"$dir/eP_MC/eP_$1\")"
