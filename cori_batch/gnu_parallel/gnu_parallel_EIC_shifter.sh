#!/bin/bash

dir="/global/project/projectdirs/alice/ftorales" #replace per user, often with ~
#dir="~"
outdir="gnu_forward_Calo_Jets"
if [[ ! -e $dir/$outdir ]]; then
    mkdir $dir/$outdir
fi

source $dir/new_Singularity/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=$dir/new_Singularity/Singularity/install
source $dir/new_Singularity/Singularity/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_local.sh $MYINSTALL

cd $dir/new_Singularity/Singularity/macros/macros/g4simulations/
root -b -q "Fun4All_G4_EICDetector_NeedsSeed.C($1, $3,\"\", \"$dir/$outdir/eP_$2\")"
