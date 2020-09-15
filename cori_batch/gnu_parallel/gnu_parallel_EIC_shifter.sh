#!/bin/bash

dir="/global/project/projectdirs/alice/ftorales" #replace per user, often with ~
singularity_dir="backwards_calo_singularity/Singularity/" #parent dir in https://github.com/sPHENIX-Collaboration/Singularity
outdir="gnu_backwards_Calo_Jets"
if [[ ! -e $dir/$outdir ]]; then
    mkdir $dir/$outdir
fi

source $dir/$singularity_dir/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=$dir/$singularity_dir/install
source $dir/$singularity_dir/cvmfs/sphenix.sdcc.bnl.gov/x8664_sl7/opt/sphenix/core/bin/setup_local.sh $MYINSTALL

cd $dir/$singularity_dir/macros/macros/g4simulations/
#root -b -q "Fun4All_G4_EICDetector_NeedsSeed.C($1, $3,\"\", \"$dir/$outdir/eP_$2\")"
root -b -q "Fun4All_G4_EICDetector.C($1, $3,\"\", \"$dir/$outdir/eP_$2\")"
