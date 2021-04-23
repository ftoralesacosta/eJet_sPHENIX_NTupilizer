#!/bin/bash

dir="/global/project/projectdirs/alice/ftorales" #replace per user, often with ~
singularity_dir="/Singularity/" #parent dir in https://github.com/sPHENIX-Collaboration/Singularity
outdir="with_calo_gnu_allsi_reco_branches"
if [[ ! -e $dir/$outdir ]]; then
    mkdir $dir/$outdir
fi

source $dir/$singularity_dir/cvmfs/sphenix.sdcc.bnl.gov/gcc-8.3/opt/sphenix/core/bin/sphenix_setup.sh
export MYINSTALL=$dir/$singularity_dir/install
source $dir/$singularity_dir/cvmfs/sphenix.sdcc.bnl.gov/gcc-8.3/opt/sphenix/core/bin/setup_local.sh $MYINSTALL
source $dir/$singularity_dir/set_g4.sh

cd $dir/$singularity_dir/fernando_lblvtx/macros/detectors/EICDetector/
#root -b -q "Fun4All_G4_EICDetector_NeedsSeed.C($1, $3,\"\", \"$dir/$outdir/eP_$2\")"
root -b -q "Fun4All_G4_EICDetector.C($1, $2,\"\",\"$dir/$outdir/eP_$3\")"
exit
