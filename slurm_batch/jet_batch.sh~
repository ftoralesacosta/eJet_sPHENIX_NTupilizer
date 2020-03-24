#!/bin/bash

cd /p/lustre2/ftorales/Singularity/

###___REPLACE "/p/lustre2/ftorales/" with your own directory___###
singularity exec -B cvmfs:/cvmfs -B /p/lustre2/ftorales:/hip_sphenix cvmfs/sphenix.sdcc.bnl.gov/singularity/rhic_sl7_ext.simg sh /hip_sphenix/Singularity/e_Jet_sPHENIX/llnl_batch/my_jetpyth_exe.sh $1
