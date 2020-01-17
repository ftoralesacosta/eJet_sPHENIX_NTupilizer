#!/bin/bash

cd /g/g12/ftorales/Singularity/

singularity exec -B cvmfs:/cvmfs -B /g/g12/ftorales:/hip_sphenix cvmfs/sphenix.sdcc.bnl.gov/singularity/rhic_sl7_ext.simg sh /hip_sphenix/batch/my_jetpyth_exe.sh $1
