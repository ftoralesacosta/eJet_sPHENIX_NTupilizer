#!/bin/sh
ls -l
date

INST=/hip_sphenix/Singularity/install
source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh $INST

cd /hip_sphenix/Singularity/macros/macros/g4simulations
pwd

echo "root -b -q Fun4All_G4_sPHENIX.C(500, \"__\", ${1})"
root -b -q Fun4All_G4_sPHENIX.C\(500,\"__\",\"${1}\"\) #Fun4All second argument is input HEPMC file. Current focus is on generating

date
ls -l

