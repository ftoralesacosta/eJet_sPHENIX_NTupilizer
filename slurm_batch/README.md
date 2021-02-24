Line 27 in sbatch.py should be changed to the base output sdirectory of your choice.
Line 5 in my_jetbatch_exe.sh should be changed tos INST=/FULL/PATH/Singularity/install
Line 5 in jet_batch.sh should contain the full path to the directory of the the Singularity container

sbatch.py generates .sh scripts with unique names based on the job number. These scripts then call jet_batch.sh passing an argument for the output file name and loading the container. Once inside the container, my_jetbatch_exe.sh is called that simply runs Fun4All (which has been edited to use myjetanalysis).

For 500 MB events, I requested 10 hours per job (usually submitting 50-100 jobs at a time). One could easily run less events per job for less time and simply submitt more jobs (assuming time is lineal with number of events). For higher pT-hat requirements, the number of events was around 350 per jobs for the same amount of time.
