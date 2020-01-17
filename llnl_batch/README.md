Line 27 in sbatch.py should be changed to the directory of your choice.

sbatch.py generates .sh scripts with unique names based on the job number. These scripts then call other two scripts in this directory to load the container and run Fun4All MC production.
