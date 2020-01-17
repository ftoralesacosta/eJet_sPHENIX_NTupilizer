from __future__ import division
from __future__ import print_function
import os
import shutil
import argparse
import subprocess

# -------------------------- ARGUMENT PARSING ------------------------ #

prog_description = 'Creates and submits job files to run jetquench chain.'
prog_epilogue    = 'And that\'s all, folks!'

# Initialize argument parser
parser = argparse.ArgumentParser(
    description=prog_description,
    usage='%(prog)s [options]',
    epilog=prog_epilogue,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Add arguments
#
# Arguments 
#
parser.add_argument('-c', '--config',type=str,default='pp',
                    help='base for config files')
parser.add_argument('-d', '--dir',type=str,
                    default='/p/lustre2/ftorales/Singularity/',
                    help='top level scratch directory')
parser.add_argument('-j', '--jobset',type=str,default='MinBias_MC',
                    help='jobset name for this jobset [MinBias_MC]')
parser.add_argument('-n', '--njob', type=int,default='10',
                    help='number of jobs to submit')
parser.add_argument('-s', '--submit',action='store_true',
                    help='submit jobs with sbatch after creating')
parser.add_argument('-v', '--verbose',action='store_true',
                    help='turn on diagnostic print statements')
parser.add_argument('-t', '--time',type=str,default = '30',
                    help='minutes for job submission')

# Parse arguments
args = parser.parse_args()

# Check validity of arguments and environment
if not os.path.isdir(args.dir):
    raise ValueError('directory {} not found.'.format(args.dir))

configfile = args.config + '.config'

# ------------- CREATE and FILL DIRECTORY TREE ---------------------- #
# Note we follow mkdir with chmod to circumvent umask restrictions

topdir = os.path.join(args.dir,args.jobset)
if (args.verbose):
    print('topdir = ', topdir)
if not os.path.isdir(topdir):
    os.mkdir(topdir) 
    os.chmod(topdir,0o2770) 

#

nvals = range(args.njob)
for n in nvals:
    nstr = '{:03}'.format(n) # pad zeros up to 1000 jobs
    jobname = args.config + '_' + nstr
    jobdir = os.path.join(topdir,jobname)
    if (args.verbose):
        print ('jobdir = ', jobdir)
    if not os.path.isdir(jobdir):
        os.mkdir(jobdir)
        os.chmod(jobdir,0o2770)

# write job script
    basefile = os.path.join(jobdir,jobname)
    batchfile = basefile + '.sh'
    outfile = basefile + '.out'
    errfile = basefile + '.err'
    f = open(batchfile,'w')
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -N 1\n')
    f.write('#SBATCH -t '+args.time+'\n')
    #f.write('#SBATCH -t 00:'+args.time+':00:00\n')
    f.write('#SBATCH -p pbatch\n')
    f.write('#SBATCH -A hizphys\n')
    f.write('#SBATCH -o '+outfile+'\n')
    f.write('#SBATCH -e '+errfile+'\n')
    f.write('#\n# Change directory and setup runtime environment\n#\n')
    f.write('cd '+jobdir+'\n')
    f.write('#\n# Run the MB MC sPHENIX chain\n#\n')
    f.write('srun -n 1 /g/g12/ftorales/batch/jet_batch.sh '+jobdir+'.root\n')
    f.write('date\n')
    f.close()
    
    # submit job
    if (args.submit):
        submit_command = 'sbatch ' + batchfile
        #submit_command = 'ls ' + batchfile
        print (submit_command)
        os.system(submit_command)
        
#        print (submit_command)
#        os.system(submit_command)
#subprocess.call(['sbatch',batchfile])
