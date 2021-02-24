These scripts can be used to run Fun4All on all threads in a single node. Any other script can be called from "parallel_EIC_shifter.sh" very easily, with the caveat that you are limited to a single node where the memory is shared.

For multi-node use, please see the mpi_parallel directory.

To run on CORI:

`salloc -N 1 -C haswell -q interactive -t 04:00:00`
`module load parallel`
`seq N_Jobs | parallel ./parallel_EIC_shifter.sh N_Events {%}`

For cori haswell, the number of jobs should be 64, and the number of events should finish within the allocated time (max 4 hours). The EIC detector similation is often faster than sPHENIX. For EIC, 64 jobs at about 1000 events maximum should be ok. sPHENIX should be no more than approximatley 300 events per job to be safe.

P.S. for short simulations, one can use the debug qos instead of interactive. You have less time, but it does not count towards your quota (and therefore won't impact future queue times).