#!/bin/sh
# embedded options to bsub - start with #BSUB
# -- choose queue --
#BSUB -q hpc
# -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=4GB]"
# -- Output File --
#BSUB -o Output_%J.txt
# -- Error File --
#BSUB -e Error_%J.txt
# -- estim6ted wall clock time (execution time): hh:mm -- 
#BSUB -W 18:00 
# -- Number of cores requested -- 
#BSUB -n 20
# -- Specify the distribution of the cores: on a single node --
#BSUB -R "span[hosts=1]"
# -- end of LSF options -- 

# -- commands you want to execute -- 
matlab -batch run_pcg_large_parallel > MatlabOutput_run_pcg_large_"$MATRIX_NAME"_"$NYSTROM".txt
