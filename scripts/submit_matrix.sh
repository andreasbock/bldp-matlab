#!/bin/sh
# embedded options to bsub - start with #BSUB
# -- choose queue --
#BSUB -q hpc
# -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=16GB]"
# -- Output File --
#BSUB -o Output_%J.txt
# -- Error File --
#BSUB -e Error_%J.txt
# -- estimated wall clock time (execution time): hh:mm -- 
#BSUB -W 08:00 
# -- Number of cores requested -- 
#BSUB -n 1
# -- Specify the distribution of the cores: on a single node --
#BSUB -R "span[hosts=1]"
# -- end of LSF options -- 

# -- commands you want to execute -- 
matlab -batch run_pcg_matrix > MatlabOutput_run_pcg_matrix_"$MATRIX_NAME"_"$NYSTROM".txt
