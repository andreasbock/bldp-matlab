#!/bin/sh

## declare an array variable
#declare -a names=("494_bus" "1138_bus") # testing
declare -a names=("Pres_Poisson" "bmwcra_1" "cfd1" "smt" "apache1" "thermal1" "crankseg_1" "crankseg_2" "bcsstk17" "bcsstk18" "consph")

## now loop through the above array
for i in "${names[@]}"
do
   echo "$i"
   bsub -env "MATRIX_NAME=$i" < scripts/submit_matrix.sh;
done
