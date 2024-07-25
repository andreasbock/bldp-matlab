#!/bin/sh

## declare an array variable
declare -a nystrom=("0" "1")
#declare -a names=("Pres_Poisson" "bmwcra_1" "cfd1" "smt" "apache1" "thermal1" "crankseg_1" "crankseg_2" "bcsstk17" "bcsstk18" "consph")
declare -a names=("494_bus", "1138_bus")

## now loop through the above array
for j in "${nystrom[@]}"
do
    echo "$j"
    for i in "${names[@]}"
    do
        echo "$i"
        bsub -J "$i" -env "MATRIX_NAME=$i, NYSTROM=$j" < scripts/submit_matrix.sh;
    done
done
