#!/bin/sh
declare -a nystrom=("0" "1")
declare -a names=("Pres_Poisson" "bmwcra_1" "cfd1" "smt" "apache1" "thermal1" "crankseg_1" "crankseg_2" "bcsstk17" "bcsstk18" "consph")

for j in "${nystrom[@]}"
do
    echo "$j"
    for i in "${names[@]}"
    do
        echo "$i"
        bsub -J "NYS=$j_$i" -env "MATRIX_NAME=$i, NYSTROM=$j" < scripts/submit_matrix.sh;
    done
done
