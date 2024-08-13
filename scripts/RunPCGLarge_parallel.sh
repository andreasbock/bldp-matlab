#!/bin/sh
declare -a nystrom=("0" "1")
declare -a names=("Dubcova1" "gyro_m" "apache1" "crankseg_1" "cfd1" "thermal1")

for j in "${nystrom[@]}"
do
    echo "$j"
    for i in "${names[@]}"
    do
        echo "$i"
        bsub -J "$i-$j" -env "MATRIX_NAME=$i, NYSTROM=$j" < scripts/submit_matrix.sh;
    done
done
