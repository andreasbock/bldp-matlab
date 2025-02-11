#!/bin/sh
declare -a nystrom=("0")
#declare -a names=("Dubcova1" "gyro_m" "apache1" "crankseg_1" "cfd1" "thermal1")
#declare -a names=("bcsstm39" "jnlbrng1" "minsurfo" "gridgena" "finan512" "G2_circuit" "Dubcova2" "s3dkq4m2" "hood" "pwtk" "af_5_k101" "bundle_adj" "nd24k")
#declare -a names=("msc01050" "bcsstk11" "bcsstk26" "bcsstk24" "bcsstk16" "s2rmt3m1" "s3rmt3m1" "msc10848" "t2dah_e" "ct20stif" "shipsec8" "hood" "offshore")

#declare -a names=("scagr7-2r" "scrs8-2r" "scsd8-2r" "sctap1-2b" "south31" "PDE1" "lp_stocfor3" "lp_nug20" "lp_nug30" "mri1" "mri2" "nug08-3rd" "fome12" "car4" "dbir1" "dbir2" "dbic1" "e18" "ex3sta1" "mod2" "route" "world" "fxm3_16" "fxm4_6" "nsct" "Maragal_7" "Rucci1" "ESOC" "Maragal_8" "sls" "208bit" "tomographic1" "JP" "hood" "shipsec8" )


#declare -a names=("Dubcova1" "gyro_m" "apache1" "crankseg_1" "cfd1" "thermal1" "south31" "hood" "lp_stocfor3" "dbic1" "route" "sctap1-2b" "mri1" "e18" "ex3sta1" "fxm3_16" "fxm4_6" "sctap1-2b")
#declare -a names=("Dubcova1" "gyro_m" "apache1" "crankseg_1" "cfd1" "thermal1" "south31" "hood" "lp_stocfor3" "dbic1" "route" "sctap1-2b" "mri1" "e18" "ex3sta1" "fxm3_16" "fxm4_6")
#declare -a names=("offshore") # RUNNING
#declare -a names=("south31" "hood" "offshore" "lp_stocfor3" "mri1" "ex3sta1" "fxm3_16" "fxm4_6")
declare -a names=("Dubcova1" "gyro_m" "apache1" "crankseg_1" "cfd1" "thermal1")
#declare -a names=("crankseg_2")

for j in "${nystrom[@]}"
do
    echo "$j"
    for i in "${names[@]}"
    do
        echo "$i"
        bsub -J "$i-$j" -env "MATRIX_NAME=$i, NYSTROM=$j" < scripts/submit_matrix.sh;
    done
done
