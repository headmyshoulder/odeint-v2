#!/bin/zsh

export LC_NUMERIC=en_US.UTF-8
declare -A times

bench="bin/clang-linux-3.2/release/threading-multi/osc_chain_1d"

export OMP_SCHEDULE=static
repeat=2

function run {
    n=$1
    steps=$2
    printf "# n=$n steps=$steps repeat=$repeat\n"
    printf '"block"\t"med"\t"mul"'
    for block in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ; do
        printf '\n%d' $block
        med=$(mpirun -np $block $bench $n $steps $repeat | tail -1 | awk '{print $4}')
        times[$build-$block]=$med
        speedup=$((${times[$build-1]}/$med))
        printf '\t%f\t%f' $med $speedup
    done
    printf '\n\n\n'
}

run 4096 1024 | tee osc_chain_speedup-short.dat
run 4194304 1 | tee osc_chain_speedup-long.dat
