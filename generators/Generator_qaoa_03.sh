#!/bin/bash

PWD=`pwd`
RUNNERS=( ../runners/runner_qaoa_rqubo_03 ../runners/runner_qaoa_rqubo_02) #Change this runner_qaoa_rqubo_02
NUM_PARTITIONS=(2)
PENALTY_COEFF_INCREASE=(1.0)
NOISES=(0.00)



for runner in ${RUNNERS[@]}; do
    for seed in in {1..1}; do
        for num_vertices in {17..33}; do
            for num_partitions in ${NUM_PARTITIONS[@]}; do
                for noise in ${NOISES[@]}; do
                    for penalty_coeff_increase in ${PENALTY_COEFF_INCREASE[@]}; do
                        COMMAND="python3.9 ${runner}.py ${seed} ${num_vertices} ${num_partitions} ${penalty_coeff_increase} ${noise}"
                        eval $COMMAND
                    done
                done
            done
        done
    done
done
