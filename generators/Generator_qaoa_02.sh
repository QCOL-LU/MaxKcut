#!/bin/bash

PWD=`pwd`
RUNNERS=(../runners/runner_qaoa_rqubo_01) #Change this 
NUM_PARTITIONS=(3 4 2)
PENALTY_COEFF_INCREASE=(1.0 2.0)
NUM_VERTICES=(7 8 6)
NOISES=(0.0)



for runner in ${RUNNERS[@]}; do  
    for seed in {1..1}; do
        for num_partitions in ${NUM_PARTITIONS[@]}; do
            for num_vertices in ${NUM_VERTICES[@]}; do
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
