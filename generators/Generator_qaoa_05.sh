#!/bin/bash

PWD=`pwd`
RUNNERS=(../runners/runner_qaoa_rqubo_05) #Change this 
NUM_PARTITIONS=(3)
PENALTY_MULTIPLIER=(1.0) # 2.0
NUM_VERTICES=(8)
GRAPH_DENSITY=(0.8) #0.2 0.4 0.6
NEG_EDGE_PERCENTAGE=(0.0) # 0.4 0.8




for runner in ${RUNNERS[@]}; do  
    for seed in {1..1}; do
        for num_partitions in ${NUM_PARTITIONS[@]}; do
            for num_vertices in ${NUM_VERTICES[@]}; do
                for graph_density in ${GRAPH_DENSITY[@]}; do
                    for neg_edge_percentage in ${NEG_EDGE_PERCENTAGE[@]}; do
                        for penalty_multiplier in ${PENALTY_MULTIPLIER[@]}; do
                            COMMAND="python3.9 ${runner}.py ${seed} ${num_vertices} ${num_partitions} ${graph_density} ${neg_edge_percentage} ${penalty_multiplier} "
                            eval $COMMAND
                        done
                    done
                done
            done
        done
    done
done
