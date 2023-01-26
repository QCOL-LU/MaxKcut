#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_jogo_qubo)  #Change this runner_peel_decompose  ind_bqo _con_bqo runner_preprocess_weighted_amilo runner_preprocess_unweighted_amilo runner_preprocess_weighted_rpmilo runner_preprocess_unweighted_rpmilo runner_preprocess_weighted_pmilo runner_preprocess_unweighted_pmilo runner_preprocess_weighted_ind_bqo runner_preprocess_unweighted_ind_bqo runner_preprocess_weighted_con_bqo runner_preprocess_unweighted_con_bqo
NUM_PARTITIONS=(3 4)
METHODS=(5)                 # 0:MISDO 1:BQO 2:A-MILO 3:RP-MILO 4:P-MILO 5:C-QUBO 6:CR-QUBO 
PENALTYS=(0)            # 0:Tight 1:Adjusted-tight 2:Naive
path="../instances/jogo_qubo_tight/"


for folder_path in "$path"*; do
    for file in "$folder_path/"*; do
        for num_partitions in ${NUM_PARTITIONS[@]}; do
            for method in ${METHODS[@]}; do
                for penalty in ${PENALTYS[@]}; do
                    for runner in ${RUNNERS[@]}; do
                        COMMAND="python3.9 ../runners/${runner}.py ${file} $(basename "$file") ${num_partitions} ${method} ${penalty}"
                        eval $COMMAND
                    done
                done
            done
        done
    done
done