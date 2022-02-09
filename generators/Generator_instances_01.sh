#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_preprocess_weighted_ind_bqo runner_preprocess_weighted_con_bqo runner_preprocess_weighted_amilo runner_preprocess_weighted_rpmilo runner_preprocess_weighted_pmilo) #Change this runner_peel_decompose  ind_bqo _con_bqo runner_preprocess_weighted_amilo runner_preprocess_unweighted_amilo runner_preprocess_weighted_rpmilo runner_preprocess_unweighted_rpmilo runner_preprocess_weighted_pmilo runner_preprocess_unweighted_pmilo runner_preprocess_weighted_ind_bqo runner_preprocess_unweighted_ind_bqo runner_preprocess_weighted_con_bqo runner_preprocess_unweighted_con_bqo
NUM_PARTITIONS=(2 3 4)
# path="instances/districting_2020/cities/"
# path="instances/wang_hijazi_instances/band/"
path="instances/wang_hijazi_instances/spinglass/"


for file in "$path"*; do
    for runner in ${RUNNERS[@]}; do
        for num_partitions in ${NUM_PARTITIONS[@]}; do
            COMMAND="python3.9 ${runner}.py ${file} $(basename "$file") ${num_partitions}"
            eval $COMMAND
        done
    done
done