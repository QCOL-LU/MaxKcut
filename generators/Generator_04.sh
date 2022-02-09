#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_bqo_03 ) #Change this runner_pmilo_paper runner_rpmilo_paper runner_ind_bqo_paper runner_con_bqo_paper
NUM_PARTITIONS=(2 3 4)
path="instances/hojny_et_al_instances/I160/"


for file in "$path"*; do
    for runner in ${RUNNERS[@]}; do
        for num_partitions in ${NUM_PARTITIONS[@]}; do
            COMMAND="python3.9 ${runner}.py ${file} $(basename "$file") ${num_partitions}"
            eval $COMMAND
        done
    done
done