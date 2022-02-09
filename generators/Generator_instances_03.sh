#!/bin/bash

PWD=`pwd`
RUNNERS=(../runners/runner_preprocess) #Change this runner_peel runner_fold runner_decompose 
NUM_PARTITIONS=(2 )
path="../instances/cities/"


for file in "$path"*; do
    for runner in ${RUNNERS[@]}; do
        for num_partitions in ${NUM_PARTITIONS[@]}; do
            COMMAND="python3.9 ${runner}.py ${file} $(basename "$file") ${num_partitions}"
            eval $COMMAND
        done
    done
done