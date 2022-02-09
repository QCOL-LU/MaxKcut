#!/bin/bash

PWD=`pwd`
RUNNERS=(runner_bqo_01) #Change this


for seed in {1..1}; do
    for runner in ${RUNNERS[@]}; do
        for numv in {30..40..10}; do #Change this
            for partition in {3..5..2}; do
                COMMAND="python3.9 ${runner}.py ${seed} ${numv} ${partition}"
                eval $COMMAND
            done
        done
    done
done