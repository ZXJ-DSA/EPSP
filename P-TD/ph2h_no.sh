#!/bin/bash
source_path=$1
algo=$2
task=$3
method=$4
for dataset in NY FLA
do
   ./build/$method $source_path $dataset $algo 64 $task 0
done
