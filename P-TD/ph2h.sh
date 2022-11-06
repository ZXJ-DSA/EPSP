#!/bin/bash
source_path=$1
algo=$2
task=$3
update=$4
for dataset in W
do
for i in 16 32 64 128 256
do
   ./build/PartiH2H $source_path $dataset $algo $i $task $update
done
done
