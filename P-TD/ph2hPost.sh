#!/bin/bash
source_path=$1
algo=$2
i=$3
for dataset in NY FLA
do
   ./build/PartiH2HPost $source_path $dataset $algo $i
done
