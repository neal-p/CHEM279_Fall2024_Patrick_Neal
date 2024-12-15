#!/bin/bash

cd ../build && make hw2_1
cd -

mkdir -p my_output
mkdir -p my_output/numerical/

for f in sample_input/numerical/*txt
do
  echo "Running ${f}"
  output=$(basename "${f}")
  ./hw2_1 $f > my_output/numerical/$output
done
