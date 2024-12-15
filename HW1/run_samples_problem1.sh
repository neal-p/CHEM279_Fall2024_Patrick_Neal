#!/bin/bash

cd ../build && make hw1_1
cd -

mkdir -p my_output
mkdir -p my_output/Energy/

for f in sample_input/Energy/*txt
do
  echo "Running ${f}"
  output=$(basename "${f}")
  ./hw1_1 $f > my_output/Energy/$output
done
