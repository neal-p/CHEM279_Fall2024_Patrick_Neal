#!/bin/bash

cd ../build && make hw3
cd -

mkdir -p my_output/

for f in sample_input/*txt
do
  echo "Running ${f}"
  output=$(basename "${f}")
  ./hw3 $f basis_set.txt > my_output/$output
done
