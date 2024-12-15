#!/bin/bash

cd ../build && make hw5
cd -

mkdir -p my_output/

for f in sample_input/*txt
do
  echo "Running ${f}"
  output=$(basename "${f}")
  ./hw5 $f basis/basis_set.txt > my_output/$output
done
