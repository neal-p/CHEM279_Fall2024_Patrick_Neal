#!/bin/bash

cd ../build && make hw4 scan
cd -

mkdir -p my_output/

for f in sample_input/*txt
do
  echo "Running ${f}"
  output=$(basename "${f}")
  ./hw4 $f basis/basis_set.txt > my_output/$output
done
