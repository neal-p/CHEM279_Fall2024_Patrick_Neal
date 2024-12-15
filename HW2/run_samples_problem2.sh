#!/bin/bash

cd ../build && make hw2_2
cd -

mkdir -p my_output
mkdir -p my_output/analytical/

for f in sample_input/analytical/*txt
do
  echo "Running ${f}"
  output=$(basename "${f}")
  ./hw2_2 $f > my_output/analytical/$output
done
