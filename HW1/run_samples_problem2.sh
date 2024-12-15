#!/bin/bash

cd ../build && cmake .. && make hw1_2
cd -

mkdir -p my_output
mkdir -p my_output/Force/

for f in sample_input/Force/*txt
do
  echo "Running ${f}"
  output=$(basename "${f}")
  ./hw1_2 "${f}" "my_output/Force/${output%.*}.csv" > "my_output/Force/${output}"
done
