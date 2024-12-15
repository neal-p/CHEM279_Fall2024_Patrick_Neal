#!/bin/bash

cd ../build && cmake .. && make hw1_3
cd -

mkdir -p my_output
mkdir -p my_output/SD_with_line_search

for f in sample_input/SD_with_line_search/*txt
do
  echo "Running ${f}"
  output=$(basename "${f}")
  ./hw1_3 "${f}" > "my_output/SD_with_line_search/${output}"
done
