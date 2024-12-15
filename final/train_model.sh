#!/bin/bash

set -e

# Default values
data_dir=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --data_dir)
      data_dir="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Check if --data_dir is provided
if [[ -z "$data_dir" ]]; then
  echo "Error: --data_dir argument is required."
  exit 1
fi

# Check if the file exists
if [[ ! -d "$data_dir" ]]; then
  echo "Error: File '$data_dir' does not exist."
  exit 1
fi

# Output parsed arguments
echo "Using data from $data_dir"


IMAGE_TAG=chem279/final

mkdir -p "${data_dir}_model"
docker run --rm -v $(readlink -e "${data_dir}_model"):/workdir/results -v $(readlink -e $data_dir):/workdir/data_dir -v $(readlink -e model):/workdir/model $IMAGE_TAG python3 model/train.py

