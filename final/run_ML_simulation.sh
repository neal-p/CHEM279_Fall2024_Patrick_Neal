#!/bin/bash

set -e

# Default values
model_dir=""
xyz=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --model_dir)
      model_dir="$2"
      shift 2
      ;;
    --xyz)
      xyz="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Check if --model_dir is provided
if [[ -z "$model_dir" ]]; then
  echo "Error: --model_dir argument is required."
  exit 1
fi

# Check if the file exists
if [[ ! -d "$model_dir" ]]; then
  echo "Error: File '$model_dir' does not exist."
  exit 1
fi

# Check if --model_dir is provided
if [[ -z "$xyz" ]]; then
  echo "Error: --xyz argument is required."
  exit 1
fi

# Check if the file exists
if [[ ! -f "$xyz" ]]; then
  echo "Error: File '$xyz' does not exist."
  exit 1
fi



# Output parsed arguments
echo "Using model from $model_dir with $xyz as starting geometry"


IMAGE_TAG=chem279/final

mkdir -p "${model_dir}_ML_simulation"
cp $xyz "${model_dir}_ML_simulation"

docker run --rm -v $(readlink -e utils):/workdir/utils -v $(readlink -e "${model_dir}_ML_simulation"):/workdir/results -v $(readlink -e $model_dir):/workdir/model_dir -v $(readlink -e model):/workdir/model $IMAGE_TAG python3 model/simulation.py "results/$xyz"

