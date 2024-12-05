#!/bin/bash

IMAGE_TAG=chem279/final
MOUNT_DIR=$(pwd)

docker run --rm -v $MOUNT_DIR:/workdir  -v  /home/cjpeh/msse/CHEM279_Final/prepare_training_data:/workdir/prepare_training_data $IMAGE_TAG python3 $@

