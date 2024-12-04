#!/bin/bash

IMAGE_TAG=chem279/final
MOUNT_DIR=$(pwd)

docker run --rm -v $MOUNT_DIR:/workdir $IMAGE_TAG obabel $@
