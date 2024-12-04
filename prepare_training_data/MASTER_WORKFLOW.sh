#!/bin/bash

# MAIN WORKFLOW FOR PREPARING TRAINING DATA

# NEEDS:
#    1. Input molecule xyz (this is optimized to a minimum)
#    2. Number of non-equillibrium structures to generate


INPUT_XYZ=$1

MOPAC_OUTPUT=${INPUT_XYZ%.*}.out
MOLDEN_FILE=${INPUT_XYZ%.*}.molden

./mopac_opt_freq.sh $INPUT_XYZ
./obabel.sh -i mopout $MOPAC_OUTPUT -o molden -O $MOLDEN_FILE
./docker_python.sh wigner_sample.py -i $MOLDEN_FILE -n $2
./obabel.sh -m "wigner-${INPUT_XYZ%.*}-298.15.xyz" -o xyz -O sample_structure.xyz

for i in sample_structure*xyz;
do  ./mopac_grad.sh $i
done

./docker_python.sh parse.py 
