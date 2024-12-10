#!/bin/bash

# 1. Read in geometry
# 2. Optimize to minimum
# 3. Calculate frequencies
# 4. Generate Wigner sampled structures
# 5. For each sampled structure, rotate and stretch bond
# 6. Compute energy for each 
# 7. Compile results

set -e

# Default values
xyz=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
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

# Check if --xyz is provided
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
echo "Reading XYZ: $xyz"

############################################
############################################
############################################

NAME="${xyz%.*}"
IMAGE="chem279/final"

mkdir -p $NAME
cd $NAME
rm -rf *
cp "../$xyz" ./

MOPAC_OPT_KEYWORDS="MNDO XYZ NOREOR OPT PRECISE"
MOPAC_FREQ_KEYWORDS="MNDO XYZ NOREOR LARGE DFORCE"
MOPAC_GRAD_KEYWORDS="MNDO XYZ NOREOR"

###########################################
# Opt
#####

OPT_INPUT="${NAME}_opt.mop"
OPT_OUTPUT="${NAME}_opt.out"
OPT_XYZ="${NAME}_opt.xyz"

echo "$MOPAC_OPT_KEYWORDS" > $OPT_INPUT
echo " " >> $OPT_INPUT
echo " " >> $OPT_INPUT
tail -n +3 $xyz >> $OPT_INPUT

docker run --rm -v "$(pwd)":/workdir $IMAGE mopac $OPT_INPUT
docker run --rm -v "$(pwd)":/workdir $IMAGE obabel $OPT_OUTPUT -o xyz -O $OPT_XYZ

##########################################
# Freq
######

FREQ_XYZ="${NAME}_freq.xyz"
cp $OPT_XYZ $FREQ_XYZ

FREQ_INPUT="${NAME}_freq.mop"
FREQ_OUTPUT="${NAME}_freq.out"
FREQ_MOLDEN="${NAME}_freq.molden"

echo "$MOPAC_FREQ_KEYWORDS" > $FREQ_INPUT
echo " " >> $FREQ_INPUT
echo " " >> $FREQ_INPUT
tail -n +3 $OPT_XYZ >> $FREQ_INPUT

docker run --rm -v "$(pwd)":/workdir $IMAGE mopac $FREQ_INPUT
docker run --rm -v "$(pwd)":/workdir $IMAGE obabel $FREQ_OUTPUT -o molden -O $FREQ_MOLDEN

##########################################
# Wigner and interpolation
##########################

cp ../wigner_sample.py ./
cp -r ../utils ./
docker run --rm -v "$(pwd)":/workdir $IMAGE python3 wigner_sample.py $FREQ_XYZ

##########################################
# Compute singlepoints
######################

for i in *sample*xyz
  do IND_NAME="${i%.*}"
  IND_INPUT="${IND_NAME}.mop"
  IND_OUTPUT="${IND_NAME}.out"

  echo "$MOPAC_GRAD_KEYWORDS" > $IND_INPUT
  echo " " >> $IND_INPUT
  echo " " >> $IND_INPUT

  docker run --rm -v "$(pwd)":/workdir $IMAGE mopac $IND_INPUT
done

#########################################
# Parse
#######

docker run --rm -v "$(pwd)":/workdir $IMAGE python3 ../parser.py

#########################################
# Done!
#######

echo "DONE"


