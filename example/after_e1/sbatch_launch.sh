#!/bin/bash

# Launch script for epoch {epoch_num}

set -e

for ((i=1; i<17; i++)); do
  cd rep$(printf "%02d" $i)
  gmx mdrun  -deffnm mdrun
  cd ..
done
