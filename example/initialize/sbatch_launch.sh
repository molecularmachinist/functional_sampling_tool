#!/bin/bash

## This line will not be copied

## The {epoch_num} below will be substituted with the epoch number for each epoch
# Launch script for epoch {epoch_num}

set -e

for ((i=1; i<17; i++)); do
  cd rep$(printf "%02d" $i)
  gmx mdrun  -deffnm mdrun
  cd ..
done
