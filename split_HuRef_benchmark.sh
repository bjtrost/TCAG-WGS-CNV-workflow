#!/usr/bin/env bash

# This script takes the HuRef benchmark file and splits it into separate files for each method. These separate files
# can then be used as input to CNV_overlap.py

#Usage:
# split_HuRef_benchmark.sh File_S1.txt

file=$1 # File containing HuRef benchmark

for method in Affymetrix6.0 Agilent24M Illumina1M NimbleGen42M SplitRead assembly.comparison cg.cnv cg.mp
do
    printf "Chromosome\tStart position\tEnd position\tSize\tType\tMethod of CNV detection\n" > $file.$method
    awk -v m=$method -F $'\t' '$6 == m' $file >> $file.$method
done
