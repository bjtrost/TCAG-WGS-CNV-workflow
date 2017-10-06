#!/usr/bin/env bash

BAM_file=$1
region=$2  # In format "chr1:30000-40000"

if [ ! -e $BAM_file.depth ]; then
    samtools depth $BAM_file > $BAM_file.depth
fi

if [ ! -e $BAM_file.depth.idx ]; then
    index_samtools_depth.py $BAM_file.depth
fi

get_normalized_depth.py $BAM_file.depth $region
