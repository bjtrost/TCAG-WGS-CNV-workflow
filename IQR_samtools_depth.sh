#!/usr/bin/env bash

# To calculate IQR based only on chromosome 1, use the following command:
# IQR_samtools_depth.sh mysample.bam 1

# To calculate IQR based on the entire genome, use the following command:
# IQR_samtools_depth.sh mysample.bam

BAM_FILE=$1
CHR=$2

if [[ $CHR != "" ]]; then
    samtools depth -r $CHR $BAM_FILE | IQR_samtools_depth.py
else
    samtools depth $BAM_FILE | IQR_samtools_depth.py
fi
