#!/usr/bin/env python

import argparse
import CNVworkflowlib

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("samtools_depth_filename", type=str, nargs="?", default="-")
args = parser.parse_args()
#####################################

def find_cumulative_val(v, total, percent):
    cumulative_sum = 0

    for i in sorted(list(v.keys())):
	    cumulative_sum = cumulative_sum + v[i]

	    if cumulative_sum > total * (percent / 100):
		    quartile_num = i
		    break

    return(quartile_num)


samtools_depth_file = CNVworkflowlib.file_or_stdin(args.samtools_depth_filename)

depth_count = {}
total = 0

for line in samtools_depth_file:
    depth = int(line.rstrip().split("\t")[2])

    if depth not in depth_count:
        depth_count[depth] = 0

    depth_count[depth] += 1
    total += 1


upper_quartile_num = find_cumulative_val(depth_count, total, 75)
lower_quartile_num = find_cumulative_val(depth_count, total, 25)

IQR = upper_quartile_num - lower_quartile_num

print(IQR)
