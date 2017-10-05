#!/usr/bin/env python3

import CNVworkflowlib
import argparse


####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("algorithm_CNV_filename", type=str)
parser.add_argument("benchmark_CNV_filename", type=str)
args = parser.parse_args()
#####################################

algorithm_CNVs, header_line = CNVworkflowlib.read_common_format(args.algorithm_CNV_filename)

benchmark_CNVs, _ = CNVworkflowlib.read_benchmark(args.benchmark_CNV_filename)

print("{}\t{}\t{}".format(header_line, "List of overlap lengths", "List of overlapping benchmark CNVs"))

for i in range(0, len(algorithm_CNVs)):
    print(algorithm_CNVs[i]["full_line"], end="")
    overlap_list = []
    match_list = []

    for j in range(0, len(benchmark_CNVs)):
        if CNVworkflowlib.fifty_percent_reciprocal_overlap(algorithm_CNVs[i], benchmark_CNVs[j]):
            overlap_list.append(str(CNVworkflowlib.get_overlap_region(algorithm_CNVs[i], benchmark_CNVs[j]))) # List of overlap lengths
            match_list.append("{}:{}:{}:{}:{}:{}".format(benchmark_CNVs[j]["chr"], benchmark_CNVs[j]["start"], benchmark_CNVs[j]["end"], benchmark_CNVs[j]["size"], benchmark_CNVs[j]["type"], benchmark_CNVs[j]["caller"])) # List of benchmark CNVs that overlap this predicted CNV

    if overlap_list:
        print("\t{}\t{}".format("|".join(overlap_list), "|".join(match_list)))
    else:
        print("\t{}\t{}".format("0", "NotApplicable"))
