#!/usr/bin/env python3

# Compare a list of CNV calls to a list of repetitive and low-complexity regions (RLCRs)
# Determine how many bases of each CNV overlap with these unclean regions, what percentage overlap,
# and what types of regions overlap

# Usage: compare_with_RLCR_definition.py RLCR_definition_filename CNV_filename1 [CNV_filename2, ...]

import CNVworkflowlib
import argparse
import copy
from intervaltree import Interval, IntervalTree

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.set_defaults(output_percentage_overlapping=False)
parser.set_defaults(always_output_overlapping_regions=False)

parser.add_argument("RLCR_definition_filename", type=str)
parser.add_argument("-p", "--output-percentage-overlapping", action="store_true")
parser.add_argument("-o", "--always-output-overlapping-regions", action="store_true")
parser.add_argument("CNV_filenames", type=str, default=None, nargs="+")
args = parser.parse_args()
#####################################

chromosomes = CNVworkflowlib.get_chromosomes()
trees = {}
RLCR_definition_file = open(args.RLCR_definition_filename)

# Make a new interval tree for each chromosome
for chrom in chromosomes:
    trees[chrom] = IntervalTree()

# Read RLCR definition file into interval tree data structure
for line in RLCR_definition_file:
    chrom, start, end, repeat_type = line.rstrip().split("\t")

    if chrom == "Chromosome":
        continue

    chrom = chrom.replace("chr", "")
    start = int(start)
    end = int(end)
    trees[chrom][start:end+1] = repeat_type # Insert this interval into the tree
RLCR_definition_file.close()

# Make a deep copy of the tree so that we can merge overlapping intervals in the copied one
trees_merged = copy.deepcopy(trees)
for chrom in chromosomes:
    trees_merged[chrom].merge_overlaps()

for CNV_filename in args.CNV_filenames:
    CNV_file_out = open(CNV_filename + ".RLCR", "w")

    CNVs, header_line = CNVworkflowlib.read_common_format(CNV_filename)


    if args.output_percentage_overlapping:
        CNV_file_out.write("{}\t{}\t{}\t{}\n".format(header_line, "Number of bases overlapping RLCRs", "Percentage of bases overlapping", "Regions overlapping RLCRs")) # Print header as in original file, plus additional headers related to overlap with RLCRs
    else:
        CNV_file_out.write("{}\t{}\t{}\n".format(header_line, "Number of bases overlapping RLCRs", "Regions overlapping RLCRs")) # Print header as in original file, plus additional headers related to overlap with RLCRs

    for CNV in CNVs:
        # Calculate the total overlap
        overlap = 0
        for i in trees_merged[CNV["chr"]][CNV["start"]:CNV["end"]+1]: # Must add 1 to end of region due to the way intervaltree handles coordinates
            overlap += CNVworkflowlib.get_overlap([CNV["start"], CNV["end"]], [i.begin, i.end-1]) # Need to subtract 1 from i.end in order to get the correct overlap (since the IntervalTree package does not include the end in the interval)

        regions_overlapped = ":".join([i.data for i in trees[CNV["chr"]][CNV["start"]:CNV["end"]+1]]) # Get the names of the overlapping regions. Must add 1 to end of region due to the way intervaltree handles coordinates

        overlap_percentage = overlap / (CNV["end"] - CNV["start"] + 1) * 100

        if overlap_percentage < 70 and not args.always_output_overlapping_regions:
            regions_overlapped = "SIMPLE"

        if args.output_percentage_overlapping:
            CNV_file_out.write("{}\t{:d}\t{:.1f}%\t{}\n".format(CNV["full_line"], overlap, overlap_percentage, regions_overlapped))
        else:
            CNV_file_out.write("{}\t{:d}\t{}\n".format(CNV["full_line"], overlap, regions_overlapped))
    CNV_file_out.close()
