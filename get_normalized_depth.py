#!/usr/bin/env python

import argparse
import CNVworkflowlib

description = '''
Takes a samtools depth file and a genomic interval in the format chr:start-end,
and outputs the average depth in that interval, as well as the normalized depth
(the average depth of the interval divided by the average depth in the chromosome
comprising that interval)
'''

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("depth_filename", type=str, help="filename of samtools depth file", default=None)
parser.add_argument("region", type=str, help="Region to find the depth of", default=None)
args = parser.parse_args()
#####################################

chr, start, end = CNVworkflowlib.parse_region(args.region.replace(",", ""))
left_flank_depths, CNV_depths, right_flank_depths = CNVworkflowlib.get_depth_lists(chr, start, end, args.depth_filename)
both_flanks_depths = left_flank_depths + right_flank_depths

flanks_mean = CNVworkflowlib.mean(both_flanks_depths)
flanks_stdev = CNVworkflowlib.stdev(both_flanks_depths)

CNV_mean = CNVworkflowlib.mean(CNV_depths)
CNV_stdev = CNVworkflowlib.stdev(CNV_depths)

avg_chromosome_depth = CNVworkflowlib.dict_from_file(args.depth_filename + ".chrinfo", header=True, value_type="float", value_col=2)

normalized_CNV_depth = CNV_mean / avg_chromosome_depth[chr]
normalized_flanking_depth = flanks_mean / avg_chromosome_depth[chr]

# Meanings of output columns:
# column 1: mean read depth of the CNV
# column 2: mean read depth of the two flanking regions (each of the same size as the CNV)
# column 3: ratio between the mean read depth of the CNV and the mean read depth of the chromosome (col1/col6)
# column 4: ratio between the mean read depth of the flanking regions and the mean read depth of the chromosome (col2/col6)
# column 5: ratio between the mean read depth of the CNV and the mean read depth of the flanking regions (col1/col2)
# column 6: mean read depth of the chromosome

try:
    ratio = normalized_CNV_depth / normalized_flanking_depth
    print("{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(CNV_mean, flanks_mean, normalized_CNV_depth, normalized_flanking_depth, ratio, avg_chromosome_depth[chr]))
except ZeroDivisionError:
    print("{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{}\t{:.3f}".format(CNV_mean, flanks_mean, normalized_CNV_depth, normalized_flanking_depth, "N/A", avg_chromosome_depth[chr]))
