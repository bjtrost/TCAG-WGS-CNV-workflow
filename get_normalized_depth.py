#!/usr/bin/env python3

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
parser.add_argument("-o", "--generate-outer-flank-ratio", action="store_true", default=False)
args = parser.parse_args()
#####################################

chrom, start, end = CNVworkflowlib.parse_region(args.region.replace(",", ""))
chrom = chrom.replace("chr", "")
avg_chromosome_depth = CNVworkflowlib.dict_from_file(args.depth_filename + ".chrinfo", header=True, value_type="float", value_col=2)

depth_lists = CNVworkflowlib.get_depth_lists(chrom, start, end, args.depth_filename, compute_left_left_right_right=args.generate_outer_flank_ratio)


if args.generate_outer_flank_ratio:
    both_outer_flanks_depths = depth_lists["left_left_flank_depths"] + depth_lists["right_right_flank_depths"]
    outer_flanks_mean = CNVworkflowlib.mean(both_outer_flanks_depths)
    normalized_outer_flanking_depth = outer_flanks_mean / avg_chromosome_depth[chrom]

both_flanks_depths = depth_lists["left_flank_depths"] + depth_lists["right_flank_depths"]
flanks_mean = CNVworkflowlib.mean(both_flanks_depths)
CNV_mean = CNVworkflowlib.mean(depth_lists["CNV_depths"])
CNV_stdev = CNVworkflowlib.stdev(depth_lists["CNV_depths"])
flanks_stdev = CNVworkflowlib.stdev(both_flanks_depths)


try:
    normalized_CNV_stdev = CNV_stdev / CNV_mean
except:
    normalized_CNV_stdev = 0

try:
    normalized_flanks_stdev = flanks_stdev / flanks_mean
except:
    normalized_flanks_stdev = 0

normalized_CNV_depth = CNV_mean / avg_chromosome_depth[chrom]
normalized_flanking_depth = flanks_mean / avg_chromosome_depth[chrom]

# Meanings of output columns:
# column 1: mean read depth of the CNV
# column 2: mean read depth of the two flanking regions (each of the same size as the CNV)
# column 3: ratio between the mean read depth of the CNV and the mean read depth of the chromosome (col1/col6)
# column 4: ratio between the mean read depth of the flanking regions and the mean read depth of the chromosome (col2/col6)
# column 5: ratio between the mean read depth of the CNV and the mean read depth of the flanking regions (col1/col2)
# column 6: mean read depth of the chromosome


if args.generate_outer_flank_ratio:
    try:
        ratio = normalized_CNV_depth / normalized_flanking_depth
        outer_ratio = normalized_CNV_depth / normalized_outer_flanking_depth
        print("{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(CNV_mean, flanks_mean, normalized_CNV_depth, normalized_flanking_depth, ratio, avg_chromosome_depth[chrom], outer_flanks_mean, outer_ratio))
    except ZeroDivisionError:
        print("{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t{:.3f}\t{}".format(CNV_mean, flanks_mean, normalized_CNV_depth, normalized_flanking_depth, "N/A", avg_chromosome_depth[chrom], outer_flanks_mean, "N/A"))
else:
    try:
        ratio = normalized_CNV_depth / normalized_flanking_depth
        print("{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(CNV_mean, flanks_mean, normalized_CNV_depth, normalized_flanking_depth, ratio, avg_chromosome_depth[chrom], CNV_stdev, normalized_CNV_stdev, flanks_stdev, normalized_flanks_stdev))
    except ZeroDivisionError:
        print("{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(CNV_mean, flanks_mean, normalized_CNV_depth, normalized_flanking_depth, "N/A", avg_chromosome_depth[chrom], CNV_stdev, normalized_CNV_stdev, flanks_stdev, normalized_flanks_stdev))
