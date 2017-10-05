#!/usr/bin/env python

import BTlib
import argparse

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("CNV_file", type=str, help="CNV file", default=None)
args = parser.parse_args()
#####################################

f = BTlib.file_or_stdin(args.CNV_file)

size_bins = ["Total", "[1000,5000)", "[5000,10000)", "[10000,100000)", "[100000,1000000)", "[1000000,...)"]

# "CONFIRM" means found in the benchmark
# "RAW" means all calls (found in the benchmark or not)
# "SIMPLE" means less than 70% overlap with RLCRs
# "REPEAT" means >= 70%
# "COUNT1" means exactly one confirming benchmark method
# "COUNT2" means two or more confirming benchmark methods
# "DEL" means deletion
# "DUP" means duplication
categories = [
"CONFIRM",
"CONFIRM_COUNT1",
"CONFIRM_COUNT2",
"CONFIRM_DEL",
"CONFIRM_DEL_COUNT1",
"CONFIRM_DEL_COUNT2",
"CONFIRM_DEL_REPEAT",
"CONFIRM_DEL_REPEAT_COUNT1",
"CONFIRM_DEL_REPEAT_COUNT2",
"CONFIRM_DEL_SIMPLE",
"CONFIRM_DEL_SIMPLE_COUNT1",
"CONFIRM_DEL_SIMPLE_COUNT2",
"CONFIRM_DUP",
"CONFIRM_DUP_COUNT1",
"CONFIRM_DUP_COUNT2",
"CONFIRM_DUP_REPEAT",
"CONFIRM_DUP_REPEAT_COUNT1",
"CONFIRM_DUP_REPEAT_COUNT2",
"CONFIRM_DUP_SIMPLE",
"CONFIRM_DUP_SIMPLE_COUNT1",
"CONFIRM_DUP_SIMPLE_COUNT2",
"CONFIRM_REPEAT",
"CONFIRM_REPEAT_COUNT1",
"CONFIRM_REPEAT_COUNT2",
"CONFIRM_SIMPLE",
"CONFIRM_SIMPLE_COUNT1",
"CONFIRM_SIMPLE_COUNT2",
"RAW",
"RAW_DEL",
"RAW_DEL_COUNT1",
"RAW_DEL_COUNT2",
"RAW_DEL_REPEAT",
"RAW_DEL_REPEAT_COUNT1",
"RAW_DEL_REPEAT_COUNT2",
"RAW_DEL_SIMPLE",
"RAW_DEL_SIMPLE_COUNT1",
"RAW_DEL_SIMPLE_COUNT2",
"RAW_DUP",
"RAW_DUP_COUNT1",
"RAW_DUP_COUNT2",
"RAW_DUP_REPEAT",
"RAW_DUP_REPEAT_COUNT1",
"RAW_DUP_REPEAT_COUNT2",
"RAW_DUP_SIMPLE",
"RAW_DUP_SIMPLE_COUNT1",
"RAW_DUP_SIMPLE_COUNT2",
"RAW_REPEAT",
"RAW_SIMPLE"
]


counts = {} # Initialize dictionary
for size_bin in size_bins:
    counts[size_bin] = {}
    for category in categories:
        counts[size_bin][category] = 0

f.readline() # Discard header line

for line in f:
    l = line.rstrip()
    fields = l.split("\t")
    size = int(float(fields[3]))
    CNV_type = fields[4]

    if size < 1000:
        continue

    if fields[8] == "0":
        num_confirming = 0
    else:
        num_confirming = int(fields[8].count("|")) + 1

    if "NotApplicable" in l and num_confirming > 0:
        print(fields[8] + "\t" + str(num_confirming))

    repeat = fields[-1]
    size_bin = BTlib.get_size_bin(size)

    if num_confirming > 0:

        counts[size_bin]["CONFIRM"] += 1
        counts["Total"]["CONFIRM"] += 1

    if num_confirming == 1:
        counts[size_bin]["CONFIRM_COUNT1"] += 1
        counts["Total"]["CONFIRM_COUNT1"] += 1

    if num_confirming > 1:
        counts[size_bin]["CONFIRM_COUNT2"] += 1
        counts["Total"]["CONFIRM_COUNT2"] += 1

    if num_confirming > 0 and CNV_type == "DEL":
        counts[size_bin]["CONFIRM_DEL"] += 1
        counts["Total"]["CONFIRM_DEL"] += 1

    if num_confirming == 1 and CNV_type == "DEL":
        counts[size_bin]["CONFIRM_DEL_COUNT1"] += 1
        counts["Total"]["CONFIRM_DEL_COUNT1"] += 1

    if num_confirming > 1 and CNV_type == "DEL":
        counts[size_bin]["CONFIRM_DEL_COUNT2"] += 1
        counts["Total"]["CONFIRM_DEL_COUNT2"] += 1

    if num_confirming > 0 and CNV_type == "DEL" and repeat != "SIMPLE":
        counts[size_bin]["CONFIRM_DEL_REPEAT"] += 1
        counts["Total"]["CONFIRM_DEL_REPEAT"] += 1

    if num_confirming == 1 and CNV_type == "DEL" and repeat != "SIMPLE":
        counts[size_bin]["CONFIRM_DEL_REPEAT_COUNT1"] += 1
        counts["Total"]["CONFIRM_DEL_REPEAT_COUNT1"] += 1

    if num_confirming > 1 and CNV_type == "DEL" and repeat != "SIMPLE":
        counts[size_bin]["CONFIRM_DEL_REPEAT_COUNT2"] += 1
        counts["Total"]["CONFIRM_DEL_REPEAT_COUNT2"] += 1

    if num_confirming > 0 and CNV_type == "DEL" and repeat == "SIMPLE":
        counts[size_bin]["CONFIRM_DEL_SIMPLE"] += 1
        counts["Total"]["CONFIRM_DEL_SIMPLE"] += 1

    if num_confirming == 1 and CNV_type == "DEL" and repeat == "SIMPLE":
        counts[size_bin]["CONFIRM_DEL_SIMPLE_COUNT1"] += 1
        counts["Total"]["CONFIRM_DEL_SIMPLE_COUNT1"] += 1

    if num_confirming > 1 and CNV_type == "DEL" and repeat == "SIMPLE":
        counts[size_bin]["CONFIRM_DEL_SIMPLE_COUNT2"] += 1
        counts["Total"]["CONFIRM_DEL_SIMPLE_COUNT2"] += 1

    if num_confirming > 0 and CNV_type == "DUP":
        counts[size_bin]["CONFIRM_DUP"] += 1
        counts["Total"]["CONFIRM_DUP"] += 1

    if num_confirming == 1 and CNV_type == "DUP":
        counts[size_bin]["CONFIRM_DUP_COUNT1"] += 1
        counts["Total"]["CONFIRM_DUP_COUNT1"] += 1

    if num_confirming > 1 and CNV_type == "DUP":
        counts[size_bin]["CONFIRM_DUP_COUNT2"] += 1
        counts["Total"]["CONFIRM_DUP_COUNT2"] += 1

    if num_confirming > 0 and CNV_type == "DUP" and repeat != "SIMPLE":
        counts[size_bin]["CONFIRM_DUP_REPEAT"] += 1
        counts["Total"]["CONFIRM_DUP_REPEAT"] += 1

    if num_confirming == 1 and CNV_type == "DUP" and repeat != "SIMPLE":
        counts[size_bin]["CONFIRM_DUP_REPEAT_COUNT1"] += 1
        counts["Total"]["CONFIRM_DUP_REPEAT_COUNT1"] += 1

    if num_confirming > 1 and CNV_type == "DUP" and repeat != "SIMPLE":
        counts[size_bin]["CONFIRM_DUP_REPEAT_COUNT2"] += 1
        counts["Total"]["CONFIRM_DUP_REPEAT_COUNT2"] += 1

    if num_confirming > 0 and CNV_type == "DUP" and repeat == "SIMPLE":
        counts[size_bin]["CONFIRM_DUP_SIMPLE"] += 1
        counts["Total"]["CONFIRM_DUP_SIMPLE"] += 1

    if num_confirming == 1 and CNV_type == "DUP" and repeat == "SIMPLE":
        counts[size_bin]["CONFIRM_DUP_SIMPLE_COUNT1"] += 1
        counts["Total"]["CONFIRM_DUP_SIMPLE_COUNT1"] += 1

    if num_confirming > 1 and CNV_type == "DUP" and repeat == "SIMPLE":
        counts[size_bin]["CONFIRM_DUP_SIMPLE_COUNT2"] += 1
        counts["Total"]["CONFIRM_DUP_SIMPLE_COUNT2"] += 1

    if num_confirming > 0 and repeat != "SIMPLE":
        counts[size_bin]["CONFIRM_REPEAT"] += 1
        counts["Total"]["CONFIRM_REPEAT"] += 1

    if num_confirming == 1 and repeat != "SIMPLE":
        counts[size_bin]["CONFIRM_REPEAT_COUNT1"] += 1
        counts["Total"]["CONFIRM_REPEAT_COUNT1"] += 1

    if num_confirming > 1 and repeat != "SIMPLE":
        counts[size_bin]["CONFIRM_REPEAT_COUNT2"] += 1
        counts["Total"]["CONFIRM_REPEAT_COUNT2"] += 1

    if num_confirming > 0 and repeat == "SIMPLE":
        counts[size_bin]["CONFIRM_SIMPLE"] += 1
        counts["Total"]["CONFIRM_SIMPLE"] += 1

    if num_confirming == 1 and repeat == "SIMPLE":
        counts[size_bin]["CONFIRM_SIMPLE_COUNT1"] += 1
        counts["Total"]["CONFIRM_SIMPLE_COUNT1"] += 1

    if num_confirming > 1 and repeat == "SIMPLE":
        counts[size_bin]["CONFIRM_SIMPLE_COUNT2"] += 1
        counts["Total"]["CONFIRM_SIMPLE_COUNT2"] += 1

    counts[size_bin]["RAW"] += 1
    counts["Total"]["RAW"] += 1

    if repeat != "SIMPLE":
        counts[size_bin]["RAW_REPEAT"] += 1
        counts["Total"]["RAW_REPEAT"] += 1

    if repeat == "SIMPLE":
        counts[size_bin]["RAW_SIMPLE"] += 1
        counts["Total"]["RAW_SIMPLE"] += 1

    if CNV_type == "DEL":
        counts[size_bin]["RAW_DEL"] += 1
        counts["Total"]["RAW_DEL"] += 1

    if repeat != "SIMPLE" and CNV_type == "DEL":
        counts[size_bin]["RAW_DEL_REPEAT"] += 1
        counts["Total"]["RAW_DEL_REPEAT"] += 1

    if repeat == "SIMPLE" and CNV_type == "DEL":
        counts[size_bin]["RAW_DEL_SIMPLE"] += 1
        counts["Total"]["RAW_DEL_SIMPLE"] += 1

    if CNV_type == "DUP":
        counts[size_bin]["RAW_DUP"] += 1
        counts["Total"]["RAW_DUP"] += 1

    if repeat != "SIMPLE" and CNV_type == "DUP":
        counts[size_bin]["RAW_DUP_REPEAT"] += 1
        counts["Total"]["RAW_DUP_REPEAT"] += 1

    if repeat == "SIMPLE" and CNV_type == "DUP":
        counts[size_bin]["RAW_DUP_SIMPLE"] += 1
        counts["Total"]["RAW_DUP_SIMPLE"] += 1

print("\t".join(["BIN", "COUNT", "STATUS"]))

for size_bin in size_bins:
    for category in categories:
        print("\t".join([size_bin, str(counts[size_bin][category]), category]))
