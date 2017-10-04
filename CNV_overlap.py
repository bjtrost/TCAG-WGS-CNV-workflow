#!/usr/bin/env python3

import CNVworkflowlib
import argparse

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("CNV_filenames", type=str) # Comma-separated list of CNV predictions. Can be blank
parser.add_argument("benchmark_filenames", type=str) # Comma-separated list of benchmark CNVs. Can be blank
args = parser.parse_args()
#####################################

overlap_function = CNVworkflowlib.fifty_percent_reciprocal_overlap
CNV_types = ["Deletions", "Duplications"]
filenames = []
CNVs = {}

if args.CNV_filenames:
    for filename in args.CNV_filenames.split(","):
        CNVs[filename], _ = CNVworkflowlib.read_common_format(filename)
        filenames.append(filename)

if args.benchmark_filenames:
    for filename in args.benchmark_filenames.split(","):
        CNVs[filename], _ = CNVworkflowlib.read_benchmark(filename)
        filenames.append(filename)

groups = []
merged_CNV_pointers = {}


for i in range(0, len(filenames) - 1):
    for k in range(0, len(CNVs[filenames[i]])):
        for j in range(i+1, len(filenames)):
            for l in range(0, len(CNVs[filenames[j]])):
                if overlap_function(CNVs[filenames[i]][k], CNVs[filenames[j]][l]):
                    if CNVs[filenames[i]][k]["string_rep_caller"] in merged_CNV_pointers and CNVs[filenames[j]][l]["string_rep_caller"] in merged_CNV_pointers:
                        if merged_CNV_pointers[CNVs[filenames[i]][k]["string_rep_caller"]] != merged_CNV_pointers[CNVs[filenames[j]][l]["string_rep_caller"]]:
                            group1_names = [item["string_rep_caller"] for item in merged_CNV_pointers[CNVs[filenames[i]][k]["string_rep_caller"]]]
                            group2_names = [item["string_rep_caller"] for item in merged_CNV_pointers[CNVs[filenames[j]][l]["string_rep_caller"]]]
                        continue
                    elif CNVs[filenames[i]][k]["string_rep_caller"] in merged_CNV_pointers: # i is already in a group, so just add j
                        group = merged_CNV_pointers[CNVs[filenames[i]][k]["string_rep_caller"]]
                        group.append(CNVs[filenames[j]][l])
                        merged_CNV_pointers[CNVs[filenames[j]][l]["string_rep_caller"]] = group
                    elif CNVs[filenames[j]][l]["string_rep_caller"] in merged_CNV_pointers: # j is already in a group, so just add i
                        group = merged_CNV_pointers[CNVs[filenames[j]][l]["string_rep_caller"]]
                        group.append(CNVs[filenames[i]][k])
                        merged_CNV_pointers[CNVs[filenames[i]][k]["string_rep_caller"]] = group
                    else: # Neither one is in a group yet, so make a new group
                        new_group = [CNVs[filenames[i]][k], CNVs[filenames[j]][l]]
                        groups.append(new_group)
                        merged_CNV_pointers[CNVs[filenames[i]][k]["string_rep_caller"]] = new_group
                        merged_CNV_pointers[CNVs[filenames[j]][l]["string_rep_caller"]] = new_group
        if not CNVs[filenames[i]][k]["string_rep_caller"] in merged_CNV_pointers:
            new_group = [CNVs[filenames[i]][k]]
            groups.append(new_group)
            merged_CNV_pointers[CNVs[filenames[i]][k]["string_rep_caller"]] = new_group


### Now we need to find unmatched CNVs in the last tool ###
last_i = len(filenames) - 1
for k in range(0, len(CNVs[filenames[last_i]])):
    if not CNVs[filenames[last_i]][k]["string_rep_caller"] in merged_CNV_pointers:
        new_group = [CNVs[filenames[last_i]][k]]
        groups.append(new_group)
        merged_CNV_pointers[CNVs[filenames[last_i]][k]["string_rep_caller"]] = new_group

print("Chr\tStart\tEnd\tSize\tType\tNumber of methods detecting this CNV\tList of methods detecting this CNV")

##################################################################################
# Done creating file containing benchmark size distribution
##################################################################################

groups = sorted(groups, key=len, reverse=True)

for group in groups:
    found = False

    # If one of the groups is ERDS, use those breakpoints; else, use the first one
    for caller in group:
        if caller["caller"] == "ERDS":
            chrom = caller["chr"]
            start = caller["start"]
            end = caller["end"]
            size = caller["size"]
            found = True

    if not found:
        chrom = group[0]["chr"]
        start = group[0]["start"]
        end = group[0]["end"]
        size = group[0]["size"]

    group_names = [item["string_rep_caller"] for item in group]
    concat = "\t".join(group_names)
    print("\t".join([str(chrom), str(start), str(end), str(size), group[0]["type"], str(len(group)), ", ".join(sorted(group_names, key=str.lower))]))
