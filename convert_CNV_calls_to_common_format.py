#!/usr/bin/env python

# Convert calls from Canvas, cn.MOPS, CNVnator, ERDS, Genome STRiP, or RDXplorer to a common format
# Usage example:
# convert_CNV_calls_to_common_format.py input_filename name_of_caller
# name_of_caller must be one of "Canvas", "cn.MOPS", "CNVnator", "ERDS", "Genome_STRiP" (note underscore), or "RDXplorer"

import os
import re
import sys
import argparse

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input_filename", type=str)
parser.add_argument("caller", type=str)
args = parser.parse_args()
#####################################

args.caller = args.caller.replace(".", "").replace(" ", "_") # Convert cn.MOPS to cnMOPS

import_str = "import {}".format(args.caller)
exec(import_str)
run_str = "converter={}.{} (\"{}\")".format(args.caller, args.caller, args.input_filename)
print("Chr\tStart\tEnd\tSize\tType\tAlgorithm-specific filtering data\tAlgorithm\tOther information provided by algorithm")
exec(run_str)
converter.run()
