#!/usr/bin/env python

# Written by Lok Kan Lee
# Modified slightly by Brett Trost

import argparse

####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("depth_filename", type=str, help="samtools depth output file", default=None)
parser.add_argument("-i", "--interval-size", type=int, help="Size of the interval to use for indexing", default=500000)
args = parser.parse_args()
#####################################

##Import from command line
DOCfile = open(args.depth_filename, "r")
index_file = open(args.depth_filename + ".idx", "w")
chrinfo_file = open(args.depth_filename + ".chrinfo", "w")

##Write header for index_file
index_file.write("chr\tlocus\toffset\n")
chrinfo_file.write("chr\tposition of last base\taverage read depth\n")

##Write to index_file helper function
def memwriter(list_of_index):
	for lop in list_of_index:
		index_file.write(lop[0].split(":")[0]+ "\t" + lop[0].split(":")[1] + "\t" + str(lop[1])+ "\n")

chr_sum = 0
list_of_indexes = []
offset_count = 0

##Update variables helper function
def variable_update(line):
	global chr_sum
	global offset_count
	chr_sum += int(line.split("\t")[2])
	offset_count += len(line)

with DOCfile as df:
	line = df.readline()
	chr = line.split("\t")[0]
	pos = int(line.split("\t")[1])
	loc = chr + ":" + str(pos)
	chr_sum = int(line.split("\t")[2])

	list_of_indexes.append([loc, offset_count])
	offset_count += len(line)

	for line in df:
		line_chr = line.split("\t")[0]
		line_pos = int(line.split("\t")[1])
		line_loc = line_chr + ":" + str(line_pos)
		if chr != line_chr:
			avg_read_depth = chr_sum / prev_pos
			memwriter(list_of_indexes)
			chrinfo_file.write("{}\t{}\t{:.3f}\n".format(chr, prev_pos, avg_read_depth))
			chr_sum = int(line.split("\t")[2])
			chr_line_count = 1
			list_of_indexes = [[line_loc, offset_count]]
			offset_count += len(line)
			chr = line_chr
		elif line_pos % args.interval_size == 0:
			list_of_indexes.append([line_loc, offset_count])
			variable_update(line)
		else:
			variable_update(line)
		prev_pos = line_pos
	avg_read_depth = chr_sum / prev_pos
	memwriter(list_of_indexes)
	chrinfo_file.write("{}\t{}\t{:.3f}\n".format(chr, prev_pos, avg_read_depth))
