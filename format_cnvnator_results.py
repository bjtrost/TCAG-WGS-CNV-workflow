#format CNVnator output
##CNV_type       coordinates     CNV_size        normalized_RD   e-val1  e-val2  e-val3  e-val4  q0
#deletion        10:43698001-43698700    700     0.302088        0.0409004       5.72382e-13     1       1       0
import re
import sys
import os
if len(sys.argv) != 3:
	print "Usage: <cnvnator output> <formatted results>"
	sys.exit(0)

q0_cutoff = 0.5
q0_exclude = '-1'
nrd_cutoff = 0.03

coi = ["CNV_type","coordinates","CNV_size","normalized_RD","e-val1","e-val2","e-val3","e-val4","q0"]
format = {"deletion":"DEL","duplication":"DUP"}
default_header_index={"CNV_type":0,"coordinates":1,"CNV_size":2,"normalized_RD":3,"e-val1":4,"e-val2":5,"e-val3":6,"e-val4":7,"q0":8}

header_index = {}
for c in coi:
	header_index[c] = -1

in_file = sys.argv[1]
out_file = sys.argv[2]

if not os.path.isfile(in_file):
        print "erds output not found ..."
        sys.exit(1)
if os.path.isfile(out_file):
        print "Delete file", out_file, "and rerun"
	sys.exit(1)      
                  
i_file = open(in_file)
o_file = open(out_file,'w')
if not out_file.find(".txt") == -1:
	f_file = open(out_file.replace(".txt",".filtered.txt"),'w')
else:
	f_file = open(out_file+".filtered.txt",'w')

#headers
print >> o_file, "#sample\tchrm\tstart\tend\tcnv\tsize\tnormalized_rd\te-val2\tq0"
print >> f_file, "#sample\tchrm\tstart\tend\tcnv\tsize\tnormalized_rd\te-val2\tq0"

#extract sample name from in_file name	
sample = re.sub(".*\/","",re.sub(".CNVnator.wg.bin500","",re.sub(".calls.txt","",in_file)))

#format file
header_flag = 0
for line in i_file:
	line = line.replace("\n","")
	words = line.replace("#","").split("\t")
	if line[0] == "#":
		if words[0] == "#CNV_type":
			for i in range (0,len(words)):
				if header_index.has_key(words[i]):
					header_index[words[i]] = i
			for c in coi:
				if header_index[c] == -1:
					print "Column header ", c, " not found in file"
					sys.exit(0)
			
			header_flag = 1
	elif header_flag == 0:
		print "Warning: Header not found, using default.."
		for c in coi:
			header_index[c] = default_header_index[c]
		header_flag = 2
	else:	
		type = format[words[header_index["CNV_type"]]]
		coordinates = words[header_index["coordinates"]].split(":")
		size = words[header_index["CNV_size"]]
		nrd = float(words[header_index["normalized_RD"]])
		eval2 = words[header_index["e-val2"]]
		q0 = words[header_index["q0"]]
		fq0 = float(words[header_index["q0"]])
	
		tmp_str = sample + "\t" + coordinates[0] + "\t" + coordinates[1].split("-")[0] + "\t" + coordinates[1].split("-")[1] + "\t" + type + "\t" + size + "\t" + `nrd` + "\t" + eval2 + "\t" + q0
		print >> o_file, tmp_str
                #filter based on normalized read depth and q0
		if (nrd <= nrd_cutoff and q0 != q0_exclude) or (fq0 <= q0_cutoff and q0 != q0_exclude):
			print >> f_file, tmp_str
i_file.close()
o_file.close()


