#format ERDS (vcf) output        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	<sample> <sample> <sample> <sample>

import os, sys, re
tags = {}
tag_pos = {}

tags_of_interest = ["FORMAT","INFO"]
for key in tags_of_interest:
    tags[key] = {}
    
def print_tags():
    for t in tags.keys():
        print t, tags[t]
        
if len(sys.argv) != 3:
    print "<program> <erds vcf output> <formatted results>"
    sys.exit(0)
    
vcf_file_name = sys.argv[1]
out_file_name = sys.argv[2]
t_file_name = out_file_name + ".temp"

if not os.path.isfile(vcf_file_name):
        print "erds output not found ..."
        sys.exit(1)
if os.path.isfile(out_file_name):
        print "Delete file", out_file_name, "and rerun"
	sys.exit(1)      

header = []
t_file = open(t_file_name ,'w')
o_file = open(out_file_name,'w')
info_header = {}

vcf_file = open(vcf_file_name)
num_samples = 0
vcf_file = open(vcf_file_name)
for line in vcf_file:
    line = line.replace("\n","")
    if line[0:2] == "##":
        words = line.split("<")
        if len(words) < 2:
            continue
        
        tag = words[0].replace("#","").replace("=","")
        id = words[1].split(",")[0].replace("ID=","")
        if tag in tags_of_interest:
            tags[tag][id]=0

    elif line[0]=="#":
	line = line.replace("POS","START")
        words = line.split("\t")
        header = words
        
        new_header = words[0:7]
        new_header.extend(tags[words[7]].keys())
        num_samples = len(words)-9 + 1
        if num_samples == 0:
            print "Error: no samples in vcf file!!"
            sys.exit(0)
            
        for i in range (1,num_samples):
            id = 8+i
            for key in tags[words[8]].keys():
                new_header.append(words[id] + "|" + key)
        print >> t_file, "\t".join(new_header)
        
    else:
        words = line.split("\t")
        entry = []
        entry = words[0:7]
	#format "INFO" column
        info = words[7].split(";")
        info_dets = {}
        for t in tags[header[7]].keys():
		info_dets[t] = "-"
        for i in info:
		i_1 = i.split("=")
		if i_1[0]=="DB" or i_1[0] == "DS" or i_1[0] == "INV5" or i_1[0] == "INV3":
			continue
		elif i_1[0]=="PRECISE":
                        info_dets["IMPRECISE"]="PRECISE"
			continue
		elif i_1[0]=="IMPRECISE" :
			info_dets["IMPRECISE"]=i_1[0]
			continue
	
		info_dets[i_1[0]]=i_1[1]

	#put together the entry
        for t in tags[header[7]].keys():
		entry.append(info_dets[t])

        #format "FORMAT" column
	format_info = words[8].split(":")
        
        for i in range (1,num_samples):
		id = 8+i
		sample = words[id].split(":")
		s_format = {}
		for t in tags[header[8]].keys():
			s_format[t]="-"
		
		if len(sample) != 1:
			for f in range(0,len(format_info)):
				s_format[format_info[f]]=sample[f]
	
		#put together the entry
        	for t in tags[header[8]].keys():
                	entry.append(s_format[t])   
        print >> t_file, "\t".join(entry)
	
vcf_file.close()
t_file.close()

#format temp file
t_file = open(t_file_name)
##sample	CHROM	START	END	SVTYPE	SIZE	sample|CN	FILTER	sample|REFCN	SVLEN	ALT	IMPRECISE
fixed_column_index = {"CHROM":-1,"START":-1,"END":-1,"SVTYPE":-1,"SAMPLE|CN":-1,"FILTER":-1,"SAMPLE|REFCN":-1,"SVLEN":-1,"ALT":-1,"IMPRECISE":-1,"SIZE":-1}
flex_column_index = {}
sample = re.sub(".*/","",out_file_name).replace(".txt","").replace(".erds.vcf","")

for line in t_file:
	line = line.replace("\n","")	
	if line[0] == "#":
		line = line.replace("#","")
		words = line.split("\t")
		for i in range (0,len(words)):
			words[i]=re.sub(r'.*\|','SAMPLE|',words[i])
			if words[i] in fixed_column_index.keys():
				fixed_column_index[words[i]] = i

		for key in fixed_column_index.keys():
			if fixed_column_index[key] == -1 and key != "SIZE":
				print "Required column " , key, " is missing.."
				sys.exit(0)

		new_header = "#sample\tCHROM\tSTART\tEND\tSVTYPE\tSIZE\tFILTER\tSAMPLE|REFCN\tSAMPLE|CN\tIMPRECISE\tSVLEN\tALT"
		print >> o_file, new_header

	else:
		words = line.split("\t")
		
		entry = []
		chrom = words[fixed_column_index["CHROM"]]
		start = words[fixed_column_index["START"]]
		end = words[fixed_column_index["END"]]
		svtype = words[fixed_column_index["SVTYPE"]]
		filter = words[fixed_column_index["FILTER"]]
		refcn = words[fixed_column_index["SAMPLE|REFCN"]]
		samplecn = words[fixed_column_index["SAMPLE|CN"]]
		imprecise = words[fixed_column_index["IMPRECISE"]]
		svlen = words[fixed_column_index["SVLEN"]]
		alt = words[fixed_column_index["ALT"]]
		
		if end == "-":
			end = start
                        
		entry.extend([sample,words[fixed_column_index["CHROM"]],start,end,svtype,`int(end)-int(start)+1`,filter,refcn,samplecn,imprecise,svlen,alt])
		print >> o_file, "\t".join(entry)

t_file.close()
o_file.close()


