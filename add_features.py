#program annotate a merged erds results with merged cnvnator results
import os
import sys
import re
import getopt
import functions

#file names
input_file_name = ""
output_file_name = ""
sample_id = ""
annot_file_name = ""

i_map = {}
alt_id = {}
#% overlap cutoff
overlap_cutoff = 0
overlap_type = "oneway"

##########
#print usage
def usage():
	"""Print usage"""
	print "usage: program -i <merged erds results> -o <erds+ results> -a <merged cnvnator results> -s <sample_id> -c <overlap cutoff> -p <reciprocal|oneway>"

##########	
#read input file - 
##SampleID	Chr	Start	End	Size	....	
def read_input_map(input_file,i_map,sample_id):
	i_file = open(input_file)
	for line in i_file:
		if line[0] == '#' or line[0] == '>':
			continue		
		line = line.replace("\n","")
		words = line.split("\t")
		sample = sample_id
		chrm = words[1].replace("chr","")
		start=words[2].replace(",","").replace("\"","")
		end=words[3].replace(",","").replace("\"","")
		id = chrm+":"+start+"-"+end
		if i_map.has_key(sample):
			if not i_map[sample].has_key(id):
				if int(start) == int(end):
					i_map[sample][id]=[chrm,int(start)-1,int(end)]
				else:
					i_map[sample][id]=[chrm,int(start),int(end)]
		else:
			if int(start) == int(end):	
				i_map[sample] = {id:[chrm,int(start)-1,int(end),sample]}
			else:
				i_map[sample] = {id:[chrm,int(start),int(end),sample]}
	i_file.close()

##########
#check if one set of coordinates overlaps with another set
def find_overlap(map_coords, file_to_read, results, overlap_cutoff, sample_id, overlap_type):
	ref_coords = {}
	ref_coords_by_chrm = {}
	read_input_ref(file_to_read, ref_coords, ref_coords_by_chrm, sample_id)	
	
	#print "#Number of record found in ", file_to_read, ":", len(ref_coords)
	for sample in map_coords.keys():
		alt_sample_id = sample_id

		for key_map in map_coords[sample]:
			temp_key = map_coords[sample][key_map]
			chrm = temp_key[0].replace("chr","")
			s_1 = temp_key[1]
			e_1 = temp_key[2]
			cnv_len = e_1 - s_1 + 1
			count = 0
			ref_dict = {}			
			cnv_str = ""
			sample_str = ""
			
			#check if the reference set has the sample
			id = key_map
			
			if not ref_coords_by_chrm.has_key(alt_sample_id):
				if results.has_key(sample):
					results[sample][id] = [0,"*","*","*","*"]
				else:
					results[sample]={id:[0,"**","**","**","**"]}
				continue
				
			#check if the reference sample has calls for a chrm	
			if not ref_coords_by_chrm[alt_sample_id].has_key(chrm):
				if results.has_key(sample):
					results[sample][id] = [0,"#","#","#","#"]
				else:
					results[sample]={id:[0,"##","##","##","##"]}
				continue
				
			#check overlaps
			for key_ref in ref_coords_by_chrm[alt_sample_id][chrm]:
				temp_ref = ref_coords[alt_sample_id][key_ref]
				#check if the coordinates are for the same chromosome
				if not temp_ref[0] == temp_key[0]:
					continue
					
				temp = functions.reciprocal_overlap(s_1,e_1,temp_ref[1],temp_ref[2])
				if overlap_type == "reciprocal":
					if temp[0] > overlap_cutoff and temp[1] > overlap_cutoff:
						cnv_str += alt_sample_id + "|" + chrm + ":" + `temp_ref[1]` + "-" + `temp_ref[2]` + "|" + temp_ref[3] + ","
				else:
					if temp[0] > overlap_cutoff:
						cnv_str += alt_sample_id + "|" + chrm + ":" + `temp_ref[1]` + "-" + `temp_ref[2]` + "|" + temp_ref[3] + ","
							
			if cnv_str == "":
				if results.has_key(sample):
					results[sample][id] = [0,"-","-","-","-","-"]
				else:
					results[sample]={id:[0,"-","-","-","-"]}

			else:
				cnv_str = cnv_str[:-1]	
				ref_str = ""
		
				#remove redundant ids
				u_cnv = dict.fromkeys(cnv_str.split(",")).keys()
				o_data = []
				a_o_data = []
				for u in u_cnv:
					u1 = u.split("|")[1].split(":")[1].split("-")
					o_data.append([int(u1[0]),int(u1[1])])
					a_o_data.append([int(u1[0]),int(u1[1])])

				o_data.sort(functions.sort_list)
				a_o_data.sort(functions.sort_list)

			    	c_data = []
				a_c_data = []
			    	functions.cluster(o_data,c_data,s_1,e_1)
			    	functions.alt_cluster(a_o_data,a_c_data)

			    	covered = 0
				a_length = 0
			    	for c in c_data:
					covered += c[1]-c[0]+1
				for c in a_c_data:
					a_length += c[1]-c[0]+1	

				u_count = len(dict.fromkeys(cnv_str.split(",")).keys())
				if results.has_key(sample):
					results[sample][id] = [u_count,cnv_str, covered, covered/float(cnv_len), covered/float(a_length)]
				else:
					results[sample]={id:[u_count,cnv_str, covered, covered/float(cnv_len), covered/float(a_length)]}
				
##########		
#read input file
def read_input_ref(file_to_read, ref_coords, ref_coords_by_chrm, sample_id):
	i_file = open(file_to_read)
	##sample	chrm	start	end	cnv	size	normalized_rd	e-val2	q0
	for line in i_file:	
		if line[0]=="#" or line[0] == '>':
			continue
			
		line = line.replace("\n","")
		words = line.split("\t")
		sample = sample_id
		chrm = words[1].replace("chr","")
		id = chrm + ":" + words[2] + "-" + words[3]
		info = words[4]+"|"+words[6].replace("|","*")
		if ref_coords_by_chrm.has_key(sample):
			if ref_coords_by_chrm[sample].has_key(chrm):
				if not id in ref_coords_by_chrm[sample][chrm]:
					ref_coords_by_chrm[sample][chrm].append(id)
			else:
				ref_coords_by_chrm[sample][chrm]=[id]
				
		else:
			ref_coords_by_chrm[sample]={chrm:[id]}
		
		if ref_coords.has_key(sample):	
			if ref_coords[sample].has_key(id):
				temp = ref_coords[sample][id]
				temp[3] += "," + info
				ref_coords[sample][id] = temp
			else:
				#store - chrm, start, end, annotation, count
				if int(words[2])==int(words[3]):
					ref_coords[sample][id] = [chrm,int(words[2])-1,int(words[3]),info]
				else:
					ref_coords[sample][id] = [chrm,int(words[2]),int(words[3]),info]
		else:
			if int(words[2])==int(words[3]):
				ref_coords[sample]={id:[chrm,int(words[2])-1,int(words[3]),info]}
			else:
				ref_coords[sample]={id:[chrm,int(words[2]),int(words[3]),info]}
	i_file.close()

#####################################################################
#main
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:o:a:s:c:p:h", ["help"])
	except getopt.GetoptError:
		# print help information and exit:
		usage()
		sys.exit(0)
	
	#read command line options
	for o, a in opts:
		if o in ("-h", "--help"):
			usage()
			sys.exit(0)
		#input file
		if o == "-i":
			input_file_name = a			
		#file with calls to annotate with
		if o == "-a":
			annot_file_name = a
		#overlap type to use
		if o == "-p":
			if a != "oneway" and a != "reciprocal":
				print "Only two options allowed - reciprocal and oneway"
				sys.exit(0)
			else:
				overlap_type = a		
		#overlap cutoff to use
		if o == "-c":
			overlap_cutoff = int(a)			
		#sample_id
		if o == "-s":
			sample_id = a			
		#formatted output 	
		if o == "-o":
			output_file_name = a
			
	if input_file_name == "" or output_file_name == ""or annot_file_name == "" or sample_id == "":
		usage()
		sys.exit(0)
	
	#check all input files
        if not os.path.isfile(input_file_name):
                print "Merged erds output not found ..."
                sys.exit(1)
        if not os.path.isfile(annot_file_name):
                print "Merged cnvnator output not found ..."
                sys.exit(1)
	#check output files		
	if os.path.isfile(output_file_name):
                print "Delete file", output_file_name, "and rerun"
		sys.exit(0)
		
	#else open file to write to	
	out = open(output_file_name,'w')
	print output_file_name	
	#datasets
	calls = {}	
		
	#read file to annotate
	read_input_map(input_file_name, i_map, sample_id)
	
	#find overlaps
	find_overlap(i_map, annot_file_name, calls, overlap_cutoff, sample_id, overlap_type)
 		
 	header = "\tcnvn_count\tcnvn_details\tcnvn_coverage\terds_fraction\tcnvn_fraction\tcnv_type_conflict\tcnv_type_confict_coverage" 
 	
	flag = 1	
	i_file = open(input_file_name)
	
	##SampleID	Chr	Start	End	.....
	for line in i_file:
		line = line.replace("\n","")
		line = re.sub("$\t","",line)
		line = re.sub("chr","",line)
		words = line.split("\t")
		if flag == 1:
			flag = 0
			if line[0] == "#" or line[0] == '>':
				print >> out, line + header
				continue
			else:
				h = ""
				for i in words:
					h += "\t"
				print >> out , h + header
			
		sample = sample_id
		start = int(words[2].replace(",","").replace("\"",""))
		end = int(words[3].replace(",","").replace("\"",""))
		id = words[1].replace("chr","")+":"+`start`+"-"+`end`
		type = words[4]
		temp_str = line
		length = end - start 
		#
		if calls[sample][id][0]==0:
			temp_str += "\t"+`calls[sample][id][0]`+"\t"+calls[sample][id][1]+"\t"+calls[sample][id][2]+"\t"+calls[sample][id][3]+"\t"+calls[sample][id][4] + "\t-\t-"
		else:
			ovlp_details = calls[sample][id][1].split(",")
			comment = "-"
			cov_opp_temp = 0
			cov_opp = "-"
			temp = {}
			for d in ovlp_details:
				d_type=d.split("|")[2]	
				if d_type != type:
					boundary = d.split("|")[1]
					start = int(boundary.split(":")[1].split("-")[0])
					end = int(boundary.split(":")[1].split("-")[1])
					cov_opp_temp += end - start

				temp[d_type]=1
			details = "*".join(temp.keys())
			if details != type:
				comment = type + "|" + details
				cov_opp = `cov_opp_temp/float(length)`
				
			temp_str += "\t"+`calls[sample][id][0]`+"\t"+calls[sample][id][1]+"\t"+`calls[sample][id][2]`+"\t"+`calls[sample][id][3]`+"\t"+`calls[sample][id][4]`+"\t"+comment+"\t"+cov_opp
			
		print >> out, temp_str
		
	out.close()	
	i_file.close()
			
