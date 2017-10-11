
##########
#sort
def sort_list(x,y):
        return cmp(x,y)

##########
#code to calculate reciprocal overlap
def reciprocal_overlap(s_1,e_1,s_2,e_2):
	if s_2 > e_1 or s_1 > e_2:
		return [0,0]	
	else:
		#get the smaller start
		if s_2 >=s_1:
			o_start = s_2
		else:
			o_start = s_1
		#get the smaller end	
		if e_2 >= e_1:
			o_end = e_1
		else:
			o_end = e_2
				
		#calculate length of call and length of overlap
		s1_len = e_1 - s_1
		s2_len = e_2 - s_2
		o_len = o_end - o_start
		
		if 100 * o_len / (s1_len * 1.0) < 0 or 100 * o_len / (s2_len * 1.0) < 0:
			print "s_1: ", s_1, "e_1: ",e_1, "s_2:", s_2, "e_2:", e_2, "o_start:", o_start, "o_end:", o_end
			print "s1_len: ", s1_len, "s2_len: ", s2_len, " o_len: ", o_len, "% s1 length overlap: ", 100 * o_len / (s1_len * 1.0), "% s2 length overlap: ", 100 * o_len / (s2_len * 1.0)
			sys.exit(0)
			
		#return the percent overlap
		return [100 * o_len / (s1_len * 1.0),100 * o_len / (s2_len * 1.0)]
		
##########
#merge overlappping regions into cluster, note that the start and end of the cluster are trimmed 
def cluster(o_data,c_data,ref_start,ref_end):
    START = 0
    END = 0
    clusterString = ""
    #for all regions
    for data in o_data:
        start = data[0]
        end = data[1]
        region = `start`+"-"+`end`+","
        if START == 0 and END == 0:
            START = start
            END = end
            clusterString += region
            continue
        elif start <= END:
            clusterString += region
            #now we have a new cluster end
            if end > END:
                END = end
        #region doesn't overlap with the cluster
        else:
            if START < ref_start:
                START = ref_start
            if END > ref_end:
                END = ref_end
            c_data.append([START,END])
            #start new cluster
            clusterString = region
            START = start
            END = end
    #the last cluster details
    if clusterString != "":
        if START < ref_start:
            START = ref_start
        if END > ref_end:
            END = ref_end
        c_data.append([START,END])
        
##########
#merge overlappping regions into cluster, no start and end cluster trimming 
def alt_cluster(o_data,c_data):
    START = 0
    END = 0
    clusterString = ""
    #for all regions
    for data in o_data:
        start = data[0]
        end = data[1]
        region = `start`+"-"+`end`+","
        if START == 0 and END == 0:
            START = start
            END = end
            clusterString += region
            continue
        elif start <= END:
            clusterString += region
            #now we have a new cluster end
            if end > END:
                END = end
        #region doesn't overlap with the cluster
        else:
            c_data.append([START,END])
            #start new cluster
            clusterString = region
            START = start
            END = end
    #the last cluster details
    if clusterString != "":
        c_data.append([START,END])
        
##########        
#code to calculate overlap
def overlap(s_1,e_1,s_2,e_2):
  if s_2 > e_1 or s_1 > e_2:
    return [0,0]    
  else:
    #get the smaller start
    if s_2 >=s_1:
      o_start = s_2
    else:
      o_start = s_1
    #get the smaller end    
    if e_2 >= e_1:
      o_end = e_1
    else:
      o_end = e_2
                                
    #calculate length of call and length of overlap
    s1_len = e_1 - s_1
    s2_len = e_2 - s_2
    o_len = o_end - o_start
                
    if 100 * o_len / (s1_len * 1.0) < 0 or 100 * o_len / (s2_len * 1.0) < 0:
      print "s_1: ", s_1, "e_1: ",e_1, "s_2:", s_2, "e_2:", e_2, "o_start:", o_start, "o_end:", o_end
      print "s1_len: ", s1_len, "s2_len: ", s2_len, " o_len: ", o_len, "% s1 length overlap: ", 100 * o_len / (s1_len * 1.0), "% s2 length overlap: ", 100 * o_len / (s2_len * 1.0)
      sys.exit(0)
                        
    #return the percent boundary
    return [o_start,o_end]

##########
#find overlap between list of intervals and the region
def find_overlap(intervals,start,end):
  boundaries = []
  c_boundaries = []
  for i in intervals:
    ovlp = overlap(i[0],i[1],start,end)
    if ovlp == [0,0]:
      continue
    else:
      boundaries.append(ovlp)
 
  boundaries.sort(sort_list)  
  cluster(boundaries,c_boundaries,start,end)
  covered = 0
  for c in c_boundaries:
    covered += c[1]-c[0]+1

  return (covered/((end-start+1)*1.0))*100

##########
#find overlap between list of calls and the region
def find_overlap_calls(calls,start,end):
  boundaries = []
  c_boundaries = []
  for i in calls:
    ovlp = overlap(i.get_start(),i.get_end(),start,end)
    if ovlp == [0,0]:
      continue
    else:
      boundaries.append(ovlp)
  
  boundaries.sort(sort_list)
  cluster(boundaries,c_boundaries,start,end)

  covered = 0
  for c in c_boundaries:
    covered += c[1]-c[0]+1

  return (covered/((end-start+1)*1.0))*100
                                          	