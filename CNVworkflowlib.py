import statistics

def mean(list):
    try:
        mean = statistics.mean(list)
    except statistics.StatisticsError:
        mean = 0
    return(mean)

def stdev(list):
    try:
        stdev = statistics.stdev(list)
    except statistics.StatisticsError:
        stdev = 0
    return(stdev)

def depth_search(depth_filename, region):
    doc_file = open(depth_filename)
    index_file = open(depth_filename + ".idx")

    chrom, lower_limit, upper_limit = parse_region(region)

    mem_read_depth = 0
    mem_offset_pos = 0

    #Find offset number for seek()
    index_file.readline()
    for mline in index_file:
        mline = mline.split("\t")
        if chrom == mline[0] and int(mline[1]) <= lower_limit: # Get the offset that is one before the lower limit
            mem_offset_pos = int(mline[2])
        elif chrom == mline[0] and lower_limit < int(mline[1]):
            break

    #In case lower_limit == 1, or lower_limit < first offset bp pos
    if mem_offset_pos == 0:
        index_file.seek(0)
        for mline in index_file:
            mline = mline.split("\t")
            if chrom == mline[0]:
                mem_offset_pos = int(mline[2])
                break

    depth_list = []

    # Make list of read depths in requested region
    doc_file.seek(mem_offset_pos)
    for line in doc_file:
        docline_chr, docline_position, docline_depth = line.strip().split("\t")
        docline_chr = docline_chr.replace("chr", "")
        docline_position = int(docline_position)
        docline_depth = int(docline_depth)
        if chrom == docline_chr and (lower_limit <= docline_position <= upper_limit):
            depth_list.append(docline_depth)
        elif chrom == docline_chr and upper_limit < docline_position:
            break
        elif chrom != docline_chr:
            break

    return(depth_list)

# Read in a two column file, and make a dictionary where the
# key is the first column and the value is the second column
def dict_from_file(filename, delimiter="\t", key_col=0, value_col=1, header=False, value_type="string"):
    f = file_or_stdin(filename)

    if header:
        f.readline()

    d = {}

    for line in f:
        fields = line.split(delimiter)
        if value_type == "int":
            d[fields[key_col].rstrip("\n")] = int(fields[value_col].rstrip("\n"))
        elif value_type == "float":
            d[fields[key_col].rstrip("\n")] = float(fields[value_col].rstrip("\n"))
        else:
            d[fields[key_col].rstrip("\n")] = fields[value_col].rstrip("\n")

    return(d)

def get_depth_lists(chrom, CNV_start, CNV_end, depth_filename, compute_left_left_right_right=False):
    region_template = "{0}:{1}-{2}"

    chrom_size = dict_from_file(depth_filename + ".chrinfo", header=True, value_type="int")

    if CNV_end < CNV_start:
        print("Invalid region entered (start coordinate {} is less than end coordinate {}). Please try again.".format(CNV_start, CNV_end))
        quit()
    if CNV_start < 1:
        print("Invalid region entered (start coordinate is less than 1). Please try again.")
        quit()
    if chrom not in chrom_size:
        print("Invalid chromomosome {} entered. Please try again.".format(chrom))
        quit()
    if CNV_end > chrom_size[chrom]:
        CNV_end = chrom_size[chrom]
        #print("Invalid region entered (end coordinate {} is greater than chromomosome {} size of {:d}). Please try again.".format(CNV_end, chrom, chrom_size[chrom]))
    if CNV_start > chrom_size[chrom]:
        return (None, None, None)

    CNV_size = CNV_end - CNV_start + 1

    left_flank_start = max(1, CNV_start - CNV_size)
    left_flank_end = CNV_start - 1
    left_flank_size = left_flank_end - left_flank_start + 1

    right_flank_start = CNV_end + 1
    right_flank_end = min(right_flank_start + CNV_size - 1, chrom_size[chrom])
    right_flank_size = right_flank_end - right_flank_start + 1

    left_flank_region = region_template.format(chrom, left_flank_start, left_flank_end)
    CNV_region = region_template.format(chrom, CNV_start, CNV_end)
    right_flank_region = region_template.format(chrom, right_flank_start, right_flank_end)

    left_flank_depths = depth_search(depth_filename, left_flank_region)
    CNV_depths = depth_search(depth_filename, CNV_region)
    right_flank_depths = depth_search(depth_filename, right_flank_region)

    depth_lists = {}
    depth_lists["left_flank_depths"] = left_flank_depths + [0] * (left_flank_size - len(left_flank_depths))
    depth_lists["CNV_depths"] = CNV_depths + [0] * (CNV_size - len(CNV_depths))
    depth_lists["right_flank_depths"] = right_flank_depths + [0] * (right_flank_size - len(right_flank_depths))

    if compute_left_left_right_right:  # For calculating read depth for short reads, if CNV is 251-300, compare with 151-200 and 351-400
        left_left_flank_start = max(1, CNV_start - 2*CNV_size)
        left_left_flank_end = max(1, left_flank_start - 1)
        left_left_flank_size = left_left_flank_end - left_left_flank_start + 1

        right_right_flank_start = min(right_flank_end + 1, chrom_size[chrom])
        right_right_flank_end = min(right_right_flank_start + CNV_size - 1, chrom_size[chrom])
        right_right_flank_size = right_right_flank_end - right_right_flank_start + 1

        left_left_flank_region = region_template.format(chrom, left_left_flank_start, left_left_flank_end)
        right_right_flank_region = region_template.format(chrom, right_right_flank_start, right_right_flank_end)

        left_left_flank_depths = depth_search(depth_filename, left_left_flank_region)
        right_right_flank_depths = depth_search(depth_filename, right_right_flank_region)

        depth_lists["left_left_flank_depths"] = left_left_flank_depths + [0] * (left_left_flank_size - len(left_left_flank_depths))
        depth_lists["right_right_flank_depths"] = right_right_flank_depths + [0] * (right_right_flank_size - len(right_right_flank_depths))

    return(depth_lists)

def parse_region(region): # Given a string in the form chr:start-end, return a three-element list containing chr, start, end
    chrom = region.split(":")[0].replace("chr", "").upper()
    start = int(region.split(":")[1].split("-")[0])
    end = int(region.split(":")[1].split("-")[1])

    return(chrom, start, end)

# Get the size bin name for a particular size
def get_size_bin(size):
    if size < 5000:
        return("[1000,5000)")
    elif size < 10000:
        return("[5000,10000)")
    elif size < 100000:
        return("[10000,100000)")
    elif size < 1000000:
        return("[100000,1000000)")
    else:
        return("[1000000,...)")

# Returns a filehandle either from the indicated filename, or if the indicated filename is "-", from standard input
def file_or_stdin(filename):
    if filename == "-":
        return(sys.stdin)
    else:
        return(open(filename))

# Returns ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
def get_chromosomes():
    return([str(c) for c in list(range(1, 23)) + ["X", "Y"]])

def read_common_format(CNV_filename):
    CNV_file = file_or_stdin(CNV_filename)
    header_line = CNV_file.readline().rstrip()
    CNVs = []

    for line in CNV_file:
        this_CNV = {}
        this_CNV["full_line"] = line.rstrip()
        this_CNV["chr"], this_CNV["start"], this_CNV["end"], this_CNV["size"], this_CNV["type"], this_CNV["quality_info"], this_CNV["caller"], this_CNV["caller_specific_info"] = line.rstrip().split("\t")[0:8]

        this_CNV["chr"] = this_CNV["chr"].replace("chr", "").upper()

        if this_CNV["chr"] == "M":
            continue

        this_CNV["start"] = int(float(this_CNV["start"]))
        this_CNV["end"] = int(float(this_CNV["end"]))
        this_CNV["size"] = this_CNV["end"] - this_CNV["start"] + 1

        this_CNV["string_rep"] = this_CNV["type"] + ":" + this_CNV["chr"] + ":" + str(this_CNV["start"]) + "-" + str(this_CNV["end"])
        this_CNV["string_rep_caller"] = this_CNV["caller"] + ":" + this_CNV["string_rep"]

        CNVs.append(this_CNV)

    return(CNVs, header_line)

def get_overlap(a, b): # a = [5,10] b = [7,12]
    return(max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1))

def fifty_percent_reciprocal_overlap(CNV1, CNV2):
    if CNV1["chr"] != CNV2["chr"]:
        return(False)

    if CNV1["type"] != CNV2["type"]:
        return(False)

    overlap = get_overlap([CNV1["start"], CNV1["end"]], [CNV2["start"], CNV2["end"]])
    return(overlap >= (CNV1["size"] / 2) and overlap >= (CNV2["size"] / 2))

# read in a benchmark file, which has the format
#chromosome start_pos   end_pos length  type    technology
def read_benchmark(filename):
    file = file_or_stdin(filename)
    header_line = file.readline().rstrip() # Discard header line

    CNVs = []

    for line in file:
        line = line.rstrip()
        this_CNV = {}
        this_CNV["full_line"] = line.rstrip()
        fields = line.split("\t")

        this_CNV["chr"] = fields[0].replace("chr", "").upper()
        this_CNV["start"] = int(fields[1])
        this_CNV["end"] = int(fields[2])
        this_CNV["size"] = this_CNV["end"] - this_CNV["start"] + 1
        this_CNV["type"] = fields[4]
        this_CNV["caller"] = fields[5]

        this_CNV["string_rep"] = this_CNV["type"] + ":" + this_CNV["chr"] + ":" + str(this_CNV["start"]) + "-" + str(this_CNV["end"])
        this_CNV["string_rep_caller"] = this_CNV["caller"] + ":" + this_CNV["string_rep"]

        CNVs.append(this_CNV)

    return(CNVs, header_line)

def get_overlap_region(CNV1, CNV2):
    return get_overlap([CNV1["start"], CNV1["end"]], [CNV2["start"], CNV2["end"]])
