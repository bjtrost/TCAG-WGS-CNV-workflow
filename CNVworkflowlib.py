

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
