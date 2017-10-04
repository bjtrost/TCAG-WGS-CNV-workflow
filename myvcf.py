import re

# parse vcf file lines
class vcf:
	# varian type
	UNKNOWN=0
	SNV=1
	DEL=2
	INS=3
	MULTI=4

	# vcf file columns
	CHROM=0
	POS=1
	ID=2
	REF=3
	ALT=4
	QUAL=5
	FILTER=6
	INFO=7
	FORMAT=8

	@staticmethod
	def get_chrom (E): return E[vcf.CHROM]

	@staticmethod
	def get_pos   (E): return E[vcf.POS]

	@staticmethod
	def get_id    (E): return E[vcf.ID]

	@staticmethod
	def get_ref   (E): return E[vcf.REF]

	@staticmethod
	def get_alt   (E): return E[vcf.ALT]

	@staticmethod
	def get_qual  (E): return E[vcf.QUAL]

	@staticmethod
	def get_filter(E): return E[vcf.FILTER]

	@staticmethod
	def get_info  (E): return E[vcf.INFO]

	@staticmethod
	def add_info  (E, added):
		E[vcf.INFO]=E[vcf.INFO]+";"+added
	# ---

	@staticmethod
	def get_call_type (E):
		ref=E[vcf.REF]
		alt=E[vcf.ALT]

		if re.search (',', alt): return vcf.MULTI
		if len (ref) > len (alt): return vcf.DEL
		if len (alt) > len (ref): return vcf.INS

		if len (ref) == 1 and len(alt) == 1: return vcf.SNV

		return vcf.UNKNOWN
	# ---


	@staticmethod
	def get_info_field (E, name):
		X=E[vcf.INFO].split(';')

		for x in X:
			z=x.split('=')

			if z[0] == name: return z[1]

		return None
	# ---

	@staticmethod
	def get_info_field_list (E, namelist):
		if not type(namelist) == dict and \
				not type(namelist) == list and \
					not type(namelist) == tuple :
			raise Exception ("vcf.get_info_filed_list 2nd parameter needs a dict/list/tuple type")

		if type(namelist) == list or type(namelist) == tuple:
			h = { }
			for a in namelist:
				h[a] = a

			namelist=h

		X=E[vcf.INFO].split(";")

		rlist=[ ]

		for x in X:
			z=x.split("=")

			if namelist.has_key (z[0]): rlist.append (x)

		return rlist
	# ---

	@staticmethod
	def get_info_dict_field_list (E, namelist):
		if not type(namelist) == dict and \
				not type(namelist) == list and \
					not type(namelist) == tuple :
			raise Exception ("vcf.get_info_filed_list 2nd parameter needs a dict/list/tuple type")

		if type(namelist) == list or type(namelist) == tuple:
			h = { }
			for a in namelist:
				h[a] = a

			namelist=h

		X=E[vcf.INFO].split(";")

		rdict={}

		for x in X:
			z=x.split("=")

			if namelist.has_key (z[0]):
				rdict[z[0]]=z[1]

		return rdict
	# ---

	@staticmethod
	def get_allele_info (E, sample_col, field):
		'''
		Return the call information given the field name.
		FORMAT:
		GT:AB:AD:DP:GQ:MQ0:PL   0/1:0.700:14,6:18:87:0:87,0,297
		'''

		format_fields=E[vcf.FORMAT].split(':')
		data_fields=E[sample_col].split(':')

		for i, nm in enumerate (format_fields):
			if nm == field: return data_fields[i]

		return '.'
	# ---

	@staticmethod
	def is_good (E, QUAL, QD, MQ, FS, DP):
		h=vcf.get_info_dict_field_list (E, ['QD', 'MQ', 'FS', 'DP'])

		return float (E[vcf.QUAL]) >= QUAL    and \
				float (h['QD']) >= QD and \
				float (h['MQ']) >= MQ and \
				float (h['FS']) < FS  and \
				float (h['DP']) >= DP
	# ---

	@staticmethod
	def is_snp (E):
		calls=E[vcf.ALT].split (',')

		# search each call, if there is one call is NOT equal to REF,
		# then return false
		for a in calls:
			if len(E[vcf.REF]) != len(a): return False

		return True
	# ---

vcf_ver=["##fileformat=VCFv4.1"]
vcf_headers=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

human_contig_lines=[
		"##contig=<ID=chr1,length=249250621,assembly=hg19>",
		"##contig=<ID=chr2,length=243199373,assembly=hg19>",
		"##contig=<ID=chr3,length=198022430,assembly=hg19>",
		"##contig=<ID=chr4,length=191154276,assembly=hg19>",
		"##contig=<ID=chr5,length=180915260,assembly=hg19>",
		"##contig=<ID=chr6,length=171115067,assembly=hg19>",
		"##contig=<ID=chr7,length=159138663,assembly=hg19>",
		"##contig=<ID=chr8,length=146364022,assembly=hg19>",
		"##contig=<ID=chr9,length=141213431,assembly=hg19>",
		"##contig=<ID=chr10,length=135534747,assembly=hg19>",
		"##contig=<ID=chr11,length=135006516,assembly=hg19>",
		"##contig=<ID=chr12,length=133851895,assembly=hg19>",
		"##contig=<ID=chr13,length=115169878,assembly=hg19>",
		"##contig=<ID=chr14,length=107349540,assembly=hg19>",
		"##contig=<ID=chr15,length=102531392,assembly=hg19>",
		"##contig=<ID=chr16,length=90354753,assembly=hg19>",
		"##contig=<ID=chr17,length=81195210,assembly=hg19>",
		"##contig=<ID=chr18,length=78077248,assembly=hg19>",
		"##contig=<ID=chr19,length=59128983,assembly=hg19>",
		"##contig=<ID=chr20,length=63025520,assembly=hg19>",
		"##contig=<ID=chr21,length=48129895,assembly=hg19>",
		"##contig=<ID=chr22,length=51304566,assembly=hg19>",
		"##contig=<ID=chrX,length=155270560,assembly=hg19>",
		"##contig=<ID=chrY,length=59373566,assembly=hg19>",
		"##reference=file:///hpf/swscherer_projects1/zwang/data/bwa.hg19/hg19.fa"]
