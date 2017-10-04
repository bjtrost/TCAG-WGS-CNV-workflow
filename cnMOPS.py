
import SV
import os
import re

class cnMOPS:
	seqnames=0
	start=1
	end=2
	width=3
	strand=4
	sampleName=5
	median=6
	mean=7
	CN=8

	@staticmethod
	def parse_CN (CN):
		i=re.sub ('CN', '', CN)
		if int(i) > 2: return 'DUP'
		else:
			return 'DEL'

	@staticmethod
	def doparse ( input_line ):
		E=input_line.rstrip('\n').split('\t')

		chrom=E[cnMOPS.seqnames]
		start=E[cnMOPS.start]
		end=E[cnMOPS.end]
		calltype=cnMOPS.parse_CN( E[cnMOPS.CN] )
		filter=E[cnMOPS.mean]
		program='cn.MOPS'
		other=':'.join(E)
		length=int(E[cnMOPS.width])

		return SV.SV.format_line (chrom, start, end, calltype, filter, program, other, str(length))
	# ------------------------------------------------------------


	def __init__ (self, fname):
		self.fname=fname

		if not os.path.exists (self.fname):
			raise Exception ("NO such file: " + self.fname)
	# -------------------------------------------------------------

	def run ( self ):
		f = open (self.fname)
		for line in f:
			# ignore the header line
			if re.search ('seqnames', line): continue

			sv_line=cnMOPS.doparse ( line )
			print(sv_line)
		f.close()
	# -------------------------------------------------------------
