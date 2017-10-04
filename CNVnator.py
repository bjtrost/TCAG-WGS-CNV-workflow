
import SV
import os
import re

class CNVnator:
	match_calltype={'deletion':'DEL', 'duplication':'DUP'}
	CNVtype=0
	coordinates=1
	CNVsize=2
	normalizedRD=3
	eval1=4
	eval2=5
	eval3=6
	eval4=7
	q0=8

	@staticmethod
	def parse_coordinates (coordinates):
		coord=re.sub (':', '-', coordinates)

		E=coord.split ('-')

		return E[0], E[1], E[2]

	@staticmethod
	def doparse ( cnvnator_line ):
		E=cnvnator_line.rstrip('\n').split('\t')

		chrom, start, end=CNVnator.parse_coordinates(E[CNVnator.coordinates])
		calltype=CNVnator.match_calltype[E[CNVnator.CNVtype]]
		length=E[CNVnator.CNVsize]
		filter=E[CNVnator.eval1]
		program='CNVnator'
		other=":".join (E)


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
			if line[0] == '#': continue

			sv_line=CNVnator.doparse ( line )
			if sv_line is not None: print(sv_line)
		f.close()
	# -------------------------------------------------------------
