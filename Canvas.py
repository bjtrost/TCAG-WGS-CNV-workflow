
'''
Convert Canvas output to a common format
'''

import myvcf
import SV
import os
import sys
import re

class Canvas:
	@staticmethod
	def doparse ( canvas_data ):
# Y	59020307	Canvas:GAIN:Y:59020308-59028273	N	<CNV>	2	q10;L10kb	SVTYPE=CNV;END=59028273;CNVLEN=7966	RC:BC:CN:MCC	364:6:7:5
		E=canvas_data

		#if type (canvas_data) is str:
		E = canvas_data.rstrip('\n').split('\t')

		chrom=myvcf.vcf.get_chrom (E)
		start=myvcf.vcf.get_pos (E)

		X=E[2].split (':')

		type=X[1]

		if re.search ('-', X[3]):
			Y=X[3].split('-')
			start=Y[0]
			end=Y[1]
		else:
			start=X[3]
			end=X[4]

		length=int(end) - int (start) + 1
		filter=E[6]
		other=','.join (E)
		program='Canvas'

		if type=='REF':
			X=E[-1].split(":")
			Y=E[-2].split(":")

			type='REFERENCE'
			if Y[-1] == 'MCC' and Y[-2] == 'CN':
				type=type+"_"+X[-2]+"_"+X[-1]

				return None

		if type=="LOSS": type="DEL"
		if type=="GAIN": type="DUP"

		return SV.SV.format_line (chrom, start, end, type,
					filter, program, other, str(length))
	# -------------------------------------------------------------

	def __init__ ( self, fname ):
		self.fname=fname

		if not os.path.exists (self.fname):
			raise Exception ("NO such file: " + self.fname)
	# -------------------------------------------------------------

	def run ( self ):
		f = open (self.fname)
		for line in f:
			if line[0] == '#': continue

			sv_line=Canvas.doparse ( line )
			if sv_line is not None: print(sv_line)

		f.close()
	# -------------------------------------------------------------
