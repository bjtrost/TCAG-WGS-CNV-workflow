
'''
Convert GENOMESTRIP output to a common format
'''
import myvcf
import SV
import os
import re
import sys

class Genome_STRiP:

	@staticmethod
	def doparse ( genomestrip_data ):
		E=genomestrip_data

		#if type (genomestrip_data) is str:
		E = genomestrip_data.rstrip('\n').split('\t')

		#if re.search ('\.\/\.', E[-1]): return None
		#if re.search ('LowQual', E[-1]): return None

		#print >> sys.stderr, E[-1]

		chrom=myvcf.vcf.get_chrom (E)
		start=myvcf.vcf.get_pos (E)
		alt=myvcf.vcf.get_alt(E)

		if re.search (',', E[2]): return None

		# temporary solution
		F=E[-1].split (':')

		if int(F[1]) == 2: return None

		#try:
		#	length=abs(int(myvcf.vcf.get_info_field (E, "SVLEN")))
		#except Exception:
		#	print "NO LENGTH"
		#	print E
		#	sys.exit(1)

		end=myvcf.vcf.get_info_field (E, "END")

		#type=myvcf.vcf.get_info_field (E, "SVTYPE")

		type='DEL'

		if int(F[1]) > 2: type='DUP'

		filter=myvcf.vcf.get_filter (E)
		program="Genome STRiP"
		other='|'.join (E)

		length=int(end)-int(start)+1
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
			sv_line=Genome_STRiP.doparse ( line )
			if sv_line is not None: print(sv_line)
		f.close()
	# -------------------------------------------------------------
