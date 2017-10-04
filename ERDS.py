
'''
Convert ERDS output to a common format
'''
import myvcf
import SV
import os

class ERDS:

	@staticmethod
	def doparse ( erds_data ):
		E=erds_data
	
		#if type (erds_data) is str:
		E = erds_data.rstrip('\n').split('\t')

		chrom=myvcf.vcf.get_chrom (E)
		start=myvcf.vcf.get_pos (E)

		length=abs(int(myvcf.vcf.get_info_field (E, "SVLEN")))

		end=myvcf.vcf.get_info_field (E, "END")

		type=myvcf.vcf.get_info_field (E, "SVTYPE")
		filter=myvcf.vcf.get_filter (E)
		program="ERDS"
		other=E[-2]+'='+E[-1]

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

			sv_line=ERDS.doparse ( line )
			print(sv_line)
		f.close()
	# -------------------------------------------------------------
