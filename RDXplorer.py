
import os
import SV

class RDXplorer:
	calltypes=['DEL', 'DUP']

	segStart=0
	segEnd=1
	state=2
	Length=3
	copyEst=4
	segMedian=5
	zstat=6
	chrom=7
	posStart=8
	posEnd=9

	@staticmethod
	def doparse (rdx_line):
		E=rdx_line.rstrip('\n').split('\t')

		chrom=E[RDXplorer.chrom]
		start=E[RDXplorer.posStart]
		end=E[RDXplorer.posEnd]
		calltype=RDXplorer.calltypes[int(E[RDXplorer.state]) > 2]
		seglength=str(int(float(E[RDXplorer.posEnd]) - float (E[RDXplorer.posStart])+1))
		filter=E[RDXplorer.zstat]
		caller='RDXplorer'
		other='copyEst=%s;segMedian=%s;zstat=%s' % (E[RDXplorer.copyEst], E[RDXplorer.segMedian], E[RDXplorer.zstat])

		return SV.SV.format_line (chrom, start, end, calltype, filter, caller, other, seglength)
	# -------------------------------------------------------------

	def __init__ (self, fname):
		self.fname=fname

		if not os.path.exists (fname):
			raise Exception ("File: "+fname+" is NOT available")
# -------------------------------------------------------------

	def run ( self ):
		f = open (self.fname)

		header=f.readline()	# ignore the header line

		for line in f:
			sv_line=RDXplorer.doparse ( line )
			print(sv_line)

		f.close()
# -------------------------------------------------------------
