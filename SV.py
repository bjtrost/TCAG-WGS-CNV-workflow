'''
common format we use for SV project
we do not include translocation in the data.
'''

class SV:
	chrom=0         # chromsome
	start=1         # starting position [start, end] is the call
	end=2           # ending position [start, end] is the call
	calltype=3	# call type, DEL/DUP/INV/INS
	calllength=4	# length of the call, insertion=inserted length
	filter=5	# which parameter used
	caller=6	# which program makes the call
	other=7		# other information, semicolon separated, no space

	@staticmethod
	def format_line (chrom, start, end, calltype, filter, caller='NotApplicable', other='NotApplicable', length=-1):
		if not calltype == 'INS' and length == -1:
			length=int(end) - int (start) + 1
		chrom = chrom.replace("chr", "")
		return '\t'.join ([ chrom, str(start), str(end), str(length), calltype, str(filter), caller, other])

	def __init__ (self, chrom=None, start=None, end=None, length=None, calltype=None, filter=None, caller=None, other=None):
		self.chrom=chrom

		self.start=start
		if start is not None: self.start=int(start)

		self.end=end
		if end is not None: self.end=int(end)

		self.length=length

		self.calltype=calltype
		self.filter=filter
		self.caller=caller
		self.other=other
