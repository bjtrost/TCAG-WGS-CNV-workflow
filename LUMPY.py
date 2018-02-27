import re, sys
from VCF import VCF

ZYG={'0/1':'HET', '1/1':'HOMO', './.':'UNKNOWN', '0/0':'REF'}

class LUMPY (VCF):
        def __init__ (self, name):
                super (LUMPY, self).__init__(name)

	def filter_func(self, chrom, beg, end, type, E):
                #if re.search ('0/0', E[-1]): return True

		if int(end) - int (beg) < 50: return True

		if not type == 'DEL': return True

		return False

	def run (self, what=None):
		super (LUMPY, self).run('LUMPY', False, self.filter_func, what=what)
# =========================================================================

class LUMPY_filter:
	chrom_col=0 
	start_col=1 
	end_col=2 
	size_col=3 
	type_col=4 
	filter_col=5 
	prog_col=6 
	other_col=7 

	@staticmethod 
	def is_good (h): 
		if h['GT'] == '0/0': return False 

		return True
# =========================================================================

class LUMPY_filter_mix:
	def __init__(self, fname, use_default, min_pe, min_sr, min_qual): 
		self.tab_file=fname 
		self.min_pe=min_pe 
		self.min_sr=min_sr
		self.min_qual=min_qual
		self.allelic_ratio=0.2

		#print >> sys.stderr, 'min_pe =', min_pe
		#print >> sys.stderr, 'min_sr =', min_sr
		#print >> sys.stderr, 'min_qual =', min_qual

		self.support_filter=LUMPY_filter_support(fname, use_default, min_pe, min_sr, min_qual)

	def is_good (self, E):

		h={}
		F=E[-1].split ('|')
		info=F[0].split (';')
		names=F[1].split(':')
		values=F[2].split(':')

		for a in info:
			b=a.split('=')
			
			if len(b) == 1: val=''
			else: val=b[1]

			h[b[0]]=val

		for i, nm in enumerate (names):
			h[nm]=values[i]

		if not self.support_filter.is_good (E): return False

		if float(h['QUAL']) < self.min_qual: return False

		return True
	# ---------------------------------------------------------------- |

	def dofilter (self, outfile='/dev/stdout', Filter=False):
		in_f=open (self.tab_file)
		out_f=open (outfile, 'w')
		for line in in_f:
			E=line.rstrip('\n').split()

			if not self.is_good(E): continue

			print line,

		out_f.close()
		in_f.close()
# ======================================================================== |

class LUMPY_filter_support:
	def __init__(self, fname, use_default, min_pe, min_sr, min_qual): 
		self.tab_file=fname 
		self.min_pe=min_pe 
		self.min_sr=min_sr
		self.allelic_ratio=0.2

	def is_good (self, E):

		h={}
		F=E[-1].split ('|')
		info=F[0].split (';')
		names=F[1].split(':')
		values=F[2].split(':')

		for a in info:
			b=a.split('=')
			
			if len(b) == 1: val=''
			else: val=b[1]

			h[b[0]]=val

		for i, nm in enumerate (names):
			h[nm]=values[i]

		if not LUMPY_filter.is_good (h): return False

		h['SR']=int(h['SR'])
		h['PE']=int(h['PE'])
		h['SU']=int(h['SU'])
		h['AB']=float(h['AB'])

		if h['SR'] < self.min_sr and \
				h['SU'] < (self.min_sr+self.min_pe)/2: 
			return False

		if h['PE'] < self.min_pe and \
				h['SU'] < (self.min_sr+self.min_pe)/2: 
			return False

		if h['AB'] < self.allelic_ratio: return False

		return True

	def dofilter (self, outfile='/dev/stdout', Filter=False):
		in_f=open (self.tab_file)
		out_f=open (outfile, 'w')
		for line in in_f:
			E=line.rstrip('\n').split()

			if not self.is_good(E): continue

			print line,

		out_f.close()
		in_f.close()
# ======================================================================== |

class LUMPY_filter_qual:
	def __init__(self, fname, use_default, min_pe, min_sr, min_qual): 
		self.tab_file=fname 
		self.min_qual=float(min_qual)
		self.use_default=use_default

	def is_good (self, E):

                if self.use_default:
                        return True
		h={}
		F=E[-1].split ('|')
		info=F[0].split (';')
		names=F[1].split(':')
		values=F[2].split(':')

		for a in info:
			b=a.split('=')
			
			if len(b) == 1: val=''
			else: val=b[1]

			h[b[0]]=val

		for i, nm in enumerate (names):
			h[nm]=values[i]

		if not LUMPY_filter.is_good (h): return False

		if float (h['QUAL']) < self.min_qual:
			return False

		return True

	def dofilter (self, outfile='/dev/stdout', Filter=False):
		in_f=open (self.tab_file)
		out_f=open (outfile, 'w')
		for line in in_f:
			E=line.rstrip('\n').split()

			if not self.is_good(E): continue

			print line,

		out_f.close()
		in_f.close()
# ======================================================================== |

class LUMPY_filter_default:
	def __init__(self, fname, use_default, min_sr=5, min_pe=5, min_qual=0.0): 
		self.tab_file=fname 
		self.sr=min_sr 
		self.pe=min_pe 

		# There is NO default filter for LUMPY, NO_FILTER == DEFAULT
		self.use_default=use_default

	def is_good (self, E):

                if self.use_default: return True

		h={}
		F=E[-1].split ('|')
		info=F[0].split (';')
		names=F[1].split(':')
		values=F[2].split(':')

		for a in info:
			b=a.split('=')
			
			if len(b) == 1: val=''
			else: val=b[1]

			h[b[0]]=val

		for i, nm in enumerate (names):
			h[nm]=values[i]

		if not LUMPY_filter.is_good (h):
			return False

		return True

	def dofilter (self, outfile='/dev/stdout', Filter=False):
		in_f=open (self.tab_file)
		out_f=open (outfile, 'w')
		for line in in_f:
			E=line.rstrip('\n').split()

			if not self.is_good(E): continue

			print line,

		out_f.close()
		in_f.close()
# ======================================================================== |

class LUMPY_param:
	def __init__(self, fname):
		self.tab_file=fname 


	def get_param (self, E):

		h={}
		F=E[-3].split ('|')
		info=F[0].split (';')
		names=F[1].split(':')
		values=F[2].split(':')

		for nm, val in zip (names, values):
			h[nm]=val

		for a in info:
			if re.search ('=', a):
				b=a.split('=')
				h[b[0]]=b[1]

		count_confirm = 0
		if not E[-1] == 'NotApplicable':
			X=E[-1].split('|')
			count_confirm=len(X)

                zygosity=E[7].split('|')
                zygosity=zygosity[-1].split(':')

                return '\t'.join ([ ':'.join (E[0:3]),
				    h['QUAL'], 
				    h['PE'], 
				    h['SR'], 
				    h['SU'], 
				    h['AB'],
				    str(count_confirm),
                                    ZYG[zygosity[0]]] )
	# ---------------------------------------------------------------

	def run (self, outfile='/dev/stdout'):
		in_f=open (self.tab_file)
		out_f=open (outfile, 'w')

		print >> out_f, '\t'.join ( [ 'SIZE',
					      'QUAL',
					      'PE',
					      'SR',
					      'SU',
					      'AB',
					      'COUNT',
                                              'ZYG'] )
		for line in in_f:
			E=line.rstrip('\n').split()

			x=self.get_param (E)

			print >> out_f, x

		out_f.close()
		in_f.close()
# ======================================================================== |
