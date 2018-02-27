import sys
from VCF import VCF
from VCF import VCF_filter

ZYG={'0/1':'HET', '1/1':'HOMO', './.':'UNKNOWN', '0/0':'REF'}

class MANTA (VCF):
        def __init__ (self, name):
                super (MANTA, self).__init__(name)

	def filter_func(self, chrom, beg, end, type, E):
		if int(end) - int (beg) < 50: return True

		return False

	def run (self, what=None):
		super (MANTA, self).run('MANTA', False, self.filter_func, what=what)
# =========================================================================

class MANTA_filter:
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

class MANTA_filter_mix:
	def __init__ (self, fname, use_default, min_pe=5, min_sr=5, min_qual=0.0):
		self.tab_file=fname

		self.qual=min_qual
		self.min_pe=float(min_pe)
		self.min_sr=float(min_sr)
		self.allelic_ratio=0.2

                if self.qual=='DEFAULT':
                        self.qual=-1
                else:
                        self.qual=float(self.qual)

                print >> sys.stderr, 'min_pe =',   min_pe
                print >> sys.stderr, 'min_sr =',   min_sr
                print >> sys.stderr, 'min_qual =', min_qual

	def is_good (self, E):
		h={}
		F=E[MANTA_filter.other_col].split ('|')
		info=F[0].split (';')

		names=F[1].split(':')
		values=F[2].split(':')

		for a in info:
			b=a.split ('=')

			if len(b) == 1: val = b[0]
			else:
				val=b[1]

			h[b[0]]=val

		for i, nm in enumerate (names):
			h[nm] = values[i]

		if not MANTA_filter.is_good (h):
			return False

		if not h.has_key ('PR'): h['PR']="0,0"
		if not h.has_key ('SR'): h['SR']="0,0"

		h['PR']=map(lambda x: int(x), h['PR'].split(','))
		h['SR']=map(lambda x: int(x), h['SR'].split(','))

		ref_pe=int(h['PR'][0])
		alt_pe=int(h['PR'][1])

		ref_sr=int(h['SR'][0])
		alt_sr=int(h['SR'][1])

		total=1.0*(ref_pe+alt_pe+ref_sr+alt_sr)

		#if (alt_pe + alt_sr) / total < self.allelic_ratio:
		#	return False

                flag=(alt_pe >= self.min_pe or
			(alt_pe+alt_sr) >= (self.min_pe+self.min_sr)/2.0 or 
                        alt_sr >= self.min_sr)

                #print >> sys.stderr, ref_pe, ref_sr, alt_pe, alt_sr, (alt_pe + alt_sr) / total, flag

                if not flag: return False

                if self.qual > -1:
                        if float (h['QUAL']) < self.qual:
			     return False
                else:
                        # Note we use the overall filter
                        if not E[MANTA_filter.filter_col] == 'PASS':
                             return False

		return True
	# -------------------------------------------------------------- |

	def dofilter (self, outfile='/dev/stdout', Filter=False):
		in_f=open (self.tab_file)
		out_f=open (outfile, 'w')

		for line in in_f:
			E=line.rstrip('\n').split()

			if not self.is_good(E): continue

			print line,

		out_f.close()
		in_f.close()
# =======================================================================

class MANTA_filter_qual:
	def __init__ (self, fname, use_default, min_pe=5, min_sr=5, min_qual=0.0):
		self.tab_file=fname
		self.qual=float(min_qual)

	def is_good (self, E):
		h={}
		F=E[MANTA_filter.other_col].split ('|')
		info=F[0].split (';')

		names=F[1].split(':')
		values=F[2].split(':')

		for a in info:
			b=a.split ('=')

			if len(b) == 1: val = b[0]
			else:
				val=b[1]

			h[b[0]]=val

		for i, nm in enumerate (names):
			h[nm] = values[i]


		if not MANTA_filter.is_good (h):
			return False

		if float (h['QUAL']) < self.qual: 
			return False

		return True
	# -------------------------------------------------------------- |

	def dofilter (self, outfile='/dev/stdout', Filter=False):
		in_f=open (self.tab_file)
		out_f=open (outfile, 'w')

		for line in in_f:
			E=line.rstrip('\n').split()

			if not self.is_good(E): continue

			print line,

		out_f.close()
		in_f.close()
# =======================================================================


class MANTA_filter_default:
	def __init__ (self, fname, use_default, min_dv=5, min_pe=5, min_sr=5, min_qual=0):

		self.tab_file=fname
		self.dv=min_dv
		self.sr=min_sr
		self.pe=min_pe
		self.use_default=use_default

	def is_good (self, E):
		if self.use_default:
			if E[MANTA_filter.filter_col] == 'PASS': 
				return True
			else:
				return False

		return True
	# -------------------------------------------------------------- |

	def dofilter (self, outfile='/dev/stdout', Filter=False):
		in_f=open (self.tab_file)
		out_f=open (outfile, 'w')

		for line in in_f:
			E=line.rstrip('\n').split()

			if not self.is_good(E): continue

			print line,

		out_f.close()
		in_f.close()
# =======================================================================

class MANTA_param:
	chrom_col=0
	start_col=1
	end_col=2
	size_col=3
	type_col=4
	filter_col=5
	prog_col=6
	other_col=7
	diff_col=8
	confirm_col=9

	def __init__ (self, fname):
		self.tab_file=fname
	# ------------------------------------------------

	def get_param (self, E):
		h={}
		F=E[MANTA_param.other_col].split ('|')

                gt_name=F[-2]
		gt_info=F[-1]

                names=gt_name.split(':')
		GT=gt_info.split(':')

                for nm, val in zip(names, GT):
                        h[nm]=val

                pr_ref='0'
                pr_var='0'
                sr_ref='0'
                sr_var='0'

                if h.has_key ('PR'):
                        X=h['PR'].split(',')
                        pr_ref=X[0]
                        pr_var=X[1]

                if h.has_key ('SR'):
                        X=h['SR'].split(',')
                        sr_ref=X[0]
                        sr_var=X[1]

		for ele in F[0].split(';'):
			print >> sys.stderr, ele

			x=ele.split('=')
			if len (x) > 1:
				h[x[0]] = x[1]
			else:
				h[x[0]] = x[0]


                zygosity=E[7].split('|')
                zygosity=zygosity[-1].split(':')

		confirm_count = 0
		if not E[MANTA_param.confirm_col] == 'NotApplicable':
			X=E[MANTA_param.confirm_col].split('|')
			confirm_count=len(X)

                print >> sys.stderr, pr_ref, pr_var, sr_ref, sr_var
                return '\t'.join ( [':'.join (E[0:3]), h['QUAL'] ]  + \
                                    [pr_ref, pr_var, sr_ref, sr_var]+ \
                                    [str(confirm_count)] + \
                                    [ZYG[zygosity[0]]])
	# ------------------------------------------------

	def run (self, outfile='/dev/stdout'):
		in_f=open (self.tab_file)
		out_f=open (outfile, 'w')
		print >> out_f, '\t'.join ( [ 'SIZE', 
					      'QUAL',
					      'PR_REF', 
					      'PR_VAR',
					      'SR_REF', 
					      'SR_VAR',
					      'COUNT',
                                              'ZYG'] )
		for line in in_f:
			E=line.rstrip('\n').split()
			
			x=self.get_param (E)

			print >> out_f, x

		out_f.close()
		in_f.close()
# =======================================================================
