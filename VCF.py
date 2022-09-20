import re
from GenomeCoordinate import GenomeCoord

class VCFCol:
    names=['CHROM', 
           'POS', 
           'ID', 
           'REF', 
           'ALT', 
           'QUAL', 
           'FILTER', 
           'INFO', 
           'FORMAT', 
           'SAMPLE']

    CHROM=0
    POS=1
    ID=2
    REF=3
    ALT=4
    QUAL=5
    FILTER=6
    INFO=7
    FORMAT=8
    SAMPLE=9
    #---

    @staticmethod
    def get_info (E, name):
        elements = E[VCFCol.INFO].split(';')
        for ele in elements:
            if re.search (name+'=', ele):
                fields = ele.split('=')
                return fields[1]
        return None
#####################################################################

class VcfRecord:
    def __init__ ( self, E, sample, index = None ):
        self.coord =  GenomeCoord ( E[VCFCol.CHROM], E[VCFCol.POS] )
        self.E = E
        self.sample= sample # from which sample this call is from; setup from the beginning and will be passed on
        self.index = index  # from which file the record coming from

    def get_record ( self ):
        return '\t'.join(self.E)

    def __getitem__ ( self, id ):
        return self.E[id]

    def __str__ ( self ):
        return str (self.coord)

    def __eq__ ( self, other ):
        return self.coord == other.coord

    def __gt__ ( self, other ):
        return self.coord > other.coord

    def __lt__ ( self, other ):
        return self.coord < other.coord

    def __ge__ ( self, other ):
        return self.coord >= other.coord

    def __le__ ( self, other ):
        return self.coord <= other.coord
# ------------------------------------------------------------------- |

