#! /usr/bin/env python

import os
import re
import sys
import argparse

BEG=0
END=1
CN=2

VCF_CHR=0
VCF_BEG=1
VCF_END=2
VCF_SIZ=3
VCF_CNV=4
VCF_FILTER=5
VCF_PROG=6
VCF_INFO=7

RAW_NAME=2

# supposedly a[0] <= b[0];
def ov ( a, b ):
        if a[BEG] > b[BEG]:
                raise Exception ('ov wrong tuple'+str(a)+','+str(b))

        if not a[CN] == b[CN]:
                return False

        # only check one direction
        if a[BEG] <= b[BEG] <= a[END]:
                return True

        return False
# ----------------------------------------------------------------------

def get_raw_name (call):
        E=call[-1]
        F=E.split('|')

        return F[RAW_NAME]
# ----------------------------------------------------------------------

def get_cn (call):
        E=call[-1]
        F=E.split('|')
        names=F[-2].split(':')
        values=F[-1].split(':')

        h={}
        for nm, val in zip(names, values):
                h[nm] = val

        if not "CN" in h:
            print >> sys.stderr, F[-2]
            sys.exit(1)
        return int(h['CN'])
# ----------------------------------------------------------------------

def output (current, info):
        num_calls=len(current)

        if num_calls == 1:
                print('\t'.join (current[0]))
                return

        original=get_raw_name (current[0])
        for cl in current[1:]:
                original=original+'+'+get_raw_name(cl)

        # merge calls, and add merging information
        X=current[0][-1].split('|')
        X[VCF_INFO]=X[VCF_INFO]+';'+'ZWANGMERGE='+original

        E=[current[0][VCF_CHR],
           str(info[BEG]),
           str(info[END]),
           str(info[END]-info[BEG]+1),
           current[0][VCF_CNV],
           current[0][VCF_FILTER],
           current[0][VCF_PROG],
           '|'.join (X) ]

        print('\t'.join (E))
# ----------------------------------------------------------------------

def process (calllist):

        if len (calllist) == 0: return

        a=calllist[0]
        current=[ a ]
        prev=[ int(a[VCF_BEG]), int (a[VCF_END]), get_cn(a) ]

        for a in calllist[1:]:
                if ov ( prev, (int(a[VCF_BEG]), int(a[VCF_END]), get_cn(a)) ):
                        current.append (a)

                        prev[BEG]=min(prev[BEG], int(a[VCF_BEG]))
                        prev[END]=max(prev[END], int(a[VCF_END]))
                        continue

                output (current, prev)
                current=[a]
                prev=[ int(a[1]), int (a[2]), get_cn(a) ]

        # output last
        output (current, prev)


#
#
# main script starts from here



####################################
### Parse command-line arguments ###
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("input_filename", type=str)
args = parser.parse_args()
#####################################

f=open(args.input_filename)
print(f.readline().rstrip())
current_chr = None
current_set = [ ]
for line in f:
        E=line.rstrip('\n').split("\t")

        if not current_chr == E[0]:
                process (current_set)

                current_set = [ ]
                current_chr = E[0]

        current_set.append ( E )

process (current_set)
f.close()
