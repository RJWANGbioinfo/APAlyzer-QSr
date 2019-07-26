####import section
import itertools
import numpy as np
import pandas as pd
from pylab import *
import os
import sys
import ntpath
import re
import glob
import collections
import string
import math
import multiprocessing
import subprocess
import argparse
import unittest
##arg setting parts###
parser = argparse.ArgumentParser(description="summary STAR logs.")
parser.add_argument("--project", action="store", dest='project',default='', metavar='<name_of_project>', help="define the project name")
parser.add_argument("--rootdir", action="store", dest='rootdir',default='', metavar='<rootdir>', help="define the rootdir")

args=parser.parse_args()
rootdir=args.rootdir
project=args.project

samoutdir=os.path.join(rootdir, project+'/rawsam/')
reportdir=os.path.join(rootdir, project+'/tbl/')
outputfile=reportdir+"mapping.summary.txt"

samplelst=[]
inputlst=[]
uniquelst=[]
multilst=[]

logfiles=sorted(glob.glob(os.path.join(samoutdir, '*.Log.final.out')))
for logfile in logfiles:
	url = ntpath.basename(logfile)
	samplename = re.sub('\.Log.final.out$', '', url)
	samplelst.append(samplename)

	f=open( logfile ).readlines()
	for line in f:
		if "Number of input reads" in line:
			# print line
			inputline=line.split( "\t" )[1]
			input=inputline.replace('\n','')
			inputlst.append(input)
		elif "Uniquely mapped reads number" in line:
			uniqueline=line.split( "\t" )[1]
			unique=uniqueline.replace('\n','')
			uniquelst.append(unique)
		elif "Number of reads mapped to multiple loci" in line:
			multiline=line.split( "\t" )[1]
			multi=multiline.replace('\n','')
			multilst.append(multi)
dfsum = pd.DataFrame(
    {'samplename': samplelst,
     'input_reads': inputlst,
     'unique_reads': uniquelst,
	 'multi_reads': multilst
    })		
dfsum=dfsum[['samplename','input_reads','unique_reads','multi_reads']]
dfsum=dfsum.assign(total_mapped_reads=dfsum.unique_reads.map(int)+dfsum.multi_reads.map(int))
dfsum["%mapping"]=(dfsum.total_mapped_reads.map(int)/dfsum.input_reads.map(int)*100).round(2)

dfsum.to_csv(outputfile, index=None, mode='w', sep='\t')