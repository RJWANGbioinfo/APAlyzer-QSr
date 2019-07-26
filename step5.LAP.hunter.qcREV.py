#===============================================================================
#         FILE: LAP.hunter.reverse.py  
#        USAGE: python LAP.hunter.reverse.py
#  DESCRIPTION: Automatically scan all the .Aligned.sortedByCoord.out.bam files in the target folder, 
#				identify all the inital LAP loci, count # of reads support this loci, and plot it
#
# REQUIREMENTS: Pythod_Modules as 'import section'
#       AUTHOR: Ruijia Wang  rw479@njms.rutgers.edu
#      ADVISOR: Bin Tian     btian@njms.rutgers.edu
# ORGANIZATION: 
#      VERSION: 0.1
#      CREATED: 2018-Feb-05
#===============================================================================

####import section
import HTSeq
import itertools
import numpy
import pandas as pd
from pylab import *
import os
import sys
import ntpath
import re
import glob
import collections
import string
import subprocess
import argparse
import unittest

#####Define the PATH of the folder that hold your .sam, .fa and .png files#####
parser = argparse.ArgumentParser(description="Using this to capture po2 position.")
parser.add_argument("--project", action="store", dest='project',default='', metavar='<name_of_project>', help="define the project name")
parser.add_argument("--rootdir", action="store", dest='rootdir',default='', metavar='<rootdir>', help="define the rootdir")
args=parser.parse_args()
rootdir=args.rootdir
project=args.project

#####Define the PATH of the folder that hold your .sam, .fa and .png files#####
samdir=os.path.join(rootdir, project+'/rawsam/')
outdir=samdir+'/LAP_raw/'
if not os.path.exists(outdir):
    os.makedirs(outdir)


samfiles=sorted(glob.glob(os.path.join(samdir, '*.sorted.bam')))
###built counter###
countsALL = collections.Counter( )
for sample in samfiles:
	counts = collections.Counter( )
	url = ntpath.basename(sample)
	samplename = re.sub('\.sorted.bam$', '', url)	
	alignment_file = HTSeq.BAM_Reader( sample )

	for almnt in alignment_file:
		if almnt.aligned:
			if almnt.iv.strand=='+':		
				poii_id=almnt.cigar[0].ref_iv.start_d_as_pos
				counts[ poii_id ] += 1
				countsALL[ poii_id ] += 1
			elif almnt.iv.strand=='-':
				poii_id=almnt.cigar[-1].ref_iv.start_d_as_pos
				counts[ poii_id ] += 1
				countsALL[ poii_id ] += 1

#####Output txt report#####
	with open(str(outdir+samplename+'.LAP.txt'), 'w') as f:
		f.write('PoII_ID\t'+'Chr\t'+'Pos\t'+'Strand\t'+'Support_Reads_Count\n')
		for poii_id in counts:
			# if counts[ poii_id ]<200:
			if counts[ poii_id ]!=0:			
				bb=str(poii_id)
				bb=bb.replace(':','\t')
				bb=bb.replace('/','\t')
				bb=bb.replace('+','ppp')
				bb=bb.replace('-','nnn')
				bb=bb.replace('ppp','-')
				bb=bb.replace('nnn','+')
				f.write(str(poii_id)+ '\t'+bb+'\t'+str(counts[ poii_id ])+'\n')
	f.close()
	print(sample+' finished')

	
with open(str(outdir+'allsample.LAP.txt'), 'w') as f:
	f.write('PoII_ID\t'+'Chr\t'+'Pos\t'+'Strand\t'+'Support_Reads_Count\n')
	for poii_id in countsALL:
		if countsALL[ poii_id ]!=0:			
			bb=str(poii_id)
			bb=bb.replace(':','\t')
			bb=bb.replace('/','\t')
			bb=bb.replace('+','ppp')
			bb=bb.replace('-','nnn')
			bb=bb.replace('ppp','-')
			bb=bb.replace('nnn','+')
			f.write(str(poii_id)+ '\t'+bb+'\t'+str(countsALL[ poii_id ])+'\n')
f.close()
print('All sample finished')
	