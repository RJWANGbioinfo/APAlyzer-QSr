#!/bin/bash
#===============================================================================
#         FILE: step1_2_3.QC_and_mapping.qcREV.py
#
#        USAGE: python step1_2_3.QC_and_mapping.qcREV.py --<OPTIONS> 
#
#  DESCRIPTION: trim and mapping fastq data from QuantSeq-Rev
#
#      OPTIONS:  --help;                                       show the help message and exit
#                --project(required) <project_name>; define the project name;
#                --rootdir(required) <path_to_your_root_dir>; define the rootdir; 
#                --genodir(required) <path_to_your_genodir_dir>; define the genodir(STAR index) dir for mapping;
#                --refdir(required) <path_to_your_ref_dir>; define the reference dir (ref sequences for trim);#                
#                --threads(required) <# of threads>; define # of threads;
#
# REQUIREMENTS: Pythod_Modules as 'import section'; "STAR" as RNAseq mapper; sambamba;
#
#       AUTHOR: Ruijia Wang  rjwang.bioinfo@gmail.com
#      ADVISOR: Bin Tian     btian@njms.rutgers.edu
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 2019-July-21
#===============================================================================
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
import scipy
import collections
import subprocess
# import multiprocessing
import argparse
import unittest

def getname(targetdir):
	samfiles=sorted(glob.glob(os.path.join(targetdir, '*.fastq')))	
	sample_lst=[]
	for samfile in samfiles:
		url = ntpath.basename(samfile)
		samplename = re.sub('\.fastq$', '', url)	
		sample_lst.append(samplename)
	return sample_lst

##arg setting parts###
parser = argparse.ArgumentParser(description="Use this to run RNAseq-APA automatic pipeline.")
parser.add_argument("--project", action="store", dest='project',default='', metavar='<name_of_project>', help="define the project name")
parser.add_argument("--rootdir", action="store", dest='rootdir',default='', metavar='<rootdir>', help="define the rootdir")
parser.add_argument("--refdir", action="store", dest='refdir',default='', metavar='<refdir>', help="define the refdir")
parser.add_argument("--genodir", action="store", dest='genodir',default='', metavar='<genodir>', help="define the genome dir for mapping")
parser.add_argument("--threads", action="store", dest='threads',default='', metavar='<threads>', help="define the number of threads")
args=parser.parse_args()

#####python setting #####
rootdir=args.rootdir
project=args.project
REFdir=args.refdir
genoDir=args.genodir
CPUS=args.threads
CPUS=str(CPUS)

####automatic setting section####
scrdir=os.path.dirname(os.path.abspath(__file__))
mapper='_star'
tail1='_1'
tail2='_2'
rawfqmdir=os.path.join(rootdir, project+'/rawfastq/')
qcdir=os.path.join(rootdir, project+'/qccheck/')
fqmdir=os.path.join(rootdir, project+'/fastq/')
samoutdir=os.path.join(rootdir, project+'/rawsam/')
rawoutdir=os.path.join(rootdir, project+'/rawout/')
# testdir=os.path.join(rootdir, project+'/test/')
plotdir=os.path.join(rootdir, project+'/plot/')
reportdir=os.path.join(rootdir, project+'/tbl/')
utrdir=os.path.join(rootdir, project+'/tbl/3UTR/')
URdir=os.path.join(rootdir, project+'/tbl/UPS/')

if not os.path.exists(fqmdir):
    os.makedirs(fqmdir)	
if not os.path.exists(qcdir):
    os.makedirs(qcdir)	
if not os.path.exists(samoutdir):
    os.makedirs(samoutdir)	
if not os.path.exists(rawoutdir):
    os.makedirs(rawoutdir)	
# if not os.path.exists(testdir):
    # os.makedirs(testdir)	
if not os.path.exists(plotdir):
    os.makedirs(plotdir)	
if not os.path.exists(reportdir):
    os.makedirs(reportdir)	
if not os.path.exists(utrdir):
    os.makedirs(utrdir)		
if not os.path.exists(URdir):
    os.makedirs(URdir)	
	
samples=getname(rawfqmdir)

for sample in samples:
	print("-------- Step1 QC check of fastq files --------")
	cmd1="fastqc --outdir "+qcdir+' --format fastq --threads '+CPUS+' '+rawfqmdir+"/*.fastq"
	os.system(cmd1)
	
for sample in samples:	
	print("-------- Step2 trim fastq files --------")
	samplerawfq=rawfqmdir+sample+'.fastq'
	samplefq=fqmdir+sample+'.clipped.fastq'	
	cmd2="bbduk.sh in="+samplerawfq+" out="+samplefq+\
		" ref="+REFdir+"polyA.fa.gz"+','+REFdir+"truseq.fa.gz"+\
		" k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20"
	os.system(cmd2)
	
for sample in samples:	
	print("-------- Step3 mapping fastq files:"+sample+" --------")	
	samplefq=fqmdir+sample+'.clipped.fastq'
	samplesam=samoutdir+sample+'.Aligned.out.sam'
	samplebam=samoutdir+sample+'.sorted.bam'
	cmd3='STAR --runThreadN '+CPUS+' --genomeDir '+genoDir+\
	' --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 '+\
	'-outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --outSAMmultNmax 1 --alignIntronMin 20 '+\
	'--alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD '+\
	'--outSAMtype SAM --outReadsUnmapped Fastx --readFilesIn '+ samplefq+ \
	' --outFileNamePrefix '+samoutdir+'/'+sample+'.'
	cmd4="sambamba view -S -t "+CPUS+" -f bam "+samplesam+" | sambamba sort -m 50G -t "+CPUS+" --tmpdir "+samoutdir+" -o "+samplebam+" /dev/stdin"
	cmd5= 'rm '+samplesam
	os.system(cmd3)
	os.system(cmd4)
	os.system(cmd5)

