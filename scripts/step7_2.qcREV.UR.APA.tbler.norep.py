####import section
import matplotlib; matplotlib.use('agg')
import itertools
import numpy as np
import pandas as pd
import os
import sys
import ntpath
import re
import glob
import scipy
from pylab import *
import scipy.stats as stats
from scipy.stats import fisher_exact
import subprocess
import argparse
import unittest
pd.options.mode.chained_assignment = None  # default='warn'

def sum_gene_reads(datainput, colst):
	datasum=datainput.groupby('gene_symbol')[colst].sum()
	datasum=datasum.assign(gene_symbol=datasum.index)
	datasum=datasum.reset_index(drop=True)
	return datasum

def calFisher(dataset):
	x=dataset[0]
	y=dataset[1]
	a=dataset[2]
	b=dataset[3]
	oddsratio, pvalue =fisher_exact([[x, y], [a, b]])
	return pvalue
def UPS_UPDN(datain, treat, control):
	###trimming###
	RPMcols2 = [col for col in datain.columns if 'DRPM_' in col]
	Readcols2 = [col for col in datain.columns if 'num_' in col]
	datain['sum_DRPM'] = datain[RPMcols2].sum(axis=1)
	datain = datain[datain['sum_DRPM']>0]	
	datain3UTR=datain[(datain['pAtype_1']=='F') | (datain['pAtype_1']=='L') | (datain['pAtype_1']=='M')| (datain['pAtype_1']=='S')| (datain['pAtype_1']=='3')]
	datainUPS=datain[(datain['pAtype_1']=='Ic') | (datain['pAtype_1']=='intron') | (datain['pAtype_1']=='CDS')]
	datain3UTRsum=sum_gene_reads(datain3UTR, Readcols2+RPMcols2)
	datainUPSsum=sum_gene_reads(datainUPS, Readcols2+RPMcols2)
	for readscol in Readcols2:
		datain3UTRsum=datain3UTRsum[datain3UTRsum[readscol]>=PASreadscut]
	for readscol in Readcols2:
		datainUPSsum=datainUPSsum[datainUPSsum[readscol]>=PASreadscut]
	
	datain_combine=pd.merge(datainUPSsum, datain3UTRsum, on=['gene_symbol'])
	datain_combine.columns=[col.replace('_x','_pA1') for col in datain_combine.columns.tolist()]
	datain_combine.columns=[col.replace('_y','_pA2') for col in datain_combine.columns.tolist()]
	Rc_PA1COLS = [col for col in datain_combine.columns if 'num_' in col and '_pA1' in col]
	Rc_PA2COLS = [col for col in datain_combine.columns if 'num_' in col and '_pA2' in col]
	for pA1, pA2 in zip(Rc_PA1COLS, Rc_PA2COLS):
		sample=pA1.replace('_pA1','')
		sample=sample.replace('num_','')
		datain_combine['Abn_'+sample+'_pA1']=datain_combine[pA1]/(datain_combine[pA2]+datain_combine[pA1])
		datain_combine['Abn_'+sample+'_pA2']=datain_combine[pA2]/(datain_combine[pA2]+datain_combine[pA1])

	datain_combine['Delta_RA']=datain_combine['Abn_'+treat+'_pA1']-datain_combine['Abn_'+control+'_pA1']	
	datain_combine['p-value']=datain_combine[Rc_PA1COLS+Rc_PA2COLS].apply(calFisher,axis=1)
	datain_combine=datain_combine.assign(status='NC')
	datain_combine['status']=np.where((datain_combine['p-value']<0.05) & (datain_combine['Delta_RA']>0.05), 'UP', datain_combine['status'])
	datain_combine['status']=np.where((datain_combine['p-value']<0.05) & (datain_combine['Delta_RA']<-0.05), 'DN', datain_combine['status'])
	datain_combine=datain_combine.drop_duplicates()
	UPcount=len(datain_combine[datain_combine['status']=='UP'])
	DNcount=len(datain_combine[datain_combine['status']=='DN'])
	datain_combine['Log2Ratio_pA1']=log2(datain_combine['DRPM_'+treat+'_pA1']/datain_combine['DRPM_'+control+'_pA1'])
	datain_combine['Log2Ratio_pA2']=log2(datain_combine['DRPM_'+treat+'_pA2']/datain_combine['DRPM_'+control+'_pA2'])
	# print datain_combine.head()
	datain_combine=datain_combine[['gene_symbol', \
									'num_'+control+'_pA1', 'num_'+control+'_pA2',
									'DRPM_'+control+'_pA1', 'DRPM_'+control+'_pA2',
									'num_'+treat+'_pA1', 'num_'+treat+'_pA2',
									'DRPM_'+treat+'_pA1', 'DRPM_'+treat+'_pA2',																
									'Log2Ratio_pA1','Log2Ratio_pA2','Delta_RA','p-value','status']]	
	datain_combine=datain_combine.assign(RED=datain_combine.Log2Ratio_pA1-datain_combine.Log2Ratio_pA2)
	datain_combine.columns=[col.replace('_pA1','_UPS') for col in datain_combine.columns]
	datain_combine.columns=[col.replace('_pA2','_exon3') for col in datain_combine.columns]
	
	
	return datain_combine, UPcount, DNcount

#####Define the PATH of the folder that hold your .sam, .fa and .png files#####
parser = argparse.ArgumentParser(description="Using this to calculate UPS APA.")
parser.add_argument("--project", action="store", dest='project',default='', metavar='<name_of_project>', help="define the project name")
parser.add_argument("--rootdir", action="store", dest='rootdir',default='', metavar='<rootdir>', help="define the rootdir")
parser.add_argument("--genome", action="store", dest='genome',default='mm9', metavar='<mm9(default)/hg19/rn4/...>', help="define the genome version")
parser.add_argument("--control", action="store", default='', metavar='<list_of_control_group>', help="e.g., 'c('NT0h_C_3', 'NT0h_F_3', 'NT0h_M_3')'  ", nargs='+')
parser.add_argument("--treatment", action="store", default='', metavar='<list_of_treated_group>', help="e.g., 'c('AS0h_C_3', 'AS0h_F_3', 'AS0h_M_3')'  ", nargs='+')

args=parser.parse_args()
rootdir=args.rootdir
project=args.project
geno=args.genome
controls = ' '.join(args.control).split()
treats = ' '.join(args.treatment).split()

#####################  SETTING ##############
tbldir=os.path.join(rootdir, project+'/tbl/UPS/')
pasfile=os.path.join(rootdir, project+'/tbl/',geno+'.pA2gene_usage.DRPM.fix.tbl')	
PASreadscut=2

dfexp=pd.read_csv(pasfile,sep='\t')
norcols=['gene_symbol', 'pAid', 'strand','pAtype_1','gene_Biotype','chromosome','pA_pos']
Readcols2 = [col for col in dfexp.columns if 'num_' in col]
DRPMcol=[col for col in dfexp.columns if 'DRPM_' in col]
dfexp=dfexp[norcols+Readcols2+DRPMcol]
dfexp=dfexp[dfexp.gene_Biotype=='protein-coding']

for treat, control in zip(treats, controls):
	concol=[col for col in dfexp.columns if control in col]
	treacol=[col for col in dfexp.columns if treat in col]
	dfexp_sub=dfexp[norcols+concol+treacol]
	ALL_df, ALL_UPcount, ALL_DNcount=UPS_UPDN(dfexp_sub, treat, control)
	tblfile=tbldir+'UPS_APA.'+treat+'_'+control+'.cut2.tbl'
	####export results####
	ALL_df.to_csv(tblfile, index=None, mode='w', sep='\t')	
	# print treat+'_'+control, ALL_df.RED.median(skipna=True)	
	