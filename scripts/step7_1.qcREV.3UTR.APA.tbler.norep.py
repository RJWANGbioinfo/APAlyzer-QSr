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

def countAPA(dfG):
	dfGcountALL=dfG.groupby('gene_symbol')['pAid'].size()
	dfGcountALL = dfGcountALL.to_frame()
	dfGcountALL = dfGcountALL.assign(PAcount=dfGcountALL['pAid'])
	dfGcountALL = dfGcountALL.assign(gene=dfGcountALL.index)
	dfGcountALL.columns=['del','PAcount', 'gene_symbol']
	del dfGcountALL['del']
	return dfGcountALL
def calFisher(dataset):
	x=dataset[0]
	y=dataset[1]
	a=dataset[2]
	b=dataset[3]
	oddsratio, pvalue =fisher_exact([[x, y], [a, b]])
	return pvalue	
def UTR_UPDN(datain, treat, control):
	###trimming###
	comcol = ['gene_symbol', 'chromosome','pA_pos','strand','pAtype_1','gene_Biotype','gene_desc','pAid']
	treatcol=[col for col in datain.columns if treat in col]
	controlcol=[col for col in datain.columns if control in col]
	datain=datain[comcol+treatcol+controlcol]	
	DRPMcols2 = [col for col in datain.columns if 'DRPM_' in col]
	Readcols2 = [col for col in datain.columns if 'num_' in col]
	datain['sum_DRPM'] = datain[DRPMcols2].sum(axis=1)
	datain=datain[(datain['num_'+treat]>=PASreadscut) & (datain['num_'+control]>=PASreadscut)]
	#APA cutoff#
	datain_count=countAPA(datain)
	datain_count=datain_count[datain_count['PAcount']>APAcount_cut]
	datain_count=datain_count.reset_index(drop=True)
	# print datain.head()
	# print datain_count.head()
	datain=pd.merge(datain, datain_count, on='gene_symbol')
	datain.sort_values(['gene_symbol','sum_DRPM'], ascending=[True,False], inplace=True)
	#calculation#
	datain=datain.groupby('gene_symbol').head(2)
	datainN=datain[datain['strand']=='-']
	datainP=datain[datain['strand']=='+']
	datainN.sort_values(['gene_symbol','pA_pos'], ascending=[True,False], inplace=True)
	datainP.sort_values(['gene_symbol','pA_pos'], ascending=[True,True], inplace=True)
	frames = [datainP, datainN]
	datain = pd.concat(frames)
	datain_pa1 = datain.groupby('gene_symbol').agg(lambda x: x.iloc[0])
	datain_pa1=datain_pa1.assign(gene_symbol=datain_pa1.index)
	datain_pa2 = datain.groupby('gene_symbol').agg(lambda x: x.iloc[1])
	datain_pa2=datain_pa2.assign(gene_symbol=datain_pa2.index)
	datain_pa1=datain_pa1.reset_index(drop=True)
	datain_pa2=datain_pa2.reset_index(drop=True)
	datain_combine=pd.merge(datain_pa1, datain_pa2, on=['gene_symbol','chromosome','strand','gene_Biotype','gene_desc'])
	datain_combine.columns=[col.replace('_x','_pA1') for col in datain_combine.columns.tolist()]
	datain_combine.columns=[col.replace('_y','_pA2') for col in datain_combine.columns.tolist()]
	del datain_combine['sum_DRPM_pA1']
	del datain_combine['sum_DRPM_pA2']	
	DRPM_PA1COLS = [col for col in datain_combine.columns if 'DRPM_' in col and '_pA1' in col]
	DRPM_PA2COLS = [col for col in datain_combine.columns if 'DRPM_' in col and '_pA2' in col]	
	Rc_PA1COLS = [col for col in datain_combine.columns if 'Dreads_' in col and '_pA1' in col]
	Rc_PA2COLS = [col for col in datain_combine.columns if 'Dreads_' in col and '_pA2' in col]
	for pA1, pA2 in zip(DRPM_PA1COLS, DRPM_PA2COLS):
		sample=pA1.replace('_pA1','')
		sample=sample.replace('DRPM_','')
		datain_combine['Abn_'+sample+'_pA1']=datain_combine[pA1]/(datain_combine[pA2]+datain_combine[pA1])
		datain_combine['Abn_'+sample+'_pA2']=datain_combine[pA2]/(datain_combine[pA2]+datain_combine[pA1])
		datain_combine['RE_'+sample]=datain_combine[pA2]/datain_combine[pA1]
		
	datain_combine['Delta_RA']=datain_combine['Abn_'+treat+'_pA1']-datain_combine['Abn_'+control+'_pA1']
	datain_combine['p-value']=datain_combine[Rc_PA1COLS+Rc_PA2COLS].apply(calFisher,axis=1)
	stauts='pA1.pAutype_'+treat+'_'+control
	datain_combine[stauts]='NC'	
	datain_combine[stauts]=np.where((datain_combine['p-value']<0.05) & (datain_combine['Delta_RA']>0.05), 'UP', datain_combine[stauts])
	datain_combine[stauts]=np.where((datain_combine['p-value']<0.05) & (datain_combine['Delta_RA']<-0.05), 'DN', datain_combine[stauts])	
	UPcount=len(datain_combine[datain_combine[stauts]=='UP'])
	DNcount=len(datain_combine[datain_combine[stauts]=='DN'])	
	Readcols_new=[col for col in datain_combine if 'Reads_count_' in col]
	datain_combine['Log2Ratio_pA1']=log2(datain_combine['DRPM_'+treat+'_pA1']/datain_combine['DRPM_'+control+'_pA1'])
	datain_combine['Log2Ratio_pA2']=log2(datain_combine['DRPM_'+treat+'_pA2']/datain_combine['DRPM_'+control+'_pA2'])	
	datain_combine['aUTR_len']=abs(datain_combine['pA_pos_pA2']-datain_combine['pA_pos_pA1'])
	datain_combine=datain_combine[['gene_symbol','chromosome','strand','gene_Biotype','gene_desc','pA_pos_pA1', 'pA_pos_pA2', 'pAid_pA1','pAid_pA2','aUTR_len', \
									'num_'+control+'_pA1', 'num_'+control+'_pA2',
									'DRPM_'+control+'_pA1', 'DRPM_'+control+'_pA2',
									'num_'+treat+'_pA1', 'num_'+treat+'_pA2',
									'DRPM_'+treat+'_pA1', 'DRPM_'+treat+'_pA2',	
									'RE_'+control, 'RE_'+treat,									
									'Log2Ratio_pA1','Log2Ratio_pA2','Delta_RA','p-value',stauts]]
	datain_combine=datain_combine.assign(RED=datain_combine.Log2Ratio_pA2-datain_combine.Log2Ratio_pA1)								
	return datain_combine, UPcount, DNcount

#####Define the PATH of the folder that hold your .sam, .fa and .png files#####
parser = argparse.ArgumentParser(description="Using this to calculate 3'UTR APA.")
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
tbldir=os.path.join(rootdir, project+'/tbl/3UTR/')
pasfile=os.path.join(rootdir, project+'/tbl/',geno+'.pA2gene_usage.DRPM.fix.tbl')

PASreadscut=5
APAcount_cut=1
	
dfpas=pd.read_csv(pasfile,sep='\t')
readscol=[col for col in dfpas.columns if 'num_' in col]
DRPMcol=[col for col in dfpas.columns if 'DRPM_' in col]
Dreadscol=[col.replace('DRPM_','Dreads_') for col in dfpas.columns if 'DRPM_' in col]
samples=[col.replace('DRPM_','') for col in dfpas.columns if 'DRPM_' in col]
dfpas[DRPMcol]=dfpas[DRPMcol].fillna(0)

for sample in samples: 
	dfpas['Dreads_'+sample]=dfpas['DRPM_'+sample].round(0).map(int)
dfpassub=dfpas[['gene_symbol', 'chromosome','pA_pos','strand','pAtype_1','gene_Biotype','gene_desc','pAid']+readscol+Dreadscol+DRPMcol]
dfpas3UTR=dfpassub[(dfpassub['pAtype_1']=='F') |(dfpassub['pAtype_1']=='L')|(dfpassub['pAtype_1']=='M')|(dfpassub['pAtype_1']=='S')]
dfpas3UTR['gene_Biotype']=dfpas3UTR['gene_Biotype'].str.replace('protein_coding','protein-coding')
dfpas3UTR=dfpas3UTR[dfpas3UTR['gene_Biotype']=='protein-coding']

for treat, control in zip(treats, controls):
	dfout, dfout_UPcount, dfout_DNcount=UTR_UPDN(dfpas3UTR, treat, control)
	tblfile=tbldir+'3mostAPA.'+treat+'_'+control+'.DRPM.fix.tbl'
	
	####export results####
	dfout.to_csv(tblfile, index=None, mode='w', sep='\t')	
	# print treat+'_'+control, dfout.RED.median(skipna=True)
