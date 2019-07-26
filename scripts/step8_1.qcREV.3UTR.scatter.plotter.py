####import section####
import matplotlib; matplotlib.use('agg')
import itertools
import numpy as np
import pandas as pd
from matplotlib import *
from matplotlib import pyplot
import matplotlib.colors as colors
import os
import sys
import ntpath
import re
import glob
import scipy
import subprocess
from pylab import *
import scipy.stats as stats
import argparse
import unittest
pyplot.switch_backend('agg')
plt.style.use('default')

#####Define the PATH of the folder that hold your .sam, .fa and .png files#####
parser = argparse.ArgumentParser(description="Using this to calculate 3'UTR APA.")
parser.add_argument("--project", action="store", dest='project',default='', metavar='<name_of_project>', help="define the project name")
parser.add_argument("--rootdir", action="store", dest='rootdir',default='', metavar='<rootdir>', help="define the rootdir")
parser.add_argument("--control", action="store", default='', metavar='<list_of_control_group>', help="e.g., 'c('NT0h_C_3', 'NT0h_F_3', 'NT0h_M_3')'  ", nargs='+')
parser.add_argument("--treatment", action="store", default='', metavar='<list_of_treated_group>', help="e.g., 'c('AS0h_C_3', 'AS0h_F_3', 'AS0h_M_3')'  ", nargs='+')

args=parser.parse_args()
rootdir=args.rootdir
project=args.project
controls = ' '.join(args.control).split()
treats = ' '.join(args.treatment).split()

#####################  SETTING ##############
setdir=os.path.join(rootdir, project+'/tbl/3UTR/')
picdir=os.path.join(rootdir, project+'/plot/')


for treat, control in zip(treats, controls):
	title=treat+'_'+control
	tmp=pd.read_csv(setdir+'3mostAPA.'+treat+'_'+control+'.DRPM.fix.tbl',sep='\t')	
	picname1=picdir+treat+'_'+control+'.3mostAPA.DRPM.fix.png'

	typecol=[col for col in tmp.columns if 'pA1.pAutype_' in col][0]
	tmp=tmp.assign(type=tmp[typecol])	
	
	rc('font', weight='bold')
	fig =plt.figure(figsize=(10, 10),dpi=300)
	fig.subplots_adjust(hspace=0.5, wspace=0.3)
	ax=plt.subplot(2,2,1)
	df1=tmp
	typecolname='type'
	dfgray=df1[df1[typecolname]=='NC']
	dfblue=df1[df1[typecolname]=='UP']
	dfred=df1[df1[typecolname]=='DN']
	graycount=str(len(df1[df1[typecolname]=='NC']))
	bluecount=str(len(df1[df1[typecolname]=='UP']))
	redcount=str(len(df1[df1[typecolname]=='DN']))
	plt.scatter(dfgray['Log2Ratio_pA1'], dfgray['Log2Ratio_pA2'], s=5, marker=".", edgecolors='none', c='gray',label='NC')
	plt.scatter(dfred['Log2Ratio_pA1'], dfred['Log2Ratio_pA2'], s=50, marker=".", edgecolors='none', c='r', label='DN')
	plt.scatter(dfblue['Log2Ratio_pA1'], dfblue['Log2Ratio_pA2'], s=50, marker=".", edgecolors='none', c='b',label='UP')
	
	ax.set_xlim([-2.8,2.8])
	ax.set_ylim([-2.8,2.8])	
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.get_xaxis().set_tick_params(direction='out', width=1)
	ax.get_yaxis().set_tick_params(direction='out', width=1)
	titlename="3'UTR APA in "+title
	plt.title(titlename,fontsize=10, bbox={'facecolor':'0.8', 'pad':5},y=1.17)
	plt.xlabel("Log2Ratio Proximal pA")
	plt.ylabel("Log2Ratio Distal pA")
	plt.legend(prop={'size':8},bbox_to_anchor=(1.25, 1.1))
	leggray='NC(gray): '+graycount
	legblue='UP(blue): '+bluecount
	legred='DN(red): '+redcount
	alltext=leggray + ', ' + legblue + ', ' + legred
	text(0.5, 1.03, alltext, ha='center', va='center', transform=ax.transAxes,fontsize=10)
	ratiotext="UP/DN= "+str(round(float(len(df1[df1[typecolname]=='UP']))/len(df1[df1[typecolname]=='DN']),2))
	text(0.5, 1.10, ratiotext, ha='center', va='center', transform=ax.transAxes,fontsize=12)

	pyplot.savefig( picname1 , dpi=300, format='png')
	# print picname1