# Plot NRIXS

# Front matter
##############
import os
from os import fdopen, remove
from tempfile import mkstemp
from shutil import move
import glob
import re
import time
import pandas as pd
import numpy as np
from scipy import constants
from scipy.optimize import curve_fit, fsolve
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
from scipy.interpolate import spline
import math
import seaborn as sns

matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16) 

rc = {'lines.linewidth': 1, 
      'axes.labelsize': 20, 
      'axes.titlesize': 20,
      'legend.fontsize': 26,
      'xtick.direction': u'in',
      'ytick.direction': u'in'}
sns.set_style('ticks', rc=rc)

start_time = time.time()


# Input file paths and info
###########################

phoxpath = '../050_phox_Fe_man/'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis, and the velocity results will be useful too for 
# calculating the Debye Gruneisen parameter:
EOSpath = '../081_vD_Fe_PowerLaw/Results/input_values.csv'
velpath = '../081_vD_Fe_PowerLaw/Results/vD_results_Constrained_Power.csv'


# Functions
###########


# Import data
#############

# Load XRD, bulk modulus, and sound velocity results
EOS_df = pd.read_csv(EOSpath, engine='python')
vel_df = pd.read_csv(velpath, engine='python')

input_df = EOS_df.merge(vel_df, on=('Folder','Index','vphi','dvphi'))
# input_df.to_csv('Results/input_values.csv')

# Load NRIXS data
# Find the filepath of all .res NRIXS files in phox directories
respath_list = [filepath for filepath in glob.glob(phoxpath+'*/*.res')]

# Prep lists, dictionaries, and df to store input data in
folder_list = []
index_dict = dict()
NRIXS_dict = dict()

# Collect folders, indices, paths, and input values
for respath in respath_list:

	# Determine filepaths for NRIXS and in_psvl
	folder = re.findall('([A-Za-z0-9_]+)/[A-Za-z0-9_]+.res',respath)[0]
	index = re.findall('/([A-Za-z0-9_]+).res',respath)[0]
	NRIXSpath = phoxpath+folder+'/'+index+'.dat'

	# Check if each folder is hcp. Don't use it otherwise
	phase = input_df[input_df['Folder']==folder].iloc[-1]['Phase']

	# Import NRIXS
	NRIXS_df = pd.read_csv(NRIXSpath, sep='\s+', comment='#', header=None,
		names = ['E','NRIXS','dNRIXS'])

	# Store to use NRIXS later
	folder_list.append(folder)
	index_dict[folder] = index
	NRIXS_dict[folder] = NRIXS_df


# make NRIXS plot
################

# Option 1
# def NRIXSsubplot(folder,number,color):
# 	offset = 3
# 	NRIXS_df = NRIXS_dict[folder]
# 	ax.plot(NRIXS_df['E'], np.log(NRIXS_df['NRIXS'])+offset*number,marker='.',
# 		markersize=2,linewidth=0.5,color=color)
# 	ax.fill_between(NRIXS_df['E'], np.log(NRIXS_df['NRIXS']-NRIXS_df['dNRIXS'])+offset*number,
# 		np.log(NRIXS_df['NRIXS']+NRIXS_df['dNRIXS'])+offset*number, facecolor=color, alpha=0.3)
	
# fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(4,10))

# color = '#1f77b4'
# NRIXSsubplot('2011Feb_171GPa',10,color)
# NRIXSsubplot('2010Aug_151GPa',9,color)	
# NRIXSsubplot('2009Oct_133GPa',8,color)
# NRIXSsubplot('2009Oct_121GPa',7,color)
# NRIXSsubplot('2009Oct_106GPa',6,color)
# NRIXSsubplot('2009Oct_90GPa',5,color)
# NRIXSsubplot('2009Oct_77GPa',4,color)
# NRIXSsubplot('2009Oct_69GPa',3,color)
# NRIXSsubplot('2009Oct_53GPa',2,color)
# NRIXSsubplot('2009Oct_30GPa',1,color)
# color = 'green'
# NRIXSsubplot('2015Mar_Ambient',0,color)

# ax.set_xlabel(r'Energy (meV)',fontsize = 16)
# ax.set_ylabel(r'NRIXS (counts)',fontsize = 16)
# ax.set_xlim([-65,80])
# ax.set_ylim([0,35])
# ax.tick_params(direction='in',left='off',top='on')
# ax.set_yticklabels([])

# fig.savefig('NRIXS_Fe1.pdf', format='pdf', bbox_inches='tight')
# plt.close()


# # Option 5
# def NRIXSsubplot(folder,number,color):
# 	offset = 2
# 	NRIXS_df = NRIXS_dict[folder]
# 	# ax.plot(NRIXS_df['E'], np.log10(NRIXS_df['NRIXS'])+offset*number,marker='.',
# 	# 	markersize=2,linewidth=0.5,color=color,linestyle='none')
# 	# ax.fill_between(NRIXS_df['E'],
# 	# 	np.log10(NRIXS_df['NRIXS']-NRIXS_df['dNRIXS'])+offset*number,
# 	# 	np.log10(NRIXS_df['NRIXS']+NRIXS_df['dNRIXS'])+offset*number, facecolor=color, alpha=0.3)
# 	ax.errorbar(NRIXS_df['E'], NRIXS_df['NRIXS']*10**(offset*number),
# 		yerr=NRIXS_df['dNRIXS']*10**(offset*number),marker='.',
# 		markersize=1.5,linestyle='none',color=color,ecolor='lightgray',elinewidth=1,
# 		zorder=-5)
# 	ax.set_yscale('log')
	
# fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(3.5,10))

# color = 'black'
# NRIXSsubplot('2011Feb_171GPa',9,color)
# NRIXSsubplot('2010Aug_151GPa',8,color)	
# NRIXSsubplot('2009Oct_133GPa',7,color)
# NRIXSsubplot('2009Oct_121GPa',6,color)
# NRIXSsubplot('2009Oct_106GPa',5,color)
# NRIXSsubplot('2009Oct_90GPa',4,color)
# NRIXSsubplot('2009Oct_77GPa',3,color)
# NRIXSsubplot('2009Oct_69GPa',2,color)
# NRIXSsubplot('2009Oct_53GPa',1,color)
# NRIXSsubplot('2009Oct_30GPa',0,color)
# # color = 'green'
# # NRIXSsubplot('2015Mar_Ambient',0,color)

# ax.set_xlabel(r'Energy (meV)',fontsize = 16)
# ax.set_ylabel(r'NRIXS (counts)',fontsize = 16)
# ax.set_xlim([-75,90])
# ax.set_ylim([0.5E1,1E21])
# # ax.tick_params(direction='in',left='off',top='on')
# # ax.set_yticklabels([])

# fig.savefig('NRIXS_Fe5.pdf', format='pdf', bbox_inches='tight')
# plt.close()


# Option 2
# def NRIXSsubplot2(folder,number,color):
# 	offset = 150
# 	NRIXS_df = NRIXS_dict[folder]
# 	ax.plot(NRIXS_df['E'], NRIXS_df['NRIXS']+offset*number,marker='.',
# 		markersize=2,linestyle='none',color=color)
# 	ax.fill_between(NRIXS_df['E'], NRIXS_df['NRIXS']-NRIXS_df['dNRIXS']+offset*number,
# 		NRIXS_df['NRIXS']+NRIXS_df['dNRIXS']+offset*number, facecolor=color, alpha=0.3)
	
# fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(5,10))

# color = '#1f77b4'
# NRIXSsubplot2('2011Feb_171GPa',11,color)
# NRIXSsubplot2('2010Aug_151GPa',10,color)	
# NRIXSsubplot2('2009Oct_133GPa',9,color)
# NRIXSsubplot2('2009Oct_121GPa',8,color)
# NRIXSsubplot2('2009Oct_106GPa',7,color)
# NRIXSsubplot2('2009Oct_90GPa',6,color)
# NRIXSsubplot2('2009Oct_77GPa',5,color)
# NRIXSsubplot2('2009Oct_69GPa',4.5,color)
# NRIXSsubplot2('2009Oct_53GPa',3.25,color)
# NRIXSsubplot2('2009Oct_30GPa',1.5,color)
# color = 'green'
# NRIXSsubplot2('2015Mar_Ambient',0,color)

# ax.set_xlabel(r'Energy (meV)',fontsize = 16)
# ax.set_ylabel(r'NRIXS (counts)',fontsize = 16)
# ax.set_xlim([-65,80])
# ax.set_ylim([0,1800])
# ax.tick_params(direction='in',left='off',top='on')
# ax.set_yticklabels([])

# fig.savefig('NRIXS_Fe2.pdf', format='pdf', bbox_inches='tight')
# plt.close()


# Option 3

def NRIXSsubplot3(ax,folder,number,color,ymin):
	offset = 5
	NRIXS_df = NRIXS_dict[folder]
	# ax.plot(NRIXS_df['E'], NRIXS_df['NRIXS'],marker='.',
	# 	markersize=2,linestyle='none',color=color)
	# ax.fill_between(NRIXS_df['E'], NRIXS_df['NRIXS']-NRIXS_df['dNRIXS'],
	# 	NRIXS_df['NRIXS']+NRIXS_df['dNRIXS'], facecolor=color, alpha=0.3)
	ax.errorbar(NRIXS_df['E'], NRIXS_df['NRIXS'],yerr=NRIXS_df['dNRIXS'],marker='.',
		markersize=1.5,linestyle='none',color=color,ecolor='lightgray',elinewidth=1,
		zorder=-5)
	ax.set_yscale('log')
	ax.tick_params(direction='in',right='on',bottom='off')
	# ax.set_yticklabels([])
	ax.set_xlim([-65,80])
	ax.set_ylim(ymin=ymin,ymax=1.5*10**4)
	ax.tick_params(axis='both', which='major', labelsize=11)
	ax.minorticks_off()
	# ax.xaxis.set_ticks([-60,-40,-20,0,20,40,60,80])


# fig, (ax0,ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10) = plt.subplots(nrows = 11, ncols=1,
# 	figsize=(4,10))

# color = '#1f77b4'
# NRIXSsubplot3(ax0,'2011Feb_171GPa',10,color,4)
# NRIXSsubplot3(ax1,'2010Aug_151GPa', 9,color,6)	
# NRIXSsubplot3(ax2,'2009Oct_133GPa', 8,color,5)
# NRIXSsubplot3(ax3,'2009Oct_121GPa', 7,color,7)
# NRIXSsubplot3(ax4,'2009Oct_106GPa', 6,color,8)
# NRIXSsubplot3(ax5,'2009Oct_90GPa',  5,color,5)
# NRIXSsubplot3(ax6,'2009Oct_77GPa',  4,color,4)
# NRIXSsubplot3(ax7,'2009Oct_69GPa',  3,color,4)
# NRIXSsubplot3(ax8,'2009Oct_53GPa',  2,color,4)
# NRIXSsubplot3(ax9,'2009Oct_30GPa',  1,color,6)
# color = 'green'
# NRIXSsubplot3(ax10,'2015Mar_Ambient',0,color,2)

# ax0.tick_params(direction='in',top='on')
# ax10.tick_params(direction='in',bottom='on')
# ax10.set_xlabel(r'Energy (meV)',fontsize = 16)
# ax5.set_ylabel(r'NRIXS (counts)',fontsize = 16)

# # Fine-tune figure; make subplots close to each other and hide x ticks for
# # all but bottom plot.
# fig.subplots_adjust(hspace=0,wspace=0)
# # plt.setp([ax0.get_xticklabels() for a in fig.axes[:1]], visible=False);

# fig.savefig('NRIXS_Fe3.pdf', format='pdf', bbox_inches='tight')
# plt.close()


# # Option 4
# fig, ((ax5,ax11),(ax4,ax10),(ax3,ax9),(ax2,ax8),(ax1,ax7),(ax0,ax6)) = plt.subplots(nrows = 6,
# 	ncols=2, figsize=(7,10))

# color = '#1f77b4'
# NRIXSsubplot3(ax10,'2011Feb_171GPa',10,color,4)
# NRIXSsubplot3(ax9,'2010Aug_151GPa', 9,color,6)	
# NRIXSsubplot3(ax8,'2009Oct_133GPa', 8,color,5)
# NRIXSsubplot3(ax7,'2009Oct_121GPa', 7,color,7)
# NRIXSsubplot3(ax6,'2009Oct_106GPa', 6,color,8)
# NRIXSsubplot3(ax5,'2009Oct_90GPa',  5,color,5)
# NRIXSsubplot3(ax4,'2009Oct_77GPa',  4,color,4)
# NRIXSsubplot3(ax3,'2009Oct_69GPa',  3,color,4)
# NRIXSsubplot3(ax2,'2009Oct_53GPa',  2,color,4)
# NRIXSsubplot3(ax1,'2009Oct_30GPa',  1,color,6)
# color = 'green'
# NRIXSsubplot3(ax0,'2015Mar_Ambient',0,color,2)

# ax5.tick_params(direction='in',top='on')
# ax10.tick_params(direction='in',top='on')
# ax0.tick_params(direction='in',bottom='on')
# ax6.tick_params(direction='in',bottom='on')
# ax0.set_xlabel(r'Energy (meV)',fontsize = 14)
# ax0.xaxis.set_ticks([-60,-40,-20,0,20,40,60,80])
# ax6.set_xlabel(r'Energy (meV)',fontsize = 14)
# ax6.xaxis.set_ticks([-60,-40,-20,0,20,40,60,80])
# ax2.set_ylabel(r'NRIXS (counts)',fontsize = 14)
# ax8.set_ylabel(r'NRIXS (counts)',fontsize = 14)

# ax11.set_frame_on(False)
# ax11.tick_params(direction='in',left='off',bottom='off')
# ax11.set_xticklabels([])
# ax11.set_yticklabels([])


# # Fine-tune figure; make subplots close to each other and hide x ticks for
# # all but bottom plot.
# fig.subplots_adjust(hspace=0,wspace=0.3)
# # plt.setp([ax0.get_xticklabels() for a in fig.axes[:1]], visible=False);

# fig.savefig('NRIXS_Fe4.pdf', format='pdf', bbox_inches='tight')
# plt.close()

