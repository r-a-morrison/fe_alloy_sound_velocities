# Plot PDOS

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

phoxpath = '../070_phox_FeNiSi_man/'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis, and the velocity results will be useful too for 
# calculating the Debye Gruneisen parameter:
EOSpath = '../101_vD_FeNiSi_PowerLaw/Results/input_values.csv'
velpath = '../101_vD_FeNiSi_PowerLaw/Results/vD_results_Constrained_Power.csv'


# Functions
###########

def plotdos(E,dos,ddos,filename):
	fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(10,3))

	ax0.plot(E,dos,color=color)
	ax0.fill_between(E, dos-ddos, dos+ddos, alpha=0.3,facecolor=color)
	ax0.set_xlabel(r'Energy (meV)',fontsize = 16)
	ax0.set_ylabel(r'PDOS $D(E,V)$',fontsize = 16)
	ax0.set_xlim(xmin=0)
	ax0.set_ylim(ymin=-10)
	ax0.tick_params(direction='in',right='on',top='on')

	fig.savefig(filename, format='pdf')
	plt.close()


# Import data
#############

# Load XRD, bulk modulus, and sound velocity results
EOS_df = pd.read_csv(EOSpath, engine='python')
vel_df = pd.read_csv(velpath, engine='python')

input_df = EOS_df.merge(vel_df, on=('Folder','Index','vphi','dvphi'))
# input_df.to_csv('Results/input_values.csv')

# Load PDOS data
# Find the filepath of all .res NRIXS files in phox directories
respath_list = [filepath for filepath in glob.glob(phoxpath+'*/*.res')]

# Prep lists, dictionaries, and df to store input data in
folder_list = []
index_dict = dict()
dos_dict = dict()

# Collect folders, indices, paths, and input values
for respath in respath_list:

	# Determine filepaths for dos and in_psvl
	folder = re.findall('([A-Za-z0-9_]+)/[A-Za-z0-9_]+.res',respath)[0]
	index = re.findall('/([A-Za-z0-9_]+).res',respath)[0]
	dospath = phoxpath+folder+'/Output/'+index+'_dos.dat'

	# Check if each folder is hcp. Don't use it otherwise
	phase = input_df[input_df['Folder']==folder].iloc[-1]['Phase']

	# Import PDOS
	dos_df = pd.read_csv(dospath, sep='\s+', comment='@', header=None,
		names = ['E','DOS','dDOS'])

	# Store to use PDOS later
	folder_list.append(folder)
	index_dict[folder] = index
	dos_dict[folder] = dos_df


# # make PDOS plot
# ################

# def DOSsubplot(folder,number,color):
# 	offset = 150
# 	dos_df = dos_dict[folder]
# 	ax.plot(dos_df['E'], dos_df['DOS']+offset*number,color=color)
# 	ax.fill_between(dos_df['E'], dos_df['DOS']-dos_df['dDOS']+offset*number,
# 		dos_df['DOS']+dos_df['dDOS']+offset*number, facecolor=color, alpha=0.3)
	

# fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(7,10))

# color = '#1f77b4'
# DOSsubplot('2015Mar_DAC13_P6',8,color) # 86
# DOSsubplot('2015Mar_DAC13_P5',7,color) # 69
# DOSsubplot('2015Mar_DAC13_P4',6,color) # 55
# DOSsubplot('2015Mar_DAC13_P3',5,color) # 41
# DOSsubplot('2014Feb_DAC13_P2',4,color) # 37
# DOSsubplot('2014Feb_DAC13_P1',3,color) # 28
# color = 'green'
# DOSsubplot('2013Nov_DAC14_P1',2,color) # 7.1 GPa
# DOSsubplot('2014Feb_DAC15_P1',1,color) # 6.5 GPa
# DOSsubplot('2015Mar_Ambient',0,color)  # 0 GPa

# ax.set_xlabel(r'Energy (meV)',fontsize = 16)
# ax.set_ylabel(r'PDOS $D(E,V)$',fontsize = 16)
# ax.set_xlim([0,70])
# ax.set_ylim(ymin=-10)
# ax.tick_params(direction='in',left='off',top='on')
# ax.set_yticklabels([])



# fig.savefig('PDOS_FeNiSi.pdf', format='pdf', bbox_inches='tight')
# plt.close()


# make hcp PDOS plot
####################

def DOSsubplot2(folder,number,color,ecolor):
	offset = 169
	dos_df = dos_dict[folder]
	ax.errorbar(dos_df['E'], dos_df['DOS']+offset*number,yerr=dos_df['dDOS'],
		marker='.',markersize=2,color=color,ecolor=ecolor,elinewidth=0.5,
		linestyle='none',zorder=-5)

fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(4,10))

color = 'black'
ecolor = 'darkgray'
DOSsubplot2('2015Mar_DAC13_P6',5,color,ecolor) # 86
DOSsubplot2('2015Mar_DAC13_P5',4,color,ecolor) # 69
DOSsubplot2('2015Mar_DAC13_P4',3,color,ecolor) # 55
DOSsubplot2('2015Mar_DAC13_P3',2,color,ecolor) # 41
DOSsubplot2('2014Feb_DAC13_P2',1,color,ecolor) # 37
DOSsubplot2('2014Feb_DAC13_P1',0,color,ecolor) # 28

ax.set_xlabel(r'Energy (meV)',fontsize = 16)
ax.set_ylabel(r'PDOS $D(E,V)$',fontsize = 16)
ax.set_xlim([0,105])
ax.set_ylim(ymin=-10,ymax=1050)
ax.xaxis.set_ticks([0,20,40,60,80,100])
ax.tick_params(direction='in',left='off',top='on')
ax.set_yticklabels([])

fig.savefig('PDOS_hcpFeNiSi_narrow.pdf', format='pdf', bbox_inches='tight')
plt.close()

# make bcc PDOS plot
####################

def DOSsubplot2(folder,number,color,ecolor):
	offset = 240
	dos_df = dos_dict[folder]
	ax.errorbar(dos_df['E'], dos_df['DOS']+offset*number,yerr=dos_df['dDOS'],
		marker='.',markersize=2,color=color,ecolor=ecolor,elinewidth=0.5,
		linestyle='none',zorder=-5)

fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(4,6))

color = 'black'
ecolor = 'darkgray'
DOSsubplot2('2013Nov_DAC14_P1',2,color,ecolor) # 7.1 GPa
DOSsubplot2('2014Feb_DAC15_P1',1,color,ecolor) # 6.5 GPa
DOSsubplot2('2015Mar_Ambient', 0,color,ecolor)  # 0 GPa

ax.set_xlabel(r'Energy (meV)',fontsize = 16)
ax.set_ylabel(r'PDOS $D(E,V)$',fontsize = 16)
ax.set_xlim([0,105])
ax.set_ylim(ymin=-10,ymax=750)
ax.xaxis.set_ticks([0,20,40,60,80,100])
ax.tick_params(direction='in',left='off',top='on')
ax.set_yticklabels([])

fig.savefig('PDOS_bccFeNiSi_narrow.pdf', format='pdf', bbox_inches='tight')
plt.close()
