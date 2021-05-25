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

phoxpath_dict = dict()
V0_dict = dict()
Vi_dict = dict()
gamma0_dict = dict()
q_dict = dict()

study = 'bccFe'
phoxpath_dict[study] = '../050_phox_Fe_man/2015Mar_Ambient/Output/Fe_Ambient_dos.dat'

study = 'hcpFe'
phoxpath_dict[study] = '../050_phox_Fe_man/2009Oct_90GPa/Output/Fe_Murphy_P6_dos.dat'
V0_dict[study] = 22.428  # in cubic Angstroms
Vi_dict[study] = 17.10   # in cubic Angstroms
gamma0_dict[study] = 2.04
q_dict[study] = 1

study = 'bccFeNi'
phoxpath_dict[study] = '../060_phox_FeNi_man/2015Mar_Ambient/Output/FeNi_Ambient_dos.dat'

study = 'hcpFeNi'
phoxpath_dict[study] = '../060_phox_FeNi_man/2012Apr_DAC11_P6/Output/FeNi_DAC11_P6_dos.dat'
V0_dict[study] = 22.505  # in cubic Angstroms
Vi_dict[study] = 18.72   # in cubic Angstroms
gamma0_dict[study] = 2.07
q_dict[study] = 1

study = 'bccFeNiSi'
phoxpath_dict[study] = '../070_phox_FeNiSi_man/2015Mar_Ambient/Output/FeNiSi_Ambient_dos.dat'

study = 'hcpFeNiSi'
phoxpath_dict[study] = '../070_phox_FeNiSi_man/2015Mar_DAC13_P3/Output/FeNiSi_DAC13_P3_dos.dat'
V0_dict[study] = 22.952  # in cubic Angstroms
Vi_dict[study] = 19.09   # in cubic Angstroms
gamma0_dict[study] = 2.03
q_dict[study] = 1

# Functions
###########

def calcScalingParam(Vi,V,V0,gamma0,q):
	# V0, gamma0, q are fit parameters
	# Vi is volume of starting DOS
	# V is volume of target DOS
	# If scaling to 0 GPa, V=V0
	scalingparam = np.exp( gamma0 * (Vi/V0)**q * (1/q) * ((V/Vi)**q-1) )
	return scalingparam

def scalePDOS(scalingparam,dos_df):
	newdos_df = dos_df.copy()
	newdos_df['E'] = dos_df['E']/scalingparam
	newdos_df['DOS'] = scalingparam*dos_df['DOS']
	newdos_df['dDOS'] = scalingparam*dos_df['dDOS']
	return newdos_df


# Import data
#############

dos_dict = dict()

# Load PDOS data
for study in phoxpath_dict.keys():

	# Import PDOS
	dos_dict[study] = pd.read_csv(phoxpath_dict[study], sep='\s+', comment='@', header=None,
		names = ['E','DOS','dDOS'])


# make hcp PDOS plot
####################

def DOSsubplot2(dos_df,number,color,ecolor):
	offset = 175
	# ax.plot(dos_df['E'], dos_df['DOS']+offset*number,color=color)
	# ax.fill_between(dos_df['E'], dos_df['DOS']-dos_df['dDOS']+offset*number,
	# 	dos_df['DOS']+dos_df['dDOS']+offset*number, facecolor=color, alpha=0.3)
	ax.errorbar(dos_df['E'], dos_df['DOS']+offset*number,yerr=dos_df['dDOS'],
		marker='.',markersize=1.5,color=color,ecolor=ecolor,elinewidth=0.5,
		linestyle='none',zorder=-5)

def DOSsubplot3(dos_df,number,color,ecolor):
	offset = 175
	ax.plot(dos_df['E'], dos_df['DOS']+offset*number,color=color,linewidth=1)
	# ax.fill_between(dos_df['E'], dos_df['DOS']-dos_df['dDOS']+offset*number,
	# 	dos_df['DOS']+dos_df['dDOS']+offset*number, facecolor=color, alpha=0.3)
	# ax.errorbar(dos_df['E'], dos_df['DOS']+offset*number,yerr=dos_df['dDOS'],
	# 	marker='.',markersize=1.5,color=color,ecolor=ecolor,elinewidth=0.5,
	# 	linestyle='none',zorder=-5)

fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(5,7.5))

study = 'bccFe'
color = 'gray'  # #808080
ecolor = '#D3D3D3'
DOSsubplot2(dos_dict[study],0,color,ecolor)

study = 'hcpFe'
color = 'gray'  # #808080
ecolor = '#D3D3D3'
scalingparam = calcScalingParam(Vi_dict[study],V0_dict[study],V0_dict[study],
	gamma0_dict[study],q_dict[study])
DOSsubplot3(scalePDOS(scalingparam,dos_dict[study]),1,color,ecolor)
# DOSsubplot2(dos_dict[study],1,color,ecolor)

study = 'bccFeNi'
color = 'darkorange' # #FF8C00
ecolor = '#FFD4A6'
DOSsubplot2(dos_dict[study],4,color,ecolor)

study = 'hcpFeNi'
color = 'darkorange' # #FF8C00
ecolor = '#FFD4A6'
scalingparam = calcScalingParam(Vi_dict[study],V0_dict[study],V0_dict[study],
	gamma0_dict[study],q_dict[study])
DOSsubplot3(scalePDOS(scalingparam,dos_dict[study]),5,color,ecolor)
# DOSsubplot2(dos_dict[study],5,color,ecolor)

study = 'bccFeNiSi'
color = 'deepskyblue'  # #00BFFF
ecolor = '#A4ECFF'
DOSsubplot2(dos_dict[study],2,color,ecolor)

study = 'hcpFeNiSi'
color = 'deepskyblue'  # #00BFFF
ecolor = '#A4ECFF'
scalingparam = calcScalingParam(Vi_dict[study],V0_dict[study],V0_dict[study],
	gamma0_dict[study],q_dict[study])
DOSsubplot3(scalePDOS(scalingparam,dos_dict[study]),3,color,ecolor)
# DOSsubplot2(dos_dict[study],3,color,ecolor)

ax.set_xlabel(r'Energy (meV)',fontsize = 16)
ax.set_ylabel(r'PDOS $D(E,V)$',fontsize = 16)
ax.set_xlim([0,60])
ax.set_ylim(ymin=-10,ymax=1100)
ax.xaxis.set_ticks([0,20,40,60])
ax.tick_params(direction='in',left='off',top='on')
ax.set_yticklabels([])

fig.savefig('PDOS_0GPa_compare.pdf', format='pdf', bbox_inches='tight')
plt.close()

# # make bcc PDOS plot
# ####################

# def DOSsubplot2(folder,number,color,ecolor):
# 	offset = 150
# 	dos_df = dos_dict[folder]
# 	# ax.plot(dos_df['E'], dos_df['DOS']+offset*number,color=color)
# 	# ax.fill_between(dos_df['E'], dos_df['DOS']-dos_df['dDOS']+offset*number,
# 	# 	dos_df['DOS']+dos_df['dDOS']+offset*number, facecolor=color, alpha=0.3)
# 	ax.errorbar(dos_df['E'], dos_df['DOS']+offset*number,yerr=dos_df['dDOS'],
# 		marker='.',markersize=1.5,color=color,ecolor=ecolor,elinewidth=1,
# 		linestyle='none',zorder=-5)
# fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(7,10))

# color = 'black'
# ecolor = 'lightgray'
# DOSsubplot2('2011Feb_171GPa',9,color,ecolor)
# DOSsubplot2('2010Aug_151GPa',8,color,ecolor)	
# DOSsubplot2('2009Oct_133GPa',7,color,ecolor)
# DOSsubplot2('2009Oct_121GPa',6,color,ecolor)
# DOSsubplot2('2009Oct_106GPa',5,color,ecolor)
# DOSsubplot2('2009Oct_90GPa',4,color,ecolor)
# DOSsubplot2('2009Oct_77GPa',3,color,ecolor)
# DOSsubplot2('2009Oct_69GPa',2,color,ecolor)
# DOSsubplot2('2009Oct_53GPa',1,color,ecolor)
# DOSsubplot2('2009Oct_30GPa',0,color,ecolor)

# ax.set_xlabel(r'Energy (meV)',fontsize = 16)
# ax.set_ylabel(r'PDOS $D(E,V)$',fontsize = 16)
# ax.set_xlim([0,85])
# ax.set_ylim(ymin=-10)
# ax.xaxis.set_ticks([0,20,40,60,80])
# ax.tick_params(direction='in',left='off',top='on')
# ax.set_yticklabels([])

# fig.savefig('PDOS_bccFe.pdf', format='pdf', bbox_inches='tight')
# plt.close()


# Make plot for presentation
############################

def DOSsubplot2(dos_df,number,color,ecolor):
	offset = 175
	# ax.plot(dos_df['E'], dos_df['DOS']+offset*number,color=color)
	# ax.fill_between(dos_df['E'], dos_df['DOS']-dos_df['dDOS']+offset*number,
	# 	dos_df['DOS']+dos_df['dDOS']+offset*number, facecolor=color, alpha=0.3)
	ax.errorbar(dos_df['E'], dos_df['DOS']+offset*number,yerr=dos_df['dDOS'],
		marker='.',markersize=1.5,color=color,ecolor=ecolor,elinewidth=0.5,
		linestyle='none',zorder=-5)

fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(10,3.5))

study = 'hcpFeNi'
color = 'black' # #FF8C00
ecolor = 'gray'
DOSsubplot2(dos_dict[study],0,color,ecolor)

ax.set_xlabel(r'Energy (meV)',fontsize = 16)
ax.set_ylabel(r'PDOS $D(E,V)$',fontsize = 16)
ax.set_xlim([0,80])
# ax.set_ylim(ymin=-10,ymax=1100)
ax.xaxis.set_ticks([0,20,40,60,80])
ax.tick_params(direction='in',left='off',top='on')
ax.set_yticklabels([])

fig.savefig('PDOS_FeNi_41GPa_Pres.pdf', format='pdf', bbox_inches='tight')
plt.close()


# Make plot for presentation
############################

def DOSsubplot2(dos_df,number,color,ecolor1):
	offset = 175
	# ax.plot(dos_df['E'], dos_df['DOS']+offset*number,color=color)
	# ax.fill_between(dos_df['E'], dos_df['DOS']-dos_df['dDOS']+offset*number,
	# 	dos_df['DOS']+dos_df['dDOS']+offset*number, facecolor=color, alpha=0.3)
	ax.errorbar(dos_df['E'], dos_df['DOS']+offset*number,yerr=dos_df['dDOS'],
		marker='.',markersize=2.5,color=color,ecolor=ecolor,elinewidth=1,
		linestyle='none',zorder=-5)

fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(9,7))

study = 'hcpFe'
color = 'gray'  # #808080
ecolor = '#D3D3D3'
DOSsubplot2(dos_dict[study],2,color,ecolor)

study = 'hcpFeNi'
color = 'darkorange' # #FF8C00
ecolor = '#FFD4A6'
DOSsubplot2(dos_dict[study],1,color,ecolor)

study = 'hcpFeNiSi'
color = 'deepskyblue'  # #00BFFF
ecolor = '#A4ECFF'
DOSsubplot2(dos_dict[study],0,color,ecolor)

ax.set_xlabel(r'Energy (meV)',fontsize = 16)
ax.set_ylabel(r'PDOS $D(E,V)$',fontsize = 16)
ax.set_xlim([0,90])
# ax.set_ylim(ymin=-10,ymax=1100)
ax.xaxis.set_ticks([0,20,40,60,80])
ax.tick_params(direction='in',left='off',top='on')
ax.set_yticklabels([])

fig.savefig('PDOS_Pres2.pdf', format='pdf', bbox_inches='tight')
plt.close()