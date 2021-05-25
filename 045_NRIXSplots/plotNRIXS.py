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

matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

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

phoxpath = dict()
EOSpath  = dict()
velpath  = dict()

# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis:
phoxpath['Fe'] = '../050_phox_Fe_man/'
EOSpath['Fe'] = '../081_vD_Fe_PowerLaw/Results/input_values.csv'
velpath['Fe'] = '../081_vD_Fe_PowerLaw/Results/vD_results_Constrained_Power.csv'

phoxpath['FeNi'] = '../060_phox_FeNi_man/'
EOSpath['FeNi'] = '../091_vD_FeNi_PowerLaw/Results/input_values.csv'
velpath['FeNi'] = '../091_vD_FeNi_PowerLaw/Results/vD_results_Constrained_Power.csv'

phoxpath['FeNiSi'] = '../070_phox_FeNiSi_man/'
EOSpath['FeNiSi'] = '../101_vD_FeNiSi_PowerLaw/Results/input_values.csv'
velpath['FeNiSi'] = '../101_vD_FeNiSi_PowerLaw/Results/vD_results_Constrained_Power.csv'


# Functions
###########



# Import data
#############

input_df_dict = dict()  # Stores volume, P, etc. info
folders_dict = dict()   # Stores folder that NRIXS data is stored in
index_dict = dict()     # Stores file name of NRIXS data (nested dict)
NRIXS_dict = dict() # Stores dictionary of NRIXS dataframes (nested dict)

for study in phoxpath.keys(): # Iterate over each composition
	
	# Load XRD, bulk modulus, and sound velocity results
	EOS_df = pd.read_csv(EOSpath[study], engine='python')
	vel_df = pd.read_csv(velpath[study], engine='python')

	input_df = EOS_df.merge(vel_df, on=('Folder','Index','vphi','dvphi'))
	input_df_dict[study] = input_df

	# Load NRIXS data
	# Find the filepath of all .res NRIXS files in phox directories
	respath_list = [filepath for filepath in glob.glob(phoxpath[study]+'*/*.res')]

	# Prep lists, dictionaries, and df to store input data in
	folders = []
	index = dict()
	NRIXS_dfs = dict()

	# Collect folders, indices, paths, and input values
	for respath in respath_list:

		# Determine filepaths for NRIXS and in_psvl
		folder = re.findall('([A-Za-z0-9_]+)/[A-Za-z0-9_]+.res',respath)[0]
		index[folder] = re.findall('/([A-Za-z0-9_]+).res',respath)[0]
		NRIXSpath = phoxpath[study]+folder+'/'+index[folder]+'.dat'

		# Check if each folder is hcp. Don't use it otherwise
		phase = input_df[input_df['Folder']==folder].iloc[-1]['Phase']

		# Import NRIXS
		NRIXS_dfs[folder] = pd.read_csv(NRIXSpath, sep='\s+', comment='#', header=None,
			names = ['E','NRIXS','dNRIXS'])

		# Store to use NRIXS later
		folders.append(folder)

	folders_dict[study] = folders
	index_dict[study] = index
	NRIXS_dict[study] = NRIXS_dfs


# make NRIXS plot
################

# # Option 5
# def NRIXSsubplot(folder,number,color):
# 	offset = 2
# 	NRIXS_dfs = NRIXS_dict[study]
# 	NRIXS_df = NRIXS_dfs[folder]
# 	# ax.plot(NRIXS_df['E'], np.log10(NRIXS_df['NRIXS'])+offset*number,marker='.',
# 	# 	markersize=2,linewidth=0.5,color=color,linestyle='none')
# 	# ax.fill_between(NRIXS_df['E'],
# 	# 	np.log10(NRIXS_df['NRIXS']-NRIXS_df['dNRIXS'])+offset*number,
# 	# 	np.log10(NRIXS_df['NRIXS']+NRIXS_df['dNRIXS'])+offset*number, 
#	#	facecolor=color, alpha=0.3)
# 	ax.errorbar(NRIXS_df['E'], NRIXS_df['NRIXS']*10**(offset*number),
# 		yerr=NRIXS_df['dNRIXS']*10**(offset*number),marker='.',
# 		markersize=1.5,linestyle='none',color=color,ecolor='lightgray',elinewidth=1,
# 		zorder=-5)
# 	ax.set_yscale('log')
	
# fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(3.5,10))

# study = 'FeNi'

# color = 'black'
# NRIXSsubplot('2012Oct_DAC11_P14',7,color)
# NRIXSsubplot('2012Apr_DAC11_P9',6,color)
# NRIXSsubplot('2012Apr_DAC11_P8',5,color)
# NRIXSsubplot('2012Apr_DAC11_P7',4,color)
# NRIXSsubplot('2012Apr_DAC11_P6',3,color)
# NRIXSsubplot('2012Apr_DAC11_P5',2,color)
# NRIXSsubplot('2012Apr_DAC11_P4',1,color)
# NRIXSsubplot('2012Apr_DAC11_P3',0,color)
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


# Option 4
def NRIXSsubplot3(ax,folder,number,color,ymin):
	offset = 5
	NRIXS_dfs = NRIXS_dict[study]
	NRIXS_df = NRIXS_dfs[folder]
	ax.errorbar(NRIXS_df['E'], NRIXS_df['NRIXS'],yerr=NRIXS_df['dNRIXS'],marker='.',
		markersize=1,linestyle='none',color=color,ecolor='lightgray',elinewidth=0.75,
		zorder=-5)
	ax.set_yscale('log')
	ax.tick_params(direction='in',right='on',bottom='off')
	ax.set_xlim([-85,105])
	ax.set_ylim(ymin=ymin,ymax=5*10**3)
	ax.yaxis.set_ticks([1, 100])
	ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator())
	ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
	ax.tick_params(which='minor',direction='in',right='on')
	ax.tick_params(axis='both', which='major', labelsize=11)
	# ax.minorticks_on()

fig, ((ax7a,ax7b),(ax6a,ax6b),(ax5a,ax5b),(ax4a,ax4b),(ax3a,ax3b),(ax2a,ax2b),
	(ax1a,ax1b),(ax0a,ax0b)) = plt.subplots(nrows = 8,ncols=2, figsize=(6,10))

# color = '#1f77b4'
color = 'black'
ymin = 0.5

study = 'FeNi'
NRIXSsubplot3(ax7a,'2012Oct_DAC11_P14', 7,color,ymin)
NRIXSsubplot3(ax6a,'2012Apr_DAC11_P9',  6,color,ymin)
NRIXSsubplot3(ax5a,'2012Apr_DAC11_P8',  5,color,ymin)
NRIXSsubplot3(ax4a,'2012Apr_DAC11_P7',  4,color,ymin)
NRIXSsubplot3(ax3a,'2012Apr_DAC11_P6',  3,color,ymin)
NRIXSsubplot3(ax2a,'2012Apr_DAC11_P5',  2,color,ymin)
NRIXSsubplot3(ax1a,'2012Apr_DAC11_P4',  1,color,ymin)
NRIXSsubplot3(ax0a,'2012Apr_DAC11_P3',  0,color,ymin)

study = 'FeNiSi'
NRIXSsubplot3(ax5b,'2015Mar_DAC13_P6',  5,color,ymin)
NRIXSsubplot3(ax4b,'2015Mar_DAC13_P5',  4,color,ymin)
NRIXSsubplot3(ax3b,'2015Mar_DAC13_P4',  3,color,ymin)
NRIXSsubplot3(ax2b,'2015Mar_DAC13_P3',  2,color,ymin)
NRIXSsubplot3(ax1b,'2014Feb_DAC13_P2',  1,color,ymin)
NRIXSsubplot3(ax0b,'2014Feb_DAC13_P1',  0,color,ymin)

xaxis_ticks = [-80,-40,0,40,80]

ax7a.xaxis.set_ticks(xaxis_ticks)
ax7a.xaxis.set_minor_locator(AutoMinorLocator(2))
ax7a.tick_params(direction='in',top='on')
ax7a.tick_params(which='minor',direction='in',top='on',bottom='off')
ax7a.set_xticklabels([])
ax6a.set_xticklabels([])
ax5a.set_xticklabels([])
ax4a.set_xticklabels([])
ax3a.set_xticklabels([])
ax2a.set_xticklabels([])
ax1a.set_xticklabels([])
ax0a.tick_params(direction='in',bottom='on')
ax0a.xaxis.set_ticks(xaxis_ticks)
ax0a.xaxis.set_minor_locator(AutoMinorLocator(2))
ax0a.set_xlabel(r'Energy (meV)',fontsize = 14)
ax3a.set_ylabel(r'NRIXS (counts)',fontsize = 14)

ax7b.set_frame_on(False)
ax7b.tick_params(direction='in',left='off',bottom='off')
ax7b.set_xticklabels([])
ax7b.set_yticklabels([])
ax6b.set_frame_on(False)
ax6b.tick_params(direction='in',left='off',bottom='off')
ax6b.set_xticklabels([])
ax6b.set_yticklabels([])
ax5b.tick_params(direction='in',top='on')
ax5b.tick_params(which='minor',direction='in',top='on',bottom='off')
ax5b.xaxis.set_minor_locator(AutoMinorLocator(2))
ax5b.xaxis.set_ticks(xaxis_ticks)
ax5b.set_xticklabels([])
ax4b.set_xticklabels([])
ax3b.set_xticklabels([])
ax3b.set_ylabel(r'NRIXS (counts)',fontsize = 14)
ax2b.set_xticklabels([])
ax1b.set_xticklabels([])
ax0b.tick_params(direction='in',bottom='on')
ax0b.xaxis.set_ticks(xaxis_ticks)
ax0b.xaxis.set_minor_locator(AutoMinorLocator(2))
ax0b.set_xlabel(r'Energy (meV)',fontsize = 14)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(hspace=0,wspace=0.4)

fig.savefig('NRIXShcp.pdf', format='pdf', bbox_inches='tight')
plt.close()



def NRIXSsubplot3(ax,folder,number,color,ymin):
	offset = 5
	NRIXS_dfs = NRIXS_dict[study]
	NRIXS_df = NRIXS_dfs[folder]
	ax.errorbar(NRIXS_df['E'], NRIXS_df['NRIXS'],yerr=NRIXS_df['dNRIXS'],marker='.',
		markersize=1,linestyle='none',color=color,ecolor='lightgray',elinewidth=0.75,
		zorder=-5)
	ax.set_yscale('log')
	ax.tick_params(direction='in',right='on',bottom='off')
	ax.set_xlim([-85,105])
	ax.set_ylim(ymin=ymin,ymax=5*10**3)
	ax.yaxis.set_ticks([1, 100])
	ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator())
	ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
	ax.tick_params(which='minor',direction='in',right='on')
	ax.tick_params(axis='both', which='major', labelsize=11)
	
fig, ((ax4a,ax4b),(ax3a,ax3b),(ax2a,ax2b),
	(ax1a,ax1b),(ax0a,ax0b)) = plt.subplots(nrows = 5,ncols=2, figsize=(6,8))

# color = '#1f77b4'
color = 'black'
ymin = 0.5

study = 'FeNi'
NRIXSsubplot3(ax4a,'2012Apr_DAC11_P2',  4,color,ymin)
NRIXSsubplot3(ax3a,'2013Feb_DAC14_P2',  3,color,ymin)
NRIXSsubplot3(ax2a,'2012Apr_DAC11_P1',  2,color,ymin)
NRIXSsubplot3(ax1a,'2013Feb_DAC14_P1',  1,color,ymin)
NRIXSsubplot3(ax0a,'2015Mar_Ambient',  0,color,ymin)

study = 'FeNiSi'
NRIXSsubplot3(ax2b,'2013Nov_DAC14_P1',  2,color,ymin)
NRIXSsubplot3(ax1b,'2014Feb_DAC15_P1',  1,color,ymin)
NRIXSsubplot3(ax0b,'2015Mar_Ambient',  0,color,ymin)

xaxis_ticks = [-80,-40,0,40,80]

ax4a.tick_params(direction='in',top='on')
ax4a.xaxis.set_ticks(xaxis_ticks)
ax4a.tick_params(which='minor',direction='in',top='on',bottom='off')
ax4a.xaxis.set_minor_locator(AutoMinorLocator(2))
ax4a.set_xticklabels([])
ax3a.set_xticklabels([])
ax2a.set_xticklabels([])
ax1a.set_xticklabels([])
ax0a.tick_params(direction='in',bottom='on')
ax0a.xaxis.set_ticks(xaxis_ticks)
ax0a.xaxis.set_minor_locator(AutoMinorLocator(2))
ax0a.set_xlabel(r'Energy (meV)',fontsize = 14)
ax2a.set_ylabel(r'NRIXS (counts)',fontsize = 14)

ax4b.set_frame_on(False)
ax4b.tick_params(direction='in',left='off',bottom='off')
ax4b.set_xticklabels([])
ax4b.set_yticklabels([])
ax3b.set_frame_on(False)
ax3b.tick_params(direction='in',left='off',bottom='off')
ax3b.set_xticklabels([])
ax3b.set_yticklabels([])
ax2b.tick_params(direction='in',top='on')
ax2b.xaxis.set_ticks(xaxis_ticks)
ax2b.tick_params(which='minor',direction='in',top='on',bottom='off')
ax2b.xaxis.set_minor_locator(AutoMinorLocator(2))
ax2b.set_xticklabels([])
ax1b.set_xticklabels([])
ax1b.set_ylabel(r'NRIXS (counts)',fontsize = 14)
ax0b.tick_params(direction='in',bottom='on')
ax0b.xaxis.set_ticks(xaxis_ticks)
ax0b.xaxis.set_minor_locator(AutoMinorLocator(2))
ax0b.set_xlabel(r'Energy (meV)',fontsize = 14)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(hspace=0,wspace=0.4)

fig.savefig('NRIXSbcc.pdf', format='pdf', bbox_inches='tight')
plt.close()



def NRIXSsubplot3(ax,folder,number,color,ecolor,ymin):
	offset = 5
	NRIXS_dfs = NRIXS_dict[study]
	NRIXS_df = NRIXS_dfs[folder]
	ax.errorbar(NRIXS_df['E'], NRIXS_df['NRIXS'],yerr=NRIXS_df['dNRIXS'],marker='.',
		markersize=1,linestyle='none',color=color,ecolor=ecolor,elinewidth=0.75,
		zorder=-5)
	ax.set_yscale('log')
	ax.tick_params(direction='in',right='on',bottom='off')
	ax.set_xlim([-85,95])
	ax.set_ylim(ymin=ymin,ymax=5*10**3)
	ax.tick_params(axis='both', which='major', labelsize=11)
	ax.minorticks_off()

fig, ((ax2),(ax1),(ax0)) = plt.subplots(nrows = 3,ncols=1, figsize=(3,6))

# color = '#1f77b4'
ymin = 0.3

study = 'Fe'
color = 'gray'  # #808080
ecolor = '#D3D3D3'
NRIXSsubplot3(ax0,'2015Mar_Ambient',  0,color,ecolor,ymin)

study = 'FeNi'
color = 'darkorange' # #FF8C00
ecolor = '#FFD4A6'
NRIXSsubplot3(ax1,'2015Mar_Ambient',  0,color,ecolor,ymin)

study = 'FeNiSi'
color = 'deepskyblue'  # #00BFFF
ecolor = '#A4ECFF'
NRIXSsubplot3(ax2,'2015Mar_Ambient',  0,color,ecolor,ymin)

ax2.tick_params(direction='in',top='on')
ax2.tick_params(which='minor',direction='in',top='on',bottom='off')
ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
ax2.xaxis.set_ticks([-80,-40,0,40,80])
ax2.set_xticklabels([])
ax1.set_xticklabels([])
ax1.set_ylabel(r'NRIXS (counts)',fontsize = 14)
ax0.tick_params(direction='in',bottom='on')
ax0.xaxis.set_ticks([-80,-40,0,40,80])
ax0.xaxis.set_minor_locator(AutoMinorLocator(2))
ax0.set_xlabel(r'Energy (meV)',fontsize = 14)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(hspace=0,wspace=0.4)

fig.savefig('NRIXSambient.pdf', format='pdf', bbox_inches='tight')
plt.close()



# fig, ax = plt.subplots(nrows = 1,ncols=1, figsize=(2.5,2))

# # color = '#1f77b4'
# color = 'black'
# ymin = 3
# study = 'FeNi'
# folder = '2012Apr_DAC11_P6'
# index = 'FeNi_DAC11_P6'

# NRIXS_dfs = NRIXS_dict[study]
# NRIXS_df = NRIXS_dfs[folder]

# norm_df = pd.read_csv(phoxpath['FeNi']+folder+'/Output/'+index+'_psn.dat',
# 			sep='\s+', comment='@', header=None, names = ['E','NRIXS','dNRIXS'])
# # res_df = pd.read_csv(phoxpath['FeNi']+folder+'/'+index+'.res',
# # 			sep='\s+', comment='#', header=None, names = ['E','NFS','dNFS'])
# # res_df = pd.read_csv(phoxpath['FeNi']+folder+'/Output/'+index+'_rfc.dat',
# # 			sep='\s+', comment='@', header=None, names = ['E','NFS'])
# ph1_df = pd.read_csv(phoxpath['FeNi']+folder+'/Output/'+index+'_1ph.dat',
# 			sep='\s+', comment='@', header=None, names = ['E','1phonon'])
# ph2_df = pd.read_csv(phoxpath['FeNi']+folder+'/Output/'+index+'_2ph.dat',
# 			sep='\s+', comment='@', header=None, names = ['E','2phonon'])
# ph3_df = pd.read_csv(phoxpath['FeNi']+folder+'/Output/'+index+'_3ph.dat',
# 			sep='\s+', comment='@', header=None, names = ['E','3phonon'])

# ax.errorbar(norm_df['E'], norm_df['NRIXS'],yerr=norm_df['dNRIXS'],marker='.',
# 	markersize=1,linestyle='none',color=color,ecolor='lightgray',elinewidth=0.75,
# 	zorder=-5)
# # ax.plot(res_df['E'],res_df['NFS'],'red')
# ax.plot(ph1_df['E'],ph1_df['1phonon'],'blue')
# ax.plot(ph2_df['E'],ph2_df['2phonon'],'green')
# ax.plot(ph3_df['E'],ph3_df['3phonon'],'orange')
# ax.set_yscale('log')
# ax.set_xlim([-80,110])
# # ax.set_ylim(ymin=ymin,ymax=5*10**3)
# ax.tick_params(axis='both', which='major', labelsize=11)
# ax.minorticks_off()

# ax.tick_params(direction='in',right='on',bottom='on',top='on')
# ax.xaxis.set_ticks([-60,-40,-20,0,20,40,60,80])
# ax.set_xlabel(r'Energy (meV)',fontsize = 14)
# ax.set_ylabel(r'NRIXS (counts)',fontsize = 14)

# # ax8.set_ylabel(r'NRIXS (counts)',fontsize = 14)

# # Fine-tune figure; make subplots close to each other and hide x ticks for
# # all but bottom plot.
# fig.subplots_adjust(hspace=0,wspace=0.4)
# # plt.setp([ax0.get_xticklabels() for a in fig.axes[:1]], visible=False);

# fig.savefig('NRIXS_PhononModes.pdf', format='pdf', bbox_inches='tight')
# plt.close()