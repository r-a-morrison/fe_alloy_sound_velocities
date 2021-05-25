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


# Input information
######################

# Path to get volumes from hcp-Fe Murphy dataset
EOSpath = '../../081_vD_Fe_PowerLaw/Results/input_values.csv'
# Path to get re-analyzed hcp-Fe Murphy PDOS
phoxpath = '../../050_phox_Fe_man/'

# Fei Gruneisen parameter variables
gamma0 = 1.74
q = 0.78

# From Dewaele 2006, agrees with Fei 2016 to third decimal place.
V0 = 22.428
dV0 = 0.098
# Verification:
rho0 = 8.2695 # g/cc, Fei 2016 density
M = 55.845    # g/mol, for natural Fe
V0_ccpermol = M/rho0  # cm^3/mol
V0_check = V0_ccpermol*(2*10**24)/constants.N_A  # A^3
print(V0_check)


# Functions
###########

def calcScalingParam(V,Vi,V0,gamma0,q):
	xi = np.exp( gamma0*(Vi/V0)**q * (1/q) * ((V/Vi)**q-1) )
	return xi

def scaleDOS(ref_folder,xi,dos_dict):
	dos_ref_df = dos_dict[ref_folder]

	# Interpolate the reference PDOS
	fdos = interp1d(dos_ref_df['E'], dos_ref_df['DOS'], kind='cubic')
	E_ref_min = min(dos_ref_df['E'])
	E_ref_max = max(dos_ref_df['E'])

	# Scale PDOS using xi
	# If xi > 1, we need to limit the energy range so we don't call the interpolated
	#	reference function out of range
	E_min = max(E_ref_min/xi,min(dos_df['E']))
	E_max = min(E_ref_max/xi,max(dos_df['E']))
	dos_crop_df = dos_df[dos_df['E']>(E_min)]
	dos_crop_df = dos_crop_df[dos_crop_df['E']<(E_max)]
	dos_scaled = xi*fdos(xi*dos_crop_df['E'])

	# Save scaled PDOS
	dos_scaled_df = dos_crop_df[['E']].copy()
	dos_scaled_df['DOS'] = dos_scaled
	# dos_scaled_df.to_csv(ref_folder+'/scaledDOSdata/scaled2_'+folder+'.csv',index=False)

	return(dos_scaled_df)


def DOSsubplot(folder,ref_folder,number,ax):
	offset = 100
	dos_df = dos_dict[folder]
	# ax.plot(dos_df['E'], dos_df['DOS']+offset*number,color='#1f77b4')
	# ax.fill_between(dos_df['E'], dos_df['DOS']-dos_df['dDOS']+offset*number,
	# 	dos_df['DOS']+dos_df['dDOS']+offset*number, facecolor='#1f77b4', alpha=0.3)
	ax.errorbar(dos_df['E'], dos_df['DOS']+offset*number,yerr=dos_df['dDOS'],
		marker='.',markersize=1.5,color='black',ecolor='darkgray',elinewidth=0.5,
		linestyle='none',zorder=-5)
	if folder != ref_folder:
		dos_scaled_df = dos_scaled_dict[folder]
		ax.plot(dos_scaled_df['E'], dos_scaled_df['DOS']+offset*number,color='red')
		
def plotScaledPDOS(ref_folder,folder_list,dos_dict,dos_scaled_dict):
	fig, (ax) = plt.subplots(nrows = 1, ncols=1, figsize=(4,10))

	offsetnum = 0
	for folder in folder_list:
		DOSsubplot(folder,ref_folder,offsetnum,ax)
		offsetnum = offsetnum + 1

	ax.set_xlabel(r'Energy (meV)',fontsize = 16)
	ax.set_ylabel(r'PDOS $D(E,V)$',fontsize = 16)
	ax.set_xlim([0,85])
	ax.set_ylim(ymin=-10,ymax=1050)
	ax.xaxis.set_ticks([0,20,40,60,80])
	ax.tick_params(direction='in',left='off',top='on')
	ax.set_yticklabels([])

	fig.savefig(ref_folder+'/Fei_scaledPDOS_Fe_narrow.pdf', format='pdf', bbox_inches='tight')
	plt.close()


# Import data
#############

input_df = pd.read_csv(EOSpath, engine='python')
# # Only use hcp data (We measured the bcc phase, not Caitlin)
# input_df = input_df[input_df['Phase']=='hcp']
# # Data was out of order as imported. To fix that:
# input_df = input_df.sort_values('P')

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

	if phase == 'hcp':
		# Import PDOS
		dos_df = pd.read_csv(dospath, sep='\s+', comment='@', header=None,
			names = ['E','DOS','dDOS'])

		# Store to use PDOS later
		folder_list.append(folder)
		index_dict[folder] = index
		dos_dict[folder] = dos_df

# Sort folder_list by pressure (needed to plot in correct order)
sort_folder_df = pd.DataFrame(columns = ['Folder','P'])
for folder in folder_list:
	P = input_df[input_df['Folder']==folder].iloc[-1]['P']
	sort_folder_df = sort_folder_df.append(
		pd.DataFrame([[folder,P]],columns=sort_folder_df.columns))
sort_folder_df = sort_folder_df.sort_values('P')
folder_list = sort_folder_df['Folder'].values


# Plot scaled PDOS
##################

# Create a dataframe to store results in
results_df = pd.DataFrame(columns = ['Ref Folder','Ref Index','Vi','dVi',
									 'Folder','Index','V','dV','V/Vi','xi'])

# for ref_folder in ['2009Oct_30GPa']:
for ref_folder in folder_list:
	print('Reference PDOS: '+ref_folder)

	# Check if a folder for the reference PDOS exists, and make one if not
	if not os.path.exists(ref_folder):
		   os.makedirs(ref_folder)

	dos_ref_df = dos_dict[ref_folder]
	# What is the reference volume?
	Vi = input_df[input_df['Folder']==ref_folder].iloc[-1]['V']
	dVi = input_df[input_df['Folder']==ref_folder].iloc[-1]['dV']

	dos_scaled_dict = dict()

	# for folder in ['2011Feb_171GPa']:
	for folder in folder_list:
		print('\tScaling to '+folder)
		# What is the volume?
		V = input_df[input_df['Folder']==folder].iloc[-1]['V']
		dV = input_df[input_df['Folder']==folder].iloc[-1]['dV']
		V_Vi = V/Vi

		xi = calcScalingParam(V,Vi,V0,gamma0,q)

		dos_scaled_dict[folder] = scaleDOS(ref_folder,xi,dos_dict)

		results_df = results_df.append(pd.DataFrame([[
			ref_folder,index_dict[ref_folder],Vi,dVi,
			folder,index_dict[folder],V,dV,V_Vi,xi]],columns = results_df.columns))

	# Create plot of ref PDOS scaled to all other PDOS
	plotScaledPDOS(ref_folder,folder_list,dos_dict,dos_scaled_dict)

	# At this point in the nested loops, save results in case code crashes
	# Will be overwritten on each loop with updated results
	results_df = results_df.round({'V/Vi':3,'xi':4,'dxi':4})
	results_df.to_csv('Results/Fei_scalingparameters.csv',index=False)
