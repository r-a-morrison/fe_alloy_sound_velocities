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

phoxpath = '../../050_phox_Fe_man/'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis, and the velocity results will be useful too for 
# calculating the Debye Gruneisen parameter:
EOSpath = '../../081_vD_Fe_PowerLaw/Results/input_values.csv'
velpath = '../../081_vD_Fe_PowerLaw/Results/vD_results_Constrained_Power.csv'


# Functions
###########

def plotdos(E,dos,ddos,filename):
	fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(10,3))

	ax0.plot(E,dos)
	ax0.fill_between(E, dos-ddos, dos+ddos, alpha=0.3)
	ax0.set_xlabel(r'Energy (meV)',fontsize = 16)
	ax0.set_ylabel(r'PDOS $D(E,V)$',fontsize = 16)
	ax0.set_xlim(xmin=0)
	ax0.set_ylim(ymin=-10)
	ax0.tick_params(direction='in',right='on',top='on')

	fig.savefig(filename, format='pdf')
	plt.close()


def plot2dos(E,dos,ddos,E_ref,dos_ref,filename):
	fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(10,3))

	ax0.plot(E,dos)
	ax0.fill_between(E, dos-ddos, dos+ddos, alpha=0.3)
	ax0.plot(E_ref,dos_ref,color='red')
	ax0.set_xlabel(r'Energy (meV)',fontsize = 16)
	ax0.set_ylabel(r'PDOS $D(E,V)$',fontsize = 16)
	ax0.set_xlim(xmin=0)
	ax0.set_ylim(ymin=-10)
	ax0.tick_params(direction='in',right='on',top='on')

	fig.savefig(filename, format='pdf',
		bbox_inches='tight')
	plt.close()


def plotMisfit(xi,chi2,prob,ref_folder,folder):

	fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,4))
	ax0.plot(xi,chi2)
	ax0.set_xlabel(r'$\xi$',fontsize = 16)
	ax0.set_ylabel(r'$\chi^2$',fontsize = 16)
	ax0.tick_params(direction='in',right='on',top='on')
	fig.savefig(ref_folder+'/chi2/scaled2_'+folder+'_chi2.pdf', format='pdf',
		bbox_inches='tight')
	plt.close()
	
	fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,4))
	ax0.plot(xi,prob)
	ax0.set_xlabel(r'$\xi$',fontsize = 16)
	ax0.set_ylabel(r'PDF',fontsize = 16)
	ax0.tick_params(direction='in',right='on',top='on')
	fig.savefig(ref_folder+'/pdf/scaled2_'+folder+'_pdf.pdf', format='pdf',
		bbox_inches='tight')
	plt.close()

def scaleDOS(ref_folder,folder,dos_dict):
	dos_ref_df = dos_dict[ref_folder]

	# Interpolate the reference PDOS
	fdos = interp1d(dos_ref_df['E'], dos_ref_df['DOS'], kind='cubic')
	E_ref_min = min(dos_ref_df['E'])
	E_ref_max = max(dos_ref_df['E'])

	dos_df = dos_dict[folder]

	xi_range = np.arange(0.4,1.6,0.0005)
	chi2 = np.empty(np.size(xi_range))
	prob = np.empty(np.size(xi_range))

	# Grid search over scaling parameter xi
	for i, xi in enumerate(xi_range):
		# If xi > 1, we need to limit the energy range so we don't call the interpolated
		#	reference function out of range
		E_min = max(E_ref_min/xi,min(dos_df['E']))
		E_max = min(E_ref_max/xi,max(dos_df['E']))
		dos_crop_df = dos_df[dos_df['E']>(E_min)]
		dos_crop_df = dos_crop_df[dos_crop_df['E']<(E_max)]
		# The first zero messes up the misfit calculation later on
		dos_crop_df = dos_crop_df[dos_crop_df['E']!=0]

		# Scale the reference PDOS
		dos_scaled = xi*fdos(xi*dos_crop_df['E'])

		# Misfit
		# chi2[i] = sum(np.divide((dos_crop_df['DOS']-dos_scaled)**2, dos_crop_df['dDOS']**2))
		chi2[i] = sum((dos_crop_df['DOS']-dos_scaled)**2)
		# Convert to reduced chi^2
		n = len(dos_crop_df) # number of data points
		m = 1 # number of parameters
		chi2[i] = chi2[i]/(n-m)

	# Create a pdf for xi
	prob = np.exp(-chi2)
	# Normalize pdf
	xistepsize = (max(xi_range)-min(xi_range))/len(xi_range)
	normconst = 1/(np.sum(prob)*xistepsize)
	prob = normconst*prob

	# Create plots of the misfit
	plotMisfit(xi_range,chi2,prob,ref_folder,folder)

	# Set to the best fit xi
	xi = xi_range[chi2.argmin()]

	# The standard deviation from the pdf is only correct if uncertainty is accounted
	# for in the optimization of xi. This produces an uncertainty on xi of ~0.02 for all
	# xi. However, including the PDOS uncertainty causes the optimization to perform poorly
	# for some compression points. Therefore, the optimization is performed without the
	# PDOS uncertainty, and an uncertainty on xi of 0.02 is assumed. Therefore, the output
	# pdf plots are not representative of the true uncertainty on xi.

	dxi = 0.02 

	# # The following uncerainty determination methods agree well but are unreasonably small
	# # if PDOS uncertainty is neglected
	# Calculate std dev of xi
	# mean = np.average(xi_range, weights=prob)
	# var = np.average((xi_range-xi)**2, weights=prob)
	# xi_stddev = np.sqrt(var)
	# # Fit a Gaussian to prob peak to estimate uncertainty
	# Gauss_opt,Gauss_cov = curve_fit(Gaussian,xi_range,prob,p0=[1,mean,stddev])
	# print(Gauss_opt)

	# Scale PDOS using the best fit xi
	# If xi > 1, we need to limit the energy range so we don't call the interpolated
	#	reference function out of range
	E_min = max(E_ref_min/xi,min(dos_df['E']))
	E_max = min(E_ref_max/xi,max(dos_df['E']))
	dos_crop_df = dos_df[dos_df['E']>(E_min)]
	dos_crop_df = dos_crop_df[dos_crop_df['E']<(E_max)]
	dos_scaled = xi*fdos(xi*dos_crop_df['E'])

	# Create a plot comparing the PDOS
	plot2dos(dos_crop_df['E'],dos_crop_df['DOS'],dos_crop_df['dDOS'],
		dos_crop_df['E'],dos_scaled,ref_folder+'/scaledDOSplot/scaled2_'+folder+'.pdf')

	# Save scaled PDOS
	dos_scaled_df = dos_crop_df[['E']].copy()
	dos_scaled_df['DOS'] = dos_scaled
	dos_scaled_df.to_csv(ref_folder+'/scaledDOSdata/scaled2_'+folder+'.csv',index=False)

	return(xi, dxi, dos_scaled_df)


# Gaussian function for fitting
def Gaussian(x,amplitude,mean,sigma):
	return amplitude*np.exp(-(x-mean)**2/(2*sigma**2))

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
		
def plotScaledPDOS(ref_folder,folder,dos_dict,dos_scaled_dict):
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

	fig.savefig(ref_folder+'/scaledPDOS_Fe_narrow.pdf', format='pdf', bbox_inches='tight')
	plt.close()


# Import data
#############

# Load XRD, bulk modulus, and sound velocity results
EOS_df = pd.read_csv(EOSpath, engine='python')
vel_df = pd.read_csv(velpath, engine='python')

input_df = EOS_df.merge(vel_df, on=('Folder','Index','vphi','dvphi'))
input_df.to_csv('Results/input_values.csv',index=False)

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

# Scale PDOS
############

# Create a dataframe to store results in
results_df = pd.DataFrame(columns = ['Ref Folder','Ref Index','Vi','dVi',
									 'Folder','Index','V','dV','V/Vi','xi', 'dxi'])

# for ref_folder in ['2009Oct_30GPa']:
for ref_folder in folder_list:
	print('Reference PDOS: '+ref_folder)

	# Check if a folder for the reference PDOS exists, and make one if not
	if not os.path.exists(ref_folder):
		   os.makedirs(ref_folder)
	# Check if subdirectories exist, and make them if not
	if not os.path.exists(ref_folder+'/chi2'):
		   os.makedirs(ref_folder+'/chi2')
	if not os.path.exists(ref_folder+'/pdf'):
		   os.makedirs(ref_folder+'/pdf')
	if not os.path.exists(ref_folder+'/scaledDOSdata'):
		   os.makedirs(ref_folder+'/scaledDOSdata')
	if not os.path.exists(ref_folder+'/scaledDOSplot'):
		   os.makedirs(ref_folder+'/scaledDOSplot')

	dos_ref_df = dos_dict[ref_folder]
	plotdos(dos_ref_df['E'],dos_ref_df['DOS'],dos_ref_df['dDOS'],
		ref_folder+'/'+ref_folder+'_pdos.pdf')
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

		xi, dxi, dos_scaled_dict[folder] = scaleDOS(ref_folder,folder,dos_dict)

		results_df = results_df.append(pd.DataFrame([[
			ref_folder,index_dict[ref_folder],Vi,dVi,
			folder,index_dict[folder],V,dV,V_Vi,xi,dxi]],columns = results_df.columns))

	# Create plot of ref PDOS scaled to all other PDOS
	plotScaledPDOS(ref_folder,folder,dos_dict,dos_scaled_dict)

	# At this point in the nested loops to save results in case code crashes
	# Will be overwritten on each loop with updated results
	results_df = results_df.round({'V/Vi':3,'xi':4,'dxi':4})
	results_df.to_csv('Results/scalingparameters.csv',index=False)

	