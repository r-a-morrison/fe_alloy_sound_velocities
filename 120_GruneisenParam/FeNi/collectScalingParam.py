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


# Input scaling parameter results
##########################################

xi_filename = 'Results/scalingparameters.csv'
xi_df = pd.read_csv(xi_filename)

# Rename columns to avoid confusion
xi_df = xi_df.rename(columns={'Vi':'Vj', 'dVi':'dVj', 'V':'Vk','dV':'dVk',
	'V/Vi':'Vk/Vj','xi':'xi(Vk/Vj)','dxi':'dxi(Vk/Vj)'})

# Transform scaling parameters to each reference volume
#######################################################

folder_list = xi_df.drop_duplicates(subset='Ref Folder')['Ref Folder'].values

for ref_folder in folder_list:
# for ref_folder in ['2009Oct_30GPa']:
	print('Rescaling to '+ref_folder)
	# Reference volume to scale everything to
	Vi = xi_df[xi_df['Ref Folder']==ref_folder].iloc[-1]['Vj']

	xi_rescaled_df = xi_df[['Vj','Vk','xi(Vk/Vj)','dxi(Vk/Vj)']].copy()
	xi_rescaled_df['Vi'] = Vi*np.ones(len(xi_rescaled_df))

	# rescaled xi(Vk/Vi) = xi(Vk/Vj) * complementary xi(Vj/Vi)
	# Complementary xi needed to calculate rescaled xi:
	xi_rescaled_df['xi(Vj/Vi)'] = [xi_rescaled_df[(xi_rescaled_df['Vj']==Vi) &
		(xi_rescaled_df['Vk']==Vj)].iloc[-1]['xi(Vk/Vj)'] for Vj in xi_rescaled_df['Vj']]
	xi_rescaled_df['dxi(Vj/Vi)'] = [xi_rescaled_df[(xi_rescaled_df['Vj']==Vi) &
		(xi_rescaled_df['Vk']==Vj)].iloc[-1]['dxi(Vk/Vj)'] for Vj in xi_rescaled_df['Vj']]
	
	xi_rescaled_df['Vk/Vi'] = xi_rescaled_df['Vk']/xi_rescaled_df['Vi']

	# Calculate rescaled xi
	xi_rescaled_df['xi(Vk/Vi)'] = xi_rescaled_df['xi(Vk/Vj)']*xi_rescaled_df['xi(Vj/Vi)']
	# Calculate uncertainty on rescaled xi
	#	If c = a*b, dc = sqrt((b*da)^2 + (a*db)^2)
	xi_rescaled_df['dxi(Vk/Vi)'] = np.sqrt(
		(xi_rescaled_df['xi(Vj/Vi)']*xi_rescaled_df['dxi(Vk/Vj)'])**2 +
		(xi_rescaled_df['xi(Vk/Vj)']*xi_rescaled_df['dxi(Vj/Vi)'])**2)

	# Eliminate data points where Vi = Vk
	xi_rescaled_df = xi_rescaled_df[xi_rescaled_df['Vk'] != Vi]
	
	xi_rescaled_df = xi_rescaled_df.round(decimals=4)
	xi_rescaled_df.to_csv(ref_folder+'/rescaledparameters.csv',index=False)

	# Plot scaling parameters
	fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))

	ax0.errorbar(xi_rescaled_df['Vk/Vi'],xi_rescaled_df['xi(Vk/Vi)'],
		yerr=xi_rescaled_df['dxi(Vk/Vi)'],
		marker = 'o', color = 'darkorange', mfc='#ffc681', ms=6, markeredgewidth=1,
		ls='none',elinewidth=1)
	ax0.set_xlabel(r'$V/V_i$',fontsize = 16)
	ax0.set_ylabel(r'$\xi$',fontsize = 16)
	ax0.tick_params(direction='in',right='on',top='on')

	fig.savefig(ref_folder+'/scalingparam.pdf', format='pdf',
		bbox_inches='tight')
	plt.close()