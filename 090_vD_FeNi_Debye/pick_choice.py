# Substitute for psvl

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
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
from scipy.interpolate import spline

# Seaborn, useful for graphics
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

V_D_results_path = 'Results/V_D_results_Debye.csv'
fits_filename = 'fit_search_results_Debye.csv'


# Functions
###########


# Load data
###########

results_df = pd.read_csv(V_D_results_path)
folder_list = results_df['Folder'].values


# Get best picks from results_df
################################

for folder in folder_list:
	index = results_df[results_df['Folder']==folder]['Index'].iloc[0]
	
	# Load df with all possible fits
	fits_df = pd.read_csv(folder+'/'+fits_filename)
	
	# Get Debye velocity and uncertainty from pdf analysis
	V_D = results_df[results_df['Folder']==folder]['V_D'].iloc[0]
	dV_D = np.abs(results_df[results_df['Folder']==folder]['dV_D'].iloc[0])

	# Get range of V_D results that we want our selected fit to agree with
	min_V_D = V_D - dV_D/4
	max_V_D = V_D + dV_D/4

	# Pick the fits that agree with our pdf analysis
	bestfits_df = fits_df[fits_df['V_D'] >= min_V_D]
	bestfits_df = bestfits_df[bestfits_df['V_D'] <= max_V_D]
	
	# Make the minimum energy range fairly small
	Emin_min = 3.0
	Emin_max = 5.0
	bestfits_df = bestfits_df[bestfits_df['Emin'] >= Emin_min]
	bestfits_df = bestfits_df[bestfits_df['Emin'] <= Emin_max]

	# Pick the fits with the best probability based on AICc
	bestfits_df = bestfits_df.nlargest(n=30,columns=['Prob'])

	print(folder+'  '+index+':')
	print(bestfits_df)
	print('\n')