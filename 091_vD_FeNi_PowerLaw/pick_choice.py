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

vD_results_path = 'Results/vD_results_Constrained_Power.csv'
fits_filename = 'fit_search_results_Constrained_Power.csv'


# Functions
###########


# Load data
###########

results_df = pd.read_csv(vD_results_path)
folder_list = results_df['Folder'].values


# Get best picks from results_df
################################

for folder in folder_list:
	index = results_df[results_df['Folder']==folder]['Index'].iloc[0]
	
	# Load df with all possible fits
	fits_df = pd.read_csv(folder+'/'+fits_filename)
	
	# Get Debye velocity and uncertainty from pdf analysis
	vD = results_df[results_df['Folder']==folder]['vD'].iloc[0]
	dvD = np.abs(results_df[results_df['Folder']==folder]['dvD'].iloc[0])

	# Get range of vD results that we want our selected fit to agree with
	min_vD = vD - dvD/4
	max_vD = vD + dvD/4

	# Pick the fits that agree with our pdf analysis
	bestfits_df = fits_df[fits_df['vD'] >= min_vD]
	bestfits_df = bestfits_df[bestfits_df['vD'] <= max_vD]
	
	# Make the minimum energy range fairly small
	Emin_min = 3.0
	Emin_max = 5.0
	bestfits_df = bestfits_df[bestfits_df['Emin'] >= Emin_min]
	bestfits_df = bestfits_df[bestfits_df['Emin'] <= Emin_max]

	# Pick the fits with the best probability based on AICc
	bestfits_df = bestfits_df.nlargest(n=90,columns=['Prob'])

	print(folder+'  '+index+':')
	print(bestfits_df)
	print('\n')