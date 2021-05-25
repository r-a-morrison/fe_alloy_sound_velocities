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

# Define V0 (from EOS)
# Dewaele et al. 2006 EOS
V0 = 22.952
dV0 = 0.072

# Values from previous analyses, e.g. volume and vD
input_filename = 'Results/scalingparameters.csv'

# Calculate estimated Gruneisen parameter
#########################################

input_df = pd.read_csv(input_filename)

input_df['gamma_est'] = (1-1/input_df['xi'])*input_df['Vi']/(input_df['V']-input_df['Vi'])
input_df = input_df.round(4)
input_df.to_csv('Results/GruneisenResults_method2_intermediate.csv',index=False)

# Remove outliers
#################
# V=19.31 is very close to V=19.09, so the derivative dE/dV is poorly constrained.
# The result of gamma_est from these two PDOS scaled to each other is a clear outlier
# and should be removed
input_df = input_df[~((input_df['Vi'] == 19.087) & (input_df['V'] == 19.305))]
input_df = input_df[~((input_df['Vi'] == 19.305) & (input_df['V'] == 19.087))]
print(input_df)

# Average estimated Gruneisen parameter for each reference volume
#################################################################

gamma_est_df = input_df.groupby('Vi')['gamma_est'].mean()
gamma_est_df = pd.DataFrame(gamma_est_df).reset_index()
dgamma_est_df = input_df.groupby('Vi')['gamma_est'].std()
dgamma_est_df = pd.DataFrame(dgamma_est_df).reset_index()
dgamma_est_df.columns = ['Vi','dgamma_est']
gamma_est_df = gamma_est_df.merge(dgamma_est_df,on='Vi')
gamma_est_df = gamma_est_df.round(4)
gamma_est_df.to_csv('Results/GruneisenResults_method2.csv',index=False)