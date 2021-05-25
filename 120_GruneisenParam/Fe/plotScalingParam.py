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


# Input scaling parameter results and plot
##########################################

xi_filename = 'Results/scalingparameters.csv'

xi_df = pd.read_csv(xi_filename)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))

# ax0.plot(xi_df['V/Vi'],xi_df['xi'],marker = 'o')
ax0.errorbar(xi_df['V/Vi'],xi_df['xi'],yerr=xi_df['dxi'],marker = 'o',
	color = 'gray', mfc='lightgray', ms=5, markeredgewidth=1,ls='none',elinewidth=1)
ax0.set_xlabel(r'$V/V_i$',fontsize = 16)
ax0.set_ylabel(r'$\xi$',fontsize = 16)
# ax0.set_xlim(xmin=0)
# ax0.set_ylim(ymin=-10)
ax0.tick_params(direction='in',right='on',top='on')

fig.savefig('Results/scalingparam.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Check fit results
###################

V0 = 22.428

# model for xi
def xiModel(V,Vi,V0,gamma0,q):
	return np.exp( gamma0*(Vi/V0)**q * (1/q) * ((V/Vi)**q-1) )

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))

gamma0 = 2.04
q = 1.0
xi_check = xiModel(xi_df['V'],xi_df['Vi'],V0,gamma0,q)
ax0.plot(xi_df['V/Vi'],xi_check,'.',marker = 'o', ms=4)

gamma0 = 1.98
q = 1.0
xi_check = xiModel(xi_df['V'],xi_df['Vi'],V0,gamma0,q)
ax0.plot(xi_df['V/Vi'],xi_check,'.',marker = 'o', ms=4)

ax0.errorbar(xi_df['V/Vi'],xi_df['xi'],yerr=xi_df['dxi'],marker = 'o',
	color = 'gray', mfc='lightgray', ms=4, markeredgewidth=1,ls='none',elinewidth=1)
ax0.set_xlabel(r'$V/V_i$',fontsize = 16)
ax0.set_ylabel(r'$\xi$',fontsize = 16)
ax0.tick_params(direction='in',right='on',top='on')

fig.savefig('test.pdf', format='pdf',
	bbox_inches='tight')
plt.close()
