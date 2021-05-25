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

Fe_filename = 'Fe/Results/scalingparameters.csv'
FeNi_filename = 'FeNi/Results/scalingparameters.csv'
FeNiSi_filename = 'FeNiSi/Results/scalingparameters.csv'


fig, (ax0,ax1,ax2) = plt.subplots(nrows = 1, ncols=3, figsize=(10,4), sharey = True,sharex = True)

xi_df = pd.read_csv(Fe_filename)
ax0.errorbar(xi_df['V/Vi'],xi_df['xi'],yerr=xi_df['dxi'],marker = 'o',
	color = 'gray', mfc='lightgray', ms=6, markeredgewidth=1,ls='none',elinewidth=1)
ax0.set_xlabel(r'$V/V_i$',fontsize = 16)
ax0.set_ylabel(r'$\xi$',fontsize = 16)
# ax0.set_xlim(xmin=0)
# ax0.set_ylim(ymin=-10)
ax0.tick_params(direction='in',right='on',top='on')

xi_df = pd.read_csv(FeNi_filename)
ax1.errorbar(xi_df['V/Vi'],xi_df['xi'],yerr=xi_df['dxi'],marker = 'o',
	color = 'darkorange', mfc='#ffc681', ms=6, markeredgewidth=1,ls='none',elinewidth=1)
ax1.set_xlabel(r'$V/V_i$',fontsize = 16)
# ax1.set_ylabel(r'$\xi$',fontsize = 16)
# ax1.set_xlim(xmin=0)
# ax1.set_ylim(ymin=-10)
ax1.tick_params(direction='in',right='on',top='on')

xi_df = pd.read_csv(FeNiSi_filename)
ax2.errorbar(xi_df['V/Vi'],xi_df['xi'],yerr=xi_df['dxi'],marker = 'o',
	color = 'deepskyblue', mfc='#86e1ff', ms=6, markeredgewidth=1,ls='none',elinewidth=1)
ax2.set_xlabel(r'$V/V_i$',fontsize = 16)
# ax2.set_ylabel(r'$\xi$',fontsize = 16)
# ax2.set_xlim(xmin=0)
# ax2.set_ylim(ymin=-10)
ax2.tick_params(direction='in',right='on',top='on')

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/scalingparam.pdf', format='pdf',
	bbox_inches='tight')
plt.close()