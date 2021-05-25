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


# Input data
############

AICc_path = dict() # Key is tuple of composition and pressure
vD = dict()
dvD = dict()

key = ('FeNi',41)
# AICc_path[key] = '../091_vD_FeNi_PowerLaw/2012Apr_DAC11_P5/fit_search_results_Constrained_Power.csv'
AICc_path[key] = 'fit_search_results_Constrained_Power.csv'
vD[key] = 4.243
dvD[key] = 0.012


# Load data
###########

AICc_dict = dict()

for key in AICc_path.keys():
	AICc_dict[key] = pd.read_csv(AICc_path[key])


# Plot AICc vs vD, per Asimow's comment
######################################

def make_subplot(study,P,ax,xmin,xmax):
	key = (study,P)
	AICc_df = AICc_dict[key]
	ymin = min(AICc_df['AICc']) - 0.05*max(AICc_df['AICc'])
	ymax = 1.05*max(AICc_df['AICc'])
	ax.plot(AICc_df['vD']/1000,AICc_df['AICc'],'.',markersize=1)
	ax.plot(vD[key]*np.ones(2), [ymin,ymax], color='red',
		linewidth=2)
	ax.fill_between([vD[key]-dvD[key],vD[key]+dvD[key]], ymin*np.ones(2), ymax*np.ones(2),
		facecolor='red', alpha=0.35)

	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])
	ax.tick_params(direction='in',left='off')
	ax.set_xlabel(r'$v_D$ (km/s)',fontsize=18)
	ax.set_ylabel(r'AICc',fontsize=18)
	ax.yaxis.set_ticklabels([])

fig, (ax1) = plt.subplots(nrows = 1, ncols=1, figsize=(6,6))

study = 'FeNi'
P = 41 # GPa
make_subplot(study,P,ax1,3,5.5)

plt.tight_layout()
fig.savefig('AICc_vD_FeNicompare.pdf', format='pdf', transparent=True)
plt.close()


# Plot chi2 vs vD
######################################

def make_subplot(study,P,ax,xmin,xmax):
	key = (study,P)
	AICc_df = AICc_dict[key]
	ymin = min(AICc_df['chi2']) - 0.05*max(AICc_df['chi2'])
	ymax = 1.05*max(AICc_df['chi2'])
	ax.plot(AICc_df['vD']/1000,AICc_df['chi2'],'.',markersize=1)
	ax.plot(vD[key]*np.ones(2), [ymin,ymax], color='red',
		linewidth=2)
	ax.fill_between([vD[key]-dvD[key],vD[key]+dvD[key]], ymin*np.ones(2), ymax*np.ones(2),
		facecolor='red', alpha=0.35)

	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])
	ax.tick_params(direction='in',left='off')
	ax.set_xlabel(r'$v_D$ (km/s)',fontsize=18)
	ax.set_ylabel(r'$\chi^2$',fontsize=18)
	ax.yaxis.set_ticklabels([])

fig, (ax1) = plt.subplots(nrows = 1, ncols=1, figsize=(6,6))

study = 'FeNi'
P = 41 # GPa
make_subplot(study,P,ax1,3,5.5)

plt.tight_layout()
fig.savefig('chi2_vD_FeNicompare.pdf', format='pdf', transparent=True)
plt.close()
