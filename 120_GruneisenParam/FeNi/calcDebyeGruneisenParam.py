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
# Our EOS
V0 = 22.505
dV0 = 0.042

# Values from previous analyses, e.g. volume and vD
input_filename = 'Results/input_values.csv'


# Functions
###########

# model for xi
def vD_model(V,V0,vD0,gammaD0,q):
	return vD0 * (V/V0)**(1/3) * np.exp( -(gammaD0/q)*((V/V0)**q-1) )

def misfit(vD,dvD,V,V0,vD0,gammaD0,q):
	chi2 = sum((vD-vD_model(V,V0,vD0,gammaD0,q))**2/(dvD)**2)
	reducedchi2 = chi2/(len(vD)-3)
	return reducedchi2

# Creates a model to pass into curve_fit. Fix q
def make_vD_model(V,V0,vD0,gammaD0,q):
	# Model for xi
	def model(V,vD0,gammaD0):
		vD = vD0 * (V/V0)**(1/3) * np.exp( -(gammaD0/q)*((V/V0)**q-1) )
		return vD
	return model


# Load input info
#################

input_df = pd.read_csv(input_filename)


# Gruneisen parameter grid search
#################################

# Try to fit with a grid search
# Define search range of gamma0 and q
vD0_range = np.arange(4500,5500,20)
gammaD0_range = np.arange(1.0,2.5 ,0.05)
q_range = np.arange(0.1,2.0,0.1)

# Create three arrays
vD0_array = np.zeros(len(gammaD0_range)*len(q_range)*len(vD0_range))
gammaD0_array = np.zeros(len(gammaD0_range)*len(q_range)*len(vD0_range))
q_array = np.zeros(len(gammaD0_range)*len(q_range)*len(vD0_range))
i = 0
for vD0 in vD0_range:
	for gammaD0 in gammaD0_range:
		for q in q_range:
			vD0_array[i] = vD0
			gammaD0_array[i] = gammaD0
			q_array[i] = q
			i += 1

# Calculate chi^2 and store in a dataframe
chi2_array = [misfit(input_df['vD'],input_df['dvD'],input_df['V'],V0,vD0,gammaD0,q)
	for vD0,gammaD0,q in zip(vD0_array,gammaD0_array,q_array)]
fit_df = pd.DataFrame(columns=['vD0','gammaD0','q','chi2'])
fit_df['vD0'] = vD0_array
fit_df['gammaD0'] = gammaD0_array
fit_df['q'] = q_array
fit_df['chi2'] = chi2_array

fit_df = fit_df.round(6)

# Fix q = 0.8
q=0.8
fit_red_df = fit_df[fit_df['q']==q]
pivotplot = fit_red_df.pivot('vD0','gammaD0', 'chi2')
vD0_pivot = pivotplot.index.values
gammaD0_pivot = pivotplot.columns.values
# Plot chi2 as a function of gamma0 and q
fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
cf0 = ax0.pcolormesh(gammaD0_pivot,vD0_pivot,pivotplot,edgecolor='face',vmax=200)
cbar0 = plt.colorbar(cf0, ax = ax0)
cbar0.ax.set_ylabel(r'$\chi^2$',fontsize=16)
ax0.set_xlabel(r'$\gamma_{D,0}$',fontsize=16)
ax0.set_ylabel(r'$v_{D,0}$',fontsize=16)
ax0.set_title(r'$q=$'+str(q))
fig.savefig('Results/FeNi_DebyeGruneisen_misfit1.pdf', format='pdf',
	bbox_inches='tight')
plt.close()

fit, cov = curve_fit(make_vD_model(input_df['V'],V0,vD0,gammaD0,q),
	input_df['V'],input_df['vD'],p0=[4000,2],sigma=input_df['dvD'])
print('q='+str(q)+':')
print(fit)
stddev = np.sqrt(np.diag(cov))
print(stddev)
print(cov)
print(np.divide( cov, np.matmul(np.matrix(stddev).T,np.matrix(stddev))  ))

# Fix q = 1.0
q=1.0
fit_red_df = fit_df[fit_df['q']==q]
pivotplot = fit_red_df.pivot('vD0','gammaD0', 'chi2')
vD0_pivot = pivotplot.index.values
gammaD0_pivot = pivotplot.columns.values
# Plot chi2 as a function of gamma0 and q
fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
cf0 = ax0.pcolormesh(gammaD0_pivot,vD0_pivot,pivotplot,edgecolor='face',vmax=200)
cbar0 = plt.colorbar(cf0, ax = ax0)
cbar0.ax.set_ylabel(r'$\chi^2$',fontsize=16)
ax0.set_xlabel(r'$\gamma_{D,0}$',fontsize=16)
ax0.set_ylabel(r'$v_{D,0}$',fontsize=16)
ax0.set_title(r'$q=$'+str(q))
fig.savefig('Results/FeNi_DebyeGruneisen_misfit2.pdf', format='pdf',
	bbox_inches='tight')
plt.close()

fit, cov = curve_fit(make_vD_model(input_df['V'],V0,vD0,gammaD0,q),
	input_df['V'],input_df['vD'],p0=[4000,2],sigma=input_df['dvD'])
print('q='+str(q)+':')
print(fit)
stddev = np.sqrt(np.diag(cov))
print(stddev)
print(cov)
print(np.divide( cov, np.matmul(np.matrix(stddev).T,np.matrix(stddev))  ))

# Fix q = 1.2
q=1.2
fit_red_df = fit_df[fit_df['q']==q]
pivotplot = fit_red_df.pivot('vD0','gammaD0', 'chi2')
vD0_pivot = pivotplot.index.values
gammaD0_pivot = pivotplot.columns.values
# Plot chi2 as a function of gamma0 and q
fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
cf0 = ax0.pcolormesh(gammaD0_pivot,vD0_pivot,pivotplot,edgecolor='face',vmax=200)
cbar0 = plt.colorbar(cf0, ax = ax0)
cbar0.ax.set_ylabel(r'$\chi^2$',fontsize=16)
ax0.set_xlabel(r'$\gamma_{D,0}$',fontsize=16)
ax0.set_ylabel(r'$v_{D,0}$',fontsize=16)
ax0.set_title(r'$q=$'+str(q))
fig.savefig('Results/FeNi_DebyeGruneisen_misfit3.pdf', format='pdf',
	bbox_inches='tight')
plt.close()

fit, cov = curve_fit(make_vD_model(input_df['V'],V0,vD0,gammaD0,q),
	input_df['V'],input_df['vD'],p0=[4000,2],sigma=input_df['dvD'])
print('q='+str(q)+':')
print(fit)
stddev = np.sqrt(np.diag(cov))
print(stddev)
print(cov)
print(np.divide( cov, np.matmul(np.matrix(stddev).T,np.matrix(stddev))  ))

# Fix vD = 3700
fit_red_df = fit_df[fit_df['vD0']==5100]
pivotplot = fit_red_df.pivot('gammaD0', 'q', 'chi2')
gammaD0_pivot = pivotplot.index.values
q_pivot = pivotplot.columns.values
# Plot chi2 as a function of gamma0 and q
fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
cf0 = ax0.pcolormesh(q_pivot,gammaD0_pivot,pivotplot,edgecolor='face',vmax=200)
cbar0 = plt.colorbar(cf0, ax = ax0)
cbar0.ax.set_ylabel(r'$\chi^2$',fontsize=16)
ax0.set_xlabel(r'$q$',fontsize=16)
ax0.set_ylabel(r'$\gamma_{D,0}$',fontsize=16)
ax0.set_title(r'$v_{D,0}=5100$')
fig.savefig('Results/FeNi_DebyeGruneisen_misfit4.pdf', format='pdf',
	bbox_inches='tight')
plt.close()



print(time.time() - start_time)