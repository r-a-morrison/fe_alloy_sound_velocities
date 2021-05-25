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
V0 = 22.952
dV0 = 0.072

xi_filename = 'Results/scalingparameters.csv'

# Values from previous analyses, e.g. volume and vD
input_filename = 'Results/input_values.csv'


# Functions
###########

# model for xi
def xiModel(V,Vi,V0,gamma0,q):
	return np.exp( gamma0*(Vi/V0)**q * (1/q) * ((V/Vi)**q-1) )

def misfit(V,Vi,V0,xi,gamma0,q):
	# chi2 = sum((xi-xiModel(V,Vi,V0,gamma0,q))**2/(xi_df['dxi'])**2)
	chi2 = sum((xi-xiModel(V,Vi,V0,gamma0,q))**2)
	reducedchi2 = chi2/(len(xi_df)-2)
	return reducedchi2

# Creates a model to pass into curve_fit
def make_xiModel(V0,q):
	# Model for xi
	def model(Vdata,gamma0):
		V, Vi = Vdata
		xi = np.exp( gamma0*(Vi/V0)**q * (1/q) * ((V/Vi)**q-1) )
		return xi
	return model


# Load scaling parameter results and other info
##########################################

xi_df = pd.read_csv(xi_filename)
# V = np.array(xi_df['V'])
# Vi = np.array(xi_df['Vi'])
# xi = np.array(xi_df['xi'])

input_df = pd.read_csv(input_filename)

# Gruneisen parameter grid search
#################################

# Try to fit with a grid search
# Define search range of gamma0 and q
gamma0_range = np.arange(1.6,2.2 ,0.01)
q_range = np.arange(0.1,1.3,0.01)

# Create two arrays
gamma0_array = np.zeros(len(gamma0_range)*len(q_range))
q_array = np.zeros(len(gamma0_range)*len(q_range))
i = 0
for gamma0 in gamma0_range:
	for q in q_range:
		gamma0_array[i] = gamma0
		q_array[i] = q
		i += 1

# Calculate chi^2 and store in a dataframe
chi2_array = [misfit(xi_df['V'],xi_df['Vi'],V0,xi_df['xi'],gamma0,q)
	for gamma0, q in zip(gamma0_array,q_array)]
fit_df = pd.DataFrame(columns=['gamma0','q','chi2'])
fit_df['gamma0'] = gamma0_array
fit_df['q'] = q_array
fit_df['chi2'] = chi2_array

pivotplot = fit_df.pivot('gamma0', 'q', 'chi2')
gamma_array = pivotplot.index.values
q_array = pivotplot.columns.values

# Plot chi2 as a function of gamma0 and q
fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
cf0 = ax0.pcolormesh(q_array,gamma_array,pivotplot,edgecolor='face')
cbar0 = plt.colorbar(cf0, ax = ax0)
cbar0.ax.set_ylabel(r'$\chi^2$',fontsize=16)
ax0.set_xlabel('$q$',fontsize=16)
ax0.set_ylabel('$\gamma_0$',fontsize=16)
fig.savefig('Results/FeNiSi_Gruneisen_misfit.pdf', format='pdf',
	bbox_inches='tight')
plt.close()

# Plot pdf as a function of gamma0 and q
fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
cf0 = ax0.pcolormesh(q_array,gamma_array,np.exp(-pivotplot),edgecolor='face')
cbar0 = plt.colorbar(cf0, ax = ax0)
cbar0.ax.set_ylabel(r'pdf',fontsize=16)
ax0.set_xlabel('$q$',fontsize=16)
ax0.set_ylabel('$\gamma_0$',fontsize=16)
fig.savefig('Results/FeNiSi_Gruneisen_pdf.pdf', format='pdf',
	bbox_inches='tight')
plt.close()

fit_df = fit_df.round({'gamma0':2,'q':2})

# If q = 0.8
fit_q08_df = fit_df[fit_df['q']==0.8]
gamma0_when_q08 = fit_q08_df.loc[fit_q08_df['chi2'].idxmin()]['gamma0']
print(gamma0_when_q08)

fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
ax0.plot(fit_q08_df['gamma0'],fit_q08_df['chi2'])
ax0.set_xlabel('$\gamma_0$',fontsize=16)
ax0.set_ylabel('$\chi^2$',fontsize=16)
# ax0.set_ylim([350,600])
ax0.set_title('q=0.8')
fig.savefig('Results/FeNiSi_misfit_q0.8.pdf', format='pdf',
	bbox_inches='tight')
plt.close()

# # Create normalized pdf and plot
# pdf = np.exp(-fit_q08_df['chi2'])
# # Normalize pdf
# stepsize = (max(fit_q08_df['gamma0'])-min(fit_q08_df['gamma0']))/len(fit_q08_df)
# normconst = 1/(np.sum(pdf)*stepsize)
# pdf = normconst*pdf
# fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
# ax0.plot(fit_q08_df['gamma0'],pdf)
# ax0.set_xlabel('$\gamma_0$',fontsize=16)
# ax0.set_ylabel('pdf',fontsize=16)
# # ax0.set_ylim([0,10**(-160)])
# ax0.set_title('q=0.8')
# fig.savefig('Results/FeNiSi_pdf_q0.8.pdf', format='pdf',
# 	bbox_inches='tight')
# plt.close()

# # Calculate std dev of gamma0
# mean = np.average(fit_q08_df['gamma0'], weights=pdf)
# var = np.average((fit_q08_df['gamma0']-gamma0_when_q08)**2, weights=pdf)
# dgamma0_when_q08 = np.sqrt(var)
# print(dgamma0_when_q08)

q=0.8
fit, cov = curve_fit(make_xiModel(V0,q),(xi_df['V'],xi_df['Vi']),xi_df['xi'],
	p0=[2.0],sigma=xi_df['dxi'])
print('q='+str(q)+':')
print(fit)
stddev = np.sqrt(np.diag(cov))
print(stddev)


# If q = 1.0
fit_q1_df = fit_df[fit_df['q']==1.0]
gamma0_when_q1 = fit_q1_df.loc[fit_q1_df['chi2'].idxmin()]['gamma0']
print(gamma0_when_q1)

fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
ax0.plot(fit_q1_df['gamma0'],fit_q1_df['chi2'])
ax0.set_xlabel('$\gamma_0$',fontsize=16)
ax0.set_ylabel('$\chi^2$',fontsize=16)
# ax0.set_ylim([350,600])
ax0.set_title('q=1.0')
fig.savefig('Results/FeNiSi_misfit_q1.0.pdf', format='pdf',
	bbox_inches='tight')
plt.close()

# # Create normalized pdf and plot
# pdf = np.exp(-fit_q1_df['chi2'])
# # Normalize pdf
# stepsize = (max(fit_q1_df['gamma0'])-min(fit_q1_df['gamma0']))/len(fit_q1_df)
# normconst = 1/(np.sum(pdf)*stepsize)
# pdf = normconst*pdf
# fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
# ax0.plot(fit_q1_df['gamma0'],pdf)
# ax0.set_xlabel('$\gamma_0$',fontsize=16)
# ax0.set_ylabel('pdf',fontsize=16)
# # ax0.set_ylim([0,10**(-160)])
# ax0.set_title('q=1.0')
# fig.savefig('Results/FeNiSi_pdf_q1.0.pdf', format='pdf',
# 	bbox_inches='tight')
# plt.close()

# # Calculate std dev of gamma0
# mean = np.average(fit_q1_df['gamma0'], weights=pdf)
# var = np.average((fit_q1_df['gamma0']-gamma0_when_q1)**2, weights=pdf)
# dgamma0_when_q1 = np.sqrt(var)
# print(dgamma0_when_q1)

q=1.0
fit, cov = curve_fit(make_xiModel(V0,q),(xi_df['V'],xi_df['Vi']),xi_df['xi'],
	p0=[2.0],sigma=xi_df['dxi'])
print('q='+str(q)+':')
print(fit)
stddev = np.sqrt(np.diag(cov))
print(stddev)



# If q = 1.2
fit_q12_df = fit_df[fit_df['q']==1.2]
gamma0_when_q12 = fit_q12_df.loc[fit_q12_df['chi2'].idxmin()]['gamma0']
print(gamma0_when_q12)

fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
ax0.plot(fit_q12_df['gamma0'],fit_q12_df['chi2'])
ax0.set_xlabel('$\gamma_0$',fontsize=16)
ax0.set_ylabel('$\chi^2$',fontsize=16)
# ax0.set_ylim([350,600])
ax0.set_title('q=1.2')
fig.savefig('Results/FeNiSi_misfit_q1.2.pdf', format='pdf',
	bbox_inches='tight')
plt.close()

# # Create normalized pdf and plot
# pdf = np.exp(-fit_q12_df['chi2'])
# # Normalize pdf
# stepsize = (max(fit_q12_df['gamma0'])-min(fit_q12_df['gamma0']))/len(fit_q12_df)
# normconst = 1/(np.sum(pdf)*stepsize)
# pdf = normconst*pdf
# fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(6,4.5))
# ax0.plot(fit_q12_df['gamma0'],pdf)
# ax0.set_xlabel('$\gamma_0$',fontsize=16)
# ax0.set_ylabel('pdf',fontsize=16)
# # ax0.set_ylim([0,10**(-160)])
# ax0.set_title('q=1.2')
# fig.savefig('Results/FeNiSi_pdf_q1.2.pdf', format='pdf',
# 	bbox_inches='tight')
# plt.close()

# # Calculate std dev of gamma0
# mean = np.average(fit_q12_df['gamma0'], weights=pdf)
# var = np.average((fit_q12_df['gamma0']-gamma0_when_q12)**2, weights=pdf)
# dgamma0_when_q12 = np.sqrt(var)
# print(dgamma0_when_q12)

q=1.2
fit, cov = curve_fit(make_xiModel(V0,q),(xi_df['V'],xi_df['Vi']),xi_df['xi'],
	p0=[2.0],sigma=xi_df['dxi'])
print('q='+str(q)+':')
print(fit)
stddev = np.sqrt(np.diag(cov))
print(stddev)
