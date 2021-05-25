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

phoxpath = dict()
EOSpath = dict()

phoxpath['Fe'] = '../050_phox_Fe_man/Results/phox_valsFromPDOS.csv'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis:
EOSpath['Fe'] = '../081_vD_Fe_PowerLaw/Results/input_values.csv'

phoxpath['FeNi'] = '../060_phox_FeNi_man/Results/phox_valsFromPDOS.csv'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis:
EOSpath['FeNi'] = '../091_vD_FeNi_PowerLaw/Results/input_values.csv'

phoxpath['FeNiSi'] = '../070_phox_FeNiSi_man/Results/phox_valsFromPDOS.csv'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis:
EOSpath['FeNiSi'] = '../101_vD_FeNiSi_PowerLaw/Results/input_values.csv'


# Functions
###########

def fitFvibLinear(V,a,b):
	Fvib = a*V + b
	return Fvib

def calcPvibLinear(V,a,b):
	Pvib = -a*np.ones(len(V))     # in mev/(atom A^3)
	Pvib = Pvib*2*1.6022*10**(-1)  # Unit conversion to GPa
	return Pvib

def calcdPvibLinear(V,cov):
	stddevs = np.sqrt(np.diag(cov))
	da = stddevs[0]
	dPvib = da*np.ones(len(V))
	return dPvib

def fitFvibQuad(V,a,b,c):
	Fvib = a*V**2 + b*V + c
	return Fvib

def calcPvibQuad(V,a,b,c):
	Pvib = -(2*a*V + b)             # in mev/(atom A^3)
	Pvib = Pvib*2*1.6022*10**(-1)  # Unit conversion to GPa
	return Pvib

def calcdPvibQuad(V,cov):
	# We only need the covariance between a and b
	cov = np.matrix(cov[:2,:2])
	# Calculate the Jacobian to propagate errors
	J = np.matrix([-2*V, -1*np.ones(len(V))]).T
	# Calculate the covariance matrix for Pvib
	Pvib_cov = J@(cov@J.T)
	# Reduce down to uncertainty for each Pvib
	dPvib = np.sqrt(np.diag(Pvib_cov))
	dPvib = dPvib*2*1.6022*10**(-1)  # Unit conversion to GPa
	return dPvib


# Import data
#############

input_dict = dict()
for study in ['Fe','FeNi','FeNiSi']:
	phox_df = pd.read_csv(phoxpath[study], engine='python')
	EOS_df  = pd.read_csv(EOSpath[study], engine='python')

	input_df = phox_df.merge(EOS_df, on=('Folder','Index'))
	# Correct for Pandas tendency to keep trailing rounding precision vals
	input_df = input_df.round(5)
	input_df.to_csv('Results/'+study+'_input_values.csv',index=False)
	input_dict[study] = input_df

# Fit Fvib with linear and polynomial function accounting for uncertainties
###########################################################################

lin_fit = dict()
lin_cov = dict()
quad_fit = dict()
quad_cov = dict()

for study in ['Fe','FeNi','FeNiSi']:
	print(study+':')
	input_df = input_dict[study]

	# Only fit hcp data
	input_df = input_df[input_df['Phase']=='hcp']

	lin_fit[study], lin_cov[study] = curve_fit(fitFvibLinear,input_df['V'],
		input_df['Fvib'], sigma=input_df['dFvib'])
	# lin_fit[study], lin_cov[study] = curve_fit(fitFvibLinear,input_df['V'],
	# 	input_df['Fvib'])
	# print(fit)
	# print(cov)
	quad_fit[study], quad_cov[study] = curve_fit(fitFvibQuad,input_df['V'],
		input_df['Fvib'], sigma=input_df['dFvib'])
	# quad_fit[study], quad_cov[study] = curve_fit(fitFvibQuad,input_df['V'],
	# 	input_df['Fvib'])
	# print(fit)
	# print(cov)


# Plot Fvib
###########

V_array = np.arange(15,21,0.1)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

input_df = input_dict['Fe']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],bcc_df['Fvib'],yerr=bcc_df['dFvib'],marker = 'o',
	color = 'gray', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['Fvib'],yerr=hcp_df['dFvib'],marker = 'o',
	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
# ax0.plot(V_array,fitFvibLinear(V_array,*lin_fit['Fe']),
# 	'--',color='gray',linewidth=1)
ax0.plot(V_array,fitFvibQuad(V_array,*quad_fit['Fe']),
	'-',color='gray',linewidth=1)

input_df = input_dict['FeNi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],bcc_df['Fvib'],yerr=bcc_df['dFvib'],marker = 'o',
	color = 'darkorange', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['Fvib'],yerr=hcp_df['dFvib'],marker = 'o',
	color = 'darkorange', mfc='#ffc681', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
# ax0.plot(V_array,fitFvibLinear(V_array,*lin_fit['FeNi']),
# 	'--',color='darkorange',linewidth=1)
ax0.plot(V_array,fitFvibQuad(V_array,*quad_fit['FeNi']),
	'-',color='darkorange',linewidth=1)

input_df = input_dict['FeNiSi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],bcc_df['Fvib'],yerr=bcc_df['dFvib'],marker = 'o',
	color = 'deepskyblue', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['Fvib'],yerr=hcp_df['dFvib'],marker = 'o',
	color = 'deepskyblue', mfc='#86e1ff', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
# ax0.plot(V_array,fitFvibLinear(V_array,*lin_fit['FeNiSi']),
# 	'--',color='deepskyblue', linewidth=1)
ax0.plot(V_array,fitFvibQuad(V_array,*quad_fit['FeNiSi']),
	'-',color='deepskyblue', linewidth=1)

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$F_{vib}$ (meV/atom)',fontsize = 16)
# ax0.set_xlim([15,21])
# ax0.set_ylim([0,75])
ax0.tick_params(direction='in',right='on',top='on')

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Fvib.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot Pvib
###########

V_array = np.arange(15,21,0.1)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

input_df = input_dict['Fe']
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.plot(V_array,calcPvibQuad(V_array,*quad_fit['Fe']),
	'-',color='gray',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*quad_fit['Fe']),
	yerr=calcdPvibQuad(hcp_df['V'],quad_cov['Fe']),marker = 'o',
	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)

input_df = input_dict['FeNi']
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.plot(V_array,calcPvibQuad(V_array,*quad_fit['FeNi']),
	'-',color='darkorange',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*quad_fit['FeNi']),
	yerr=calcdPvibQuad(hcp_df['V'],quad_cov['FeNi']),marker = 'o',
	color = 'darkorange', mfc='#ffc681', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)

input_df = input_dict['FeNiSi']
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.plot(V_array,calcPvibQuad(V_array,*quad_fit['FeNiSi']),
	'-',color='deepskyblue',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*quad_fit['FeNiSi']),
	yerr=calcdPvibQuad(hcp_df['V'],quad_cov['FeNiSi']),marker = 'o',
	color = 'deepskyblue', mfc='#86e1ff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$P^h_{vib}$ (GPa)',fontsize = 16)
# ax0.set_xlim([15,23])
ax0.set_ylim([0,4])
ax0.tick_params(direction='in',right='on',top='on')

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Pvib.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# # Plot Pvib
# ###########

# V_array = np.arange(15,21,0.1)

# fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

# input_df = input_dict['Fe']
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibLinear(V_array,*lin_fit['Fe']),
# 	'-',color='gray',linewidth=1)
# ax0.errorbar(hcp_df['V'],calcPvibLinear(hcp_df['V'],*lin_fit['Fe']),
# 	yerr=calcdPvibLinear(hcp_df['V'],lin_cov['Fe']),marker = 'o',
# 	color = 'gray', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
# 	capsize = 1)

# input_df = input_dict['FeNi']
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibLinear(V_array,*lin_fit['FeNi']),
# 	'-',color='darkorange',linewidth=1)
# ax0.errorbar(hcp_df['V'],calcPvibLinear(hcp_df['V'],*lin_fit['FeNi']),
# 	yerr=calcdPvibLinear(hcp_df['V'],lin_cov['FeNi']),marker = 'o',
# 	color = 'darkorange', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
# 	capsize = 1)

# input_df = input_dict['FeNiSi']
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibLinear(V_array,*lin_fit['FeNiSi']),
# 	'-',color='deepskyblue',linewidth=1)
# ax0.errorbar(hcp_df['V'],calcPvibLinear(hcp_df['V'],*lin_fit['FeNiSi']),
# 	yerr=calcdPvibLinear(hcp_df['V'],lin_cov['FeNiSi']),marker = 'o',
# 	color = 'deepskyblue', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
# 	capsize = 1)

# ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
# ax0.set_ylabel(r'$P^h_{vib}$ (GPa)',fontsize = 16)
# # ax0.set_xlim([15,23])
# ax0.set_ylim([0,10])
# ax0.tick_params(direction='in',right='on',top='on')

# # Fine-tune figure; make subplots close to each other and hide x ticks for
# # all but bottom plot.
# fig.subplots_adjust(wspace=0)

# fig.savefig('Plots/Pvib_lin_withdFvib.pdf', format='pdf',
# 	bbox_inches='tight')
# plt.close()



# # Fit Fvib with linear and polynomial function w/o accounting for uncertainties
# ###############################################################################

# lin_fit = dict()
# lin_cov = dict()
# quad_fit = dict()
# quad_cov = dict()

# for study in ['Fe','FeNi','FeNiSi']:
# 	print(study+':')
# 	input_df = input_dict[study]

# 	# Only fit hcp data
# 	input_df = input_df[input_df['Phase']=='hcp']

# 	lin_fit[study], lin_cov[study] = curve_fit(fitFvibLinear,input_df['V'],
# 		input_df['Fvib'])
# 	quad_fit[study], quad_cov[study] = curve_fit(fitFvibQuad,input_df['V'],
# 		input_df['Fvib'])


# # Plot Fvib
# ###########

# V_array = np.arange(15,21,0.1)

# fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

# input_df = input_dict['Fe']
# bcc_df = input_df[input_df['Phase']=='bcc']
# ax0.errorbar(bcc_df['V'],bcc_df['Fvib'],yerr=bcc_df['dFvib'],marker = 'o',
# 	color = 'gray', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.errorbar(hcp_df['V'],hcp_df['Fvib'],yerr=hcp_df['dFvib'],marker = 'o',
# 	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
# ax0.plot(V_array,fitFvibLinear(V_array,*lin_fit['Fe']),
# 	'--',color='gray',linewidth=1)
# ax0.plot(V_array,fitFvibQuad(V_array,*quad_fit['Fe']),
# 	'-',color='gray',linewidth=1)

# input_df = input_dict['FeNi']
# bcc_df = input_df[input_df['Phase']=='bcc']
# ax0.errorbar(bcc_df['V'],bcc_df['Fvib'],yerr=bcc_df['dFvib'],marker = 'o',
# 	color = 'darkorange', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.errorbar(hcp_df['V'],hcp_df['Fvib'],yerr=hcp_df['dFvib'],marker = 'o',
# 	color = 'darkorange', mfc='#ffc681', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
# ax0.plot(V_array,fitFvibLinear(V_array,*lin_fit['FeNi']),
# 	'--',color='darkorange',linewidth=1)
# ax0.plot(V_array,fitFvibQuad(V_array,*quad_fit['FeNi']),
# 	'-',color='darkorange',linewidth=1)

# input_df = input_dict['FeNiSi']
# bcc_df = input_df[input_df['Phase']=='bcc']
# ax0.errorbar(bcc_df['V'],bcc_df['Fvib'],yerr=bcc_df['dFvib'],marker = 'o',
# 	color = 'deepskyblue', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.errorbar(hcp_df['V'],hcp_df['Fvib'],yerr=hcp_df['dFvib'],marker = 'o',
# 	color = 'deepskyblue', mfc='#86e1ff', ms=7, markeredgewidth=1,ls='none',elinewidth=1)
# ax0.plot(V_array,fitFvibLinear(V_array,*lin_fit['FeNiSi']),
# 	'--',color='deepskyblue', linewidth=1)
# ax0.plot(V_array,fitFvibQuad(V_array,*quad_fit['FeNiSi']),
# 	'-',color='deepskyblue', linewidth=1)

# ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
# ax0.set_ylabel(r'$F_{vib}$ (meV/atom)',fontsize = 16)
# ax0.set_xlim([15,21])
# # ax0.set_ylim([0,75])
# ax0.tick_params(direction='in',right='on',top='on')

# # Fine-tune figure; make subplots close to each other and hide x ticks for
# # all but bottom plot.
# fig.subplots_adjust(wspace=0)

# fig.savefig('Plots/Fvib_nodFvib.pdf', format='pdf',
# 	bbox_inches='tight')
# plt.close()


# # Plot Pvib
# ###########

# V_array = np.arange(15,21,0.1)

# fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

# input_df = input_dict['Fe']
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*quad_fit['Fe']),
# 	'-',color='gray',linewidth=1)
# ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*quad_fit['Fe']),
# 	yerr=calcdPvibQuad(hcp_df['V'],quad_cov['Fe']),marker = 'o',
# 	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
# 	capsize = 1)

# input_df = input_dict['FeNi']
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*quad_fit['FeNi']),
# 	'-',color='darkorange',linewidth=1)
# ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*quad_fit['FeNi']),
# 	yerr=calcdPvibQuad(hcp_df['V'],quad_cov['FeNi']),marker = 'o',
# 	color = 'darkorange', mfc='#ffc681', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
# 	capsize = 1)

# input_df = input_dict['FeNiSi']
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*quad_fit['FeNiSi']),
# 	'-',color='deepskyblue',linewidth=1)
# ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*quad_fit['FeNiSi']),
# 	yerr=calcdPvibQuad(hcp_df['V'],quad_cov['FeNiSi']),marker = 'o',
# 	color = 'deepskyblue', mfc='#86e1ff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
# 	capsize = 1)

# ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
# ax0.set_ylabel(r'$P^h_{vib}$ (GPa)',fontsize = 16)
# # ax0.set_xlim([15,23])
# ax0.set_ylim([0,10])
# ax0.tick_params(direction='in',right='on',top='on')

# # Fine-tune figure; make subplots close to each other and hide x ticks for
# # all but bottom plot.
# fig.subplots_adjust(wspace=0)

# fig.savefig('Plots/Pvib_quad_nodFvib.pdf', format='pdf',
# 	bbox_inches='tight')
# plt.close()


# # Plot Pvib
# ###########

# V_array = np.arange(15,21,0.1)

# fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

# input_df = input_dict['Fe']
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibLinear(V_array,*lin_fit['Fe']),
# 	'-',color='gray',linewidth=1)
# ax0.errorbar(hcp_df['V'],calcPvibLinear(hcp_df['V'],*lin_fit['Fe']),
# 	yerr=calcdPvibLinear(hcp_df['V'],lin_cov['Fe']),marker = 'o',
# 	color = 'gray', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
# 	capsize = 1)

# input_df = input_dict['FeNi']
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibLinear(V_array,*lin_fit['FeNi']),
# 	'-',color='darkorange',linewidth=1)
# ax0.errorbar(hcp_df['V'],calcPvibLinear(hcp_df['V'],*lin_fit['FeNi']),
# 	yerr=calcdPvibLinear(hcp_df['V'],lin_cov['FeNi']),marker = 'o',
# 	color = 'darkorange', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
# 	capsize = 1)

# input_df = input_dict['FeNiSi']
# hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibLinear(V_array,*lin_fit['FeNiSi']),
# 	'-',color='deepskyblue',linewidth=1)
# ax0.errorbar(hcp_df['V'],calcPvibLinear(hcp_df['V'],*lin_fit['FeNiSi']),
# 	yerr=calcdPvibLinear(hcp_df['V'],lin_cov['FeNiSi']),marker = 'o',
# 	color = 'deepskyblue', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
# 	capsize = 1)

# ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
# ax0.set_ylabel(r'$P^h_{vib}$ (GPa)',fontsize = 16)
# # ax0.set_xlim([15,23])
# ax0.set_ylim([0,10])
# ax0.tick_params(direction='in',right='on',top='on')

# # Fine-tune figure; make subplots close to each other and hide x ticks for
# # all but bottom plot.
# fig.subplots_adjust(wspace=0)

# fig.savefig('Plots/Pvib_lin_nodFvib.pdf', format='pdf',
# 	bbox_inches='tight')
# plt.close()