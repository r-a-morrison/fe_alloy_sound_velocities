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


# Define input values
#####################

gamma_est_dict = dict()

gamma_est_dict['Fe'] = pd.read_csv('Fe/Results/GruneisenResults_method2.csv')
gamma_est_dict['FeNi'] = pd.read_csv('FeNi/Results/GruneisenResults_method2.csv')
gamma_est_dict['FeNiSi'] = pd.read_csv('FeNiSi/Results/GruneisenResults_method2.csv')

V0_dict = dict()
dV0_dict = dict()

# Fe V0 from Dewaele et al. 2006
V0_dict['Fe'] = 22.428
dV0_dict['Fe'] = 0.098

# Fe0.91Ni0.09 V0 from our EOS study
V0_dict['FeNi'] = 22.505
dV0_dict['FeNi'] = 0.042

# Fe0.8Ni0.1Si0.1 V0 from our EOS study
V0_dict['FeNiSi'] = 22.952
dV0_dict['FeNiSi'] = 0.072

V0 = dict()
dV0 = dict()
gamma0 = dict()
dgamma0 = dict()
q = dict()
dq = dict()

# Define input parameters for Gruneisen parameter calc: V0, gamma0, q
# Fe
fit = 'Fe vib q=0.8'
V0[fit] = V0_dict['Fe']
dV0[fit] = dV0_dict['Fe']
gamma0[fit] = 1.94
dgamma0[fit] = 0.01
q[fit] = 0.8

fit = 'Fe vib q=1.0'
V0[fit] = V0_dict['Fe']
dV0[fit] = dV0_dict['Fe']
gamma0[fit] = 2.04
dgamma0[fit] = 0.01
q[fit] = 1.0

fit = 'Fe vib q=1.2'
V0[fit] = V0_dict['Fe']
dV0[fit] = dV0_dict['Fe']
gamma0[fit] = 2.14
dgamma0[fit] = 0.01
q[fit] = 1.2

# fit = 'Fe Debye q=0.8'
# V0[fit] = V0_dict['Fe']
# dV0[fit] = dV0_dict['Fe']
# gamma0[fit] = 1.68
# dgamma0[fit] = 0.02
# q[fit] = 0.8

# fit = 'Fe Debye q=1.0'
# V0[fit] = V0_dict['Fe']
# dV0[fit] = dV0_dict['Fe']
# gamma0[fit] = 1.74
# dgamma0[fit] = 0.03
# q[fit] = 1.0

# fit = 'Fe Debye q=1.2'
# V0[fit] = V0_dict['Fe']
# dV0[fit] = dV0_dict['Fe']
# gamma0[fit] = 1.79
# dgamma0[fit] = 0.03
# q[fit] = 1.2

# Fe0.91Ni0.09
fit = 'FeNi vib q=0.8'
V0[fit] = V0_dict['FeNi']
dV0[fit] = dV0_dict['FeNi']
gamma0[fit] = 1.99
dgamma0[fit] = 0.02
q[fit] = 0.8

fit = 'FeNi vib q=1.0'
V0[fit] = V0_dict['FeNi']
dV0[fit] = dV0_dict['FeNi']
gamma0[fit] = 2.07
dgamma0[fit] = 0.02
q[fit] = 1.0

fit = 'FeNi vib q=1.2'
V0[fit] = V0_dict['FeNi']
dV0[fit] = dV0_dict['FeNi']
gamma0[fit] = 2.15
dgamma0[fit] = 0.02
q[fit] = 1.2

# fit = 'FeNi Debye q=0.8'
# V0[fit] = V0_dict['FeNi']
# dV0[fit] = dV0_dict['FeNi']
# gamma0[fit] = 1.65
# dgamma0[fit] = 0.02
# q[fit] = 0.8

# fit = 'FeNi Debye q=1.0'
# V0[fit] = V0_dict['FeNi']
# dV0[fit] = dV0_dict['FeNi']
# gamma0[fit] = 1.69
# dgamma0[fit] = 0.02
# q[fit] = 1.0

# fit = 'FeNi Debye q=1.2'
# V0[fit] = V0_dict['FeNi']
# dV0[fit] = dV0_dict['FeNi']
# gamma0[fit] = 1.72
# dgamma0[fit] = 0.02
# q[fit] = 1.2

# Fe0.8Ni0.1Si0.1
fit = 'FeNiSi vib q=0.8'
V0[fit] = V0_dict['FeNiSi']
dV0[fit] = dV0_dict['FeNiSi']
gamma0[fit] =  1.95
dgamma0[fit] =  0.04
q[fit] = 0.8

fit = 'FeNiSi vib q=1.0'
V0[fit] = V0_dict['FeNiSi']
dV0[fit] = dV0_dict['FeNiSi']
gamma0[fit] =  2.03
dgamma0[fit] = 0.05
q[fit] = 1.0

fit = 'FeNiSi vib q=1.2'
V0[fit] = V0_dict['FeNiSi']
dV0[fit] = dV0_dict['FeNiSi']
gamma0[fit] = 2.12
dgamma0[fit] = 0.05
q[fit] = 1.2

# fit = 'FeNiSi Debye q=0.8'
# V0[fit] = V0_dict['FeNiSi']
# dV0[fit] = dV0_dict['FeNiSi']
# gamma0[fit] = 1.78
# dgamma0[fit] = 0.06
# q[fit] = 0.8

# fit = 'FeNiSi Debye q=1.0'
# V0[fit] = V0_dict['FeNiSi']
# dV0[fit] = dV0_dict['FeNiSi']
# gamma0[fit] = 1.82
# dgamma0[fit] = 0.06
# q[fit] = 1.0

# fit = 'FeNiSi Debye q=1.2'
# V0[fit] = V0_dict['FeNiSi']
# dV0[fit] = dV0_dict['FeNiSi']
# gamma0[fit] = 1.85
# dgamma0[fit] = 0.06
# q[fit] = 1.2

fit = 'Merkel2000'
V0[fit] = V0_dict['Fe']
gamma0[fit] = 1.68
dgamma0[fit] = 0.2
q[fit] = 0.7
dq[fit] = 0.5

fit = 'Dubrovinksy2000'
V0[fit] = V0_dict['Fe']
gamma0[fit] = 1.68
dgamma0[fit] = 0.06
q[fit] = 0.69
dq[fit] = 0.10

def Dewaele_func(V):
	V0 = V0_dict['Fe']
	gamma0 = 1.875
	gammainf = 1.305
	x = V/V0
	beta = gamma0/(gamma0-gammainf)
	gamma = gammainf + (gamma0-gammainf)*x**beta
	return gamma

Sha_df = pd.read_csv('../006_Literature/TheoryWithSameUnits/Sha2010_hcpFe.csv')

fit = 'Fei2016'
# Their reported rho0 below is equivalent to our V0
# rho0 = 8.2695 # g/cc
# M = 55.845    # g/mol
# V0_ccpermol = M/rho0  # cm^3/mol
# V0[fit] = V0_ccpermol*(2*10**24)/constants.N_A  # A^3
V0[fit] = V0_dict['Fe']
gamma0[fit] = 1.74
q[fit] = 0.78

fit = 'Murphy2011a'
V0[fit] = V0_dict['Fe']
gamma0[fit] = 1.98
dgamma0[fit] = 0.02
q[fit] = 1


# Functions
###########

def calcGruneisen(V,V0,gamma0,q):
	gamma = gamma0*(V/V0)**q
	return gamma


# Plot comparing Fe, Fe-Ni, Fe-Ni-Si Gruneisen parameters
#########################################################

V_array = np.arange(15,22,0.01)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(8,6))

# Fe
color = 'gray'
# gamma_est_df = gamma_est_dict['Fe']
# ax0.errorbar(gamma_est_df['Vi'],gamma_est_df['gamma_est'],
# 	# xerr=gamma_est_df['dV'],
# 	yerr = gamma_est_df['dgamma_est'], marker='o',
# 	color=color,ls='None',label=r'$\gamma_{vib}$')

fit = 'Fe vib q=1.0'
ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '-',
	color = color, linewidth = 1.0)
ax0.fill_between(V_array,
	calcGruneisen(V_array,V0[fit],gamma0[fit]-dgamma0[fit],q[fit]),
	calcGruneisen(V_array,V0[fit],gamma0[fit]+dgamma0[fit],q[fit]),
	facecolor=color, alpha=0.3)

# Fe0.91Ni0.09
color = 'darkorange'
# gamma_est_df = gamma_est_dict['FeNi']
# ax0.errorbar(gamma_est_df['Vi'],gamma_est_df['gamma_est'],
# 	# xerr=gamma_est_df['dV'],
# 	yerr = gamma_est_df['dgamma_est'], marker='o',
# 	color=color,ls='None',label=r'$\gamma_{vib}$')

fit = 'FeNi vib q=1.0'
ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '-',
	color = color, linewidth = 1.0)
ax0.fill_between(V_array,
	calcGruneisen(V_array,V0[fit],gamma0[fit]-dgamma0[fit],q[fit]),
	calcGruneisen(V_array,V0[fit],gamma0[fit]+dgamma0[fit],q[fit]),
	facecolor=color, alpha=0.3)

# Fe0.8Ni0.1Si0.1
color = 'deepskyblue'
# gamma_est_df = gamma_est_dict['FeNiSi']
# ax0.errorbar(gamma_est_df['Vi'],gamma_est_df['gamma_est'],
# 	# xerr=gamma_est_df['dV'],
# 	yerr = gamma_est_df['dgamma_est'], marker='o',
# 	color=color,ls='None',label=r'$\gamma_{vib}$')

fit = 'FeNiSi vib q=1.0'
ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '-',
	color = color, linewidth = 1.0)
ax0.fill_between(V_array,
	calcGruneisen(V_array,V0[fit],gamma0[fit]-dgamma0[fit],q[fit]),
	calcGruneisen(V_array,V0[fit],gamma0[fit]+dgamma0[fit],q[fit]),
	facecolor=color, alpha=0.3)
ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$\gamma_{vib}$',fontsize = 16)
ax0.set_xlim([15,22])
# ax0.set_ylim(ymin=-10)
ax0.tick_params(direction='in',right='on',top='on')

# fig.subplots_adjust(wspace=0)

ax0.set_xlim([15,21])
ax0.set_ylim([1.25,2.0])
fig.savefig('Plots/GruneisenParam_est_2.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot comparing Fe Gruneisen parameters
########################################

results_df = pd.read_csv('Fe/Results/GruneisenResults.csv')
V_array = np.arange(15,22,0.01)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(8,6))

# Fe
color = 'gray'
fit = 'Fe vib q=1.0'
ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '-',
	color = color, linewidth = 1.0, label='Fe (re-analyzed from Murphy et al. 2011a)')
ax0.fill_between(V_array,
	calcGruneisen(V_array,V0[fit],gamma0[fit]-dgamma0[fit],q[fit]),
	calcGruneisen(V_array,V0[fit],gamma0[fit]+dgamma0[fit],q[fit]),
	facecolor=color, alpha=0.3)

# Fe0.91Ni0.09
color = 'darkorange'
fit = 'FeNi vib q=1.0'
ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '-',
	color = color, linewidth = 1.0, label='Fe$_{0.91}$Ni$_{0.09}$ (this study)')
ax0.fill_between(V_array,
	calcGruneisen(V_array,V0[fit],gamma0[fit]-dgamma0[fit],q[fit]),
	calcGruneisen(V_array,V0[fit],gamma0[fit]+dgamma0[fit],q[fit]),
	facecolor=color, alpha=0.3)

# Fe0.8Ni0.1Si0.1
color = 'deepskyblue'
gamma_est_df = gamma_est_dict['FeNiSi']
fit = 'FeNiSi vib q=1.0'
ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '-',
	color = color, linewidth = 1.0, label='Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$ (this study)')
ax0.fill_between(V_array,
	calcGruneisen(V_array,V0[fit],gamma0[fit]-dgamma0[fit],q[fit]),
	calcGruneisen(V_array,V0[fit],gamma0[fit]+dgamma0[fit],q[fit]),
	facecolor=color, alpha=0.3)

fit = 'Fei2016'
ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '-.',
	color = 'limegreen', linewidth = 1.25, label = 'Fe (Fei et al. 2016)')
# Sha and Cohen 2010
ax0.plot(Sha_df['V'],Sha_df['gamma'], 'o', color = 'green', markerfacecolor='None',
	linewidth = 1.25, markersize = 5, label = 'Fe (Sha and Cohen 2010)')
# Dewaele et al. 2006
ax0.plot(V_array,Dewaele_func(V_array), '--', color = 'purple', linewidth = 1.25,
	label = 'Fe (Dewaele et al. 2006)')
fit = 'Merkel2000'
ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), ':',
	color = 'blue', linewidth = 1.25, label = 'Fe (Merkel et al. 2000)')
fit = 'Dubrovinksy2000'
ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '--',
	color = 'red', linewidth = 1.25, label = 'Fe (Dubrovinksy et al. 2000)')
# fit = 'Murphy2011a'
# ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '--',
# 	color = 'black', linewidth = 1.25, label = 'Murphy et al. 2011a')

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$\gamma_{vib}$',fontsize = 16)
ax0.set_xlim([15,22])
ax0.set_ylim([1.2,2.1])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=2,fontsize=11)

# fig.subplots_adjust(wspace=0)

fig.savefig('Plots/GruneisenParam.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot comparing Fe Gruneisen parameters
########################################

results_df = pd.read_csv('Fe/Results/GruneisenResults.csv')
V_array = np.arange(15,22,0.01)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(8,6))

# Fe
color = 'gray'

fit = 'Fe vib q=1.0'
h0, = ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '-',
	color = color,linewidth = 1.0,label=r'This study (Equation 17)')
ax0.fill_between(V_array,
	calcGruneisen(V_array,V0[fit],gamma0[fit]-dgamma0[fit],q[fit]),
	calcGruneisen(V_array,V0[fit],gamma0[fit]+dgamma0[fit],q[fit]),
	facecolor=color, alpha=0.3)

# gamma_est_df = gamma_est_dict['Fe']
# h1, = ax0.plot(gamma_est_df['Vi'],gamma_est_df['gamma_est'],
# 	marker='o',color=color,ls='None',label=r'This study (equation 4.20)')
# ax0.errorbar(gamma_est_df['Vi'],gamma_est_df['gamma_est'],
# 	# xerr=gamma_est_df['dV'],
# 	yerr = gamma_est_df['dgamma_est'], marker='o',
# 	color=color,ls='None')

fit = 'Fei2016'
h2, = ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '-.',
	color = 'limegreen', linewidth = 1.25, label = 'Fei et al. 2016')
fit = 'Murphy2011a'
h3, = ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '--',
	color = 'black', linewidth = 1.25, label = 'Murphy et al. 2011a')

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$\gamma_{vib}$',fontsize = 16)
ax0.set_xlim([15,22])
# ax0.set_ylim(ymin=-10)
ax0.tick_params(direction='in',right='on',top='on')

hlabels = [h0,h2,h3]
ax0.legend(loc=2,fontsize=11,handles = hlabels)

# fig.subplots_adjust(wspace=0)

fig.savefig('Plots/GruneisenParam_Fe.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot comparing Fe Gruneisen parameters again for some dumb reason
########################################

results_df = pd.read_csv('Fe/Results/GruneisenResults.csv')
V_array = np.arange(15,22,0.01)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(8,6))

# Fe
color = 'black'

fit = 'Fe vib q=1.0'
h0, = ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '-',
	color = color,linewidth = 1.0,label=r'This study (Equation 17)')
ax0.fill_between(V_array,
	calcGruneisen(V_array,V0[fit],gamma0[fit]-dgamma0[fit],q[fit]),
	calcGruneisen(V_array,V0[fit],gamma0[fit]+dgamma0[fit],q[fit]),
	facecolor=color, alpha=0.3)

# gamma_est_df = gamma_est_dict['Fe']
# h1, = ax0.plot(gamma_est_df['Vi'],gamma_est_df['gamma_est'],
# 	marker='o',color=color,ls='None',label=r'This study (equation 4.20)')
# ax0.errorbar(gamma_est_df['Vi'],gamma_est_df['gamma_est'],
# 	# xerr=gamma_est_df['dV'],
# 	yerr = gamma_est_df['dgamma_est'], marker='o',
# 	color=color,ls='None')

# fit = 'Fei2016'
# h2, = ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '-.',
# 	color = 'limegreen', linewidth = 1.25, label = 'Fei et al. 2016')
fit = 'Murphy2011a'
h3, = ax0.plot(V_array,calcGruneisen(V_array,V0[fit],gamma0[fit],q[fit]), '--',
	color = 'red', linewidth = 1.25, label = 'Murphy et al. 2011a')

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$\gamma_{vib}$',fontsize = 16)
ax0.set_xlim([15,22])
# ax0.set_ylim(ymin=-10)
ax0.tick_params(direction='in',right='on',top='on')

hlabels = [h0,h3]
ax0.legend(loc=2,fontsize=11,handles = hlabels)

# fig.subplots_adjust(wspace=0)

fig.savefig('Plots/GruneisenParam_Fe2.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot comparing Pvib for Fe, Fe-Ni, Fe-Ni-Si
#############################################

V_array = np.arange(15,22,0.01)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(8,6))

# Fe
results_df = pd.read_csv('Fe/Results/GruneisenResults.csv')
color = 'gray'
ax0.errorbar(results_df['V'],results_df['Pvib'],
	xerr=results_df['dV'],yerr = results_df['dPvib'], marker='o',
	color=color,ls='None',label=r'Fe')

# Fe0.91Ni0.09
results_df = pd.read_csv('FeNi/Results/GruneisenResults.csv')
color = 'darkorange'
ax0.errorbar(results_df['V'],results_df['Pvib'],
	xerr=results_df['dV'],yerr = results_df['dPvib'], marker='o',
	color=color,ls='None',label=r'Fe$_{0.91}$Ni$_{0.09}$')

# Fe0.8Ni0.1Si0.1
color = 'deepskyblue'
results_df = pd.read_csv('FeNiSi/Results/GruneisenResults.csv')
ax0.errorbar(results_df['V'],results_df['Pvib'],
	xerr=results_df['dV'],yerr = results_df['dPvib'], marker='o',
	color=color,ls='None',label=r'Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$')

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$P_{vib}$',fontsize = 16)
# ax0.set_xlim([15,22])
# ax0.set_ylim(ymin=-10)
ax0.tick_params(direction='in',right='on',top='on')

# fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Pvib.pdf', format='pdf',
	bbox_inches='tight')
plt.close()