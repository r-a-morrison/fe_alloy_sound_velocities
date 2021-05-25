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
refDatapath = dict()
EOSpath = dict()

phoxpath['Fe'] = '../050_phox_Fe_man/Results/phox_valsFromPDOS.csv'
refDatapath['Fe'] = '../050_phox_Fe_man/Results/phox_valsFromRefData.csv'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis:
EOSpath['Fe'] = '../081_vD_Fe_PowerLaw/Results/input_values.csv'

phoxpath['FeNi'] = '../060_phox_FeNi_man/Results/phox_valsFromPDOS.csv'
refDatapath['FeNi'] = '../060_phox_FeNi_man/Results/phox_valsFromRefData.csv'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis:
EOSpath['FeNi'] = '../091_vD_FeNi_PowerLaw/Results/input_values.csv'

phoxpath['FeNiSi'] = '../070_phox_FeNiSi_man/Results/phox_valsFromPDOS.csv'
refDatapath['FeNiSi'] = '../070_phox_FeNiSi_man/Results/phox_valsFromRefData.csv'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis:
EOSpath['FeNiSi'] = '../101_vD_FeNiSi_PowerLaw/Results/input_values.csv'


# Functions
###########


# Import data
#############

input_dict = dict()
refData_dict = dict()

for study in ['Fe','FeNi','FeNiSi']:
	phox_df = pd.read_csv(phoxpath[study], engine='python')
	refData_df = pd.read_csv(refDatapath[study], engine='python')
	EOS_df  = pd.read_csv(EOSpath[study], engine='python')

	input_df = phox_df.merge(EOS_df, on=('Folder','Index'))
	# Correct for Pandas tendency to keep trailing rounding precision vals
	input_df = input_df.round(5)
	input_df.to_csv('Results/'+study+'_input_values.csv',index=False)
	input_dict[study] = input_df

	refData_df = refData_df.merge(EOS_df, on=('Folder','Index'))
	# Correct for Pandas tendency to keep trailing rounding precision vals
	refData_df = refData_df.round(5)
	refData_df.to_csv('Results/'+study+'_refData_values.csv',index=False)
	refData_dict[study] = refData_df
 

# Plot Kinetic Energy
#####################

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,6))

markersize = 7
input_df = input_dict['Fe']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],3*bcc_df['KE'],xerr=bcc_df['dV'],yerr=3*bcc_df['dKE'],marker = 'o',
	color = 'gray', mfc='white', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label='Fe')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],3*hcp_df['KE'],xerr=hcp_df['dV'],yerr=3*hcp_df['dKE'],marker = 'o',
	color = 'gray', mfc='lightgray', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label='Fe')

input_df = input_dict['FeNi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],3*bcc_df['KE'],xerr=bcc_df['dV'],yerr=3*bcc_df['dKE'],marker = 'o',
	color = 'darkorange', mfc='white', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label=r'Fe$_{0.91}$Ni$_{0.09}$')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],3*hcp_df['KE'],xerr=hcp_df['dV'],yerr=3*hcp_df['dKE'],marker = 'o',
	color = 'darkorange', mfc='#ffc681', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label=r'Fe$_{0.91}$Ni$_{0.09}$')

input_df = input_dict['FeNiSi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],3*bcc_df['KE'],xerr=bcc_df['dV'],yerr=3*bcc_df['dKE'],marker = 'o',
	color = 'deepskyblue', mfc='white', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label=r'Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],3*hcp_df['KE'],xerr=hcp_df['dV'],yerr=3*hcp_df['dKE'],marker = 'o',
	color = 'deepskyblue', mfc='#86e1ff', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label=r'Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$')

markersize = 6
marker = '^'
refData_df = refData_dict['Fe']
bcc_df = refData_df[refData_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],3*bcc_df['KE'],xerr=bcc_df['dV'],yerr=3*bcc_df['dKE'],
	marker = marker, color = 'gray', mfc='white', ms=markersize, 
	markeredgewidth=1,ls='none',elinewidth=1, label=r'',zorder=5)
hcp_df = refData_df[refData_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],3*hcp_df['KE'],xerr=hcp_df['dV'],yerr=3*hcp_df['dKE'],
	marker = marker, color = 'gray', mfc='lightgray', ms=markersize, 
	markeredgewidth=1,ls='none',elinewidth=1, label=r'',zorder=5)

refData_df = refData_dict['FeNi']
bcc_df = refData_df[refData_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],3*bcc_df['KE'],xerr=bcc_df['dV'],yerr=3*bcc_df['dKE'],
	marker = marker, color = 'darkorange', mfc='white', ms=markersize, 
	markeredgewidth=1,ls='none',elinewidth=1, label=r'',zorder=5)
hcp_df = refData_df[refData_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],3*hcp_df['KE'],xerr=hcp_df['dV'],yerr=3*hcp_df['dKE'],
	marker = marker, color = 'darkorange', mfc='#ffc681', ms=markersize, 
	markeredgewidth=1,ls='none',elinewidth=1, label=r'',zorder=5)

refData_df = refData_dict['FeNiSi']
bcc_df = refData_df[refData_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],3*bcc_df['KE'],xerr=bcc_df['dV'],yerr=3*bcc_df['dKE'],
	marker = marker, color = 'deepskyblue', mfc='white', ms=markersize, 
	markeredgewidth=1,ls='none',elinewidth=1, label=r'',zorder=5)
hcp_df = refData_df[refData_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],3*hcp_df['KE'],xerr=hcp_df['dV'],yerr=3*hcp_df['dKE'],
	marker = marker, color = 'deepskyblue', mfc='#86e1ff', ms=markersize, 
	markeredgewidth=1,ls='none',elinewidth=1, label=r'',zorder=5)

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$E_K$ (meV/atom)',fontsize = 16)
ax0.set_xlim([15,24])
ax0.set_ylim([40.5,51])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=1,fontsize=11)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/KE.pdf', format='pdf',
	bbox_inches='tight')
plt.close()



# Plot Kinetic Energy vs P for comparison
#####################

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,6))

# Mao 2001
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Mao2001_bccFe.csv')
ax0.plot(lit_df['P'],lit_df['Ek'],'x',color='gray', mfc='white',ms=5,
	label='Fe (Mao et al. 2001)')
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Mao2001_hcpFe.csv')
ax0.plot(lit_df['P'],lit_df['Ek'],'x',color='gray', mfc='white',ms=5,
	label='')

# Gleason 2013
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Gleason2013_hcpFe.csv')
ax0.errorbar(lit_df['P'],lit_df['Ek'],xerr=lit_df['dP'],yerr=lit_df['dEk'],
	marker='^',color='gray', mfc='white',ms=5,linestyle='None',elinewidth=1,
	label='Fe (Gleason et al. 2013)')

# Lin 2005
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Lin2005_hcpFe.csv')
lit_df = lit_df[lit_df['T']==300] # Only use 300 K data
ax0.errorbar(lit_df['P'],3*lit_df['Ek'],xerr=lit_df['dP'],yerr=3*lit_df['dEk'],
	marker='d',color='gray', mfc='white',ms=5,linestyle='None',elinewidth=1,
	label='Fe (Lin et al. 2005)')

# Lin 2003a
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Lin2003_hcpFeNi.csv')
ax0.errorbar(lit_df['P'],3*lit_df['Ek'],xerr=lit_df['dP'],yerr=3*lit_df['dEk'],
	marker='s',color='darkorange', mfc='white',ms=5,linestyle='None',elinewidth=1,
	label='Fe$_{0.92}$Ni$_{0.08}$ (Lin et al. 2003a)')
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Lin2003_hcpFeSi.csv',
	header=1)
ax0.errorbar(lit_df['P'],3*lit_df['Ek'],xerr=lit_df['dP'],yerr=3*lit_df['dEk'],
	marker='s',color='deepskyblue', mfc='white',ms=5,linestyle='None',elinewidth=1,
	label='Fe$_{0.85}$Si$_{0.15}$ (Lin et al. 2003a)')

markersize = 7
input_df = input_dict['Fe']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['P'],3*bcc_df['KE'],xerr=bcc_df['dP'],yerr=3*bcc_df['dKE'],marker = 'o',
	color = 'gray', mfc='white', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label='')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['P'],3*hcp_df['KE'],xerr=hcp_df['dP'],yerr=3*hcp_df['dKE'],marker = 'o',
	color = 'gray', mfc='lightgray', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label='')

input_df = input_dict['FeNi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['P'],3*bcc_df['KE'],xerr=bcc_df['dP'],yerr=3*bcc_df['dKE'],marker = 'o',
	color = 'darkorange', mfc='white', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label=r'')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['P'],3*hcp_df['KE'],xerr=hcp_df['dP'],yerr=3*hcp_df['dKE'],marker = 'o',
	color = 'darkorange', mfc='#ffc681', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label=r'')

input_df = input_dict['FeNiSi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['P'],3*bcc_df['KE'],xerr=bcc_df['dP'],yerr=3*bcc_df['dKE'],marker = 'o',
	color = 'deepskyblue', mfc='white', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label=r'')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['P'],3*hcp_df['KE'],xerr=hcp_df['dP'],yerr=3*hcp_df['dKE'],marker = 'o',
	color = 'deepskyblue', mfc='#86e1ff', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label=r'')

ax0.set_xlabel(r'$P$ (GPa)',fontsize = 16)
ax0.set_ylabel(r'$E_K$ (meV/atom)',fontsize = 16)
# ax0.set_xlim([15,24])
ax0.set_ylim([40.5,51])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=4,fontsize=12)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/KE_vs_P.pdf', format='pdf',
	bbox_inches='tight')
plt.close()



# Plot Kinetic Energy for Murphy comparison
#####################

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,6))

# Murphy 2013
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Murphy2013_hcpFe.csv')
ax0.errorbar(lit_df['V'],lit_df['KE'],xerr=lit_df['dV'],yerr=lit_df['dKE'],
	marker='o',color='red', mfc='red',ms=7,linestyle='None',elinewidth=1,
	label='Murphy et al. 2013, as published')


markersize = 7
input_df = input_dict['Fe']
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],3*hcp_df['KE'],xerr=hcp_df['dV'],yerr=3*hcp_df['dKE'],marker = 'o',
	color = 'black', mfc='black', ms=markersize, markeredgewidth=1,
	ls='none',elinewidth=1, label='Murphy et al. 2013, re-analyzed')


ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$E_K$ (meV/atom)',fontsize = 16)
# ax0.set_xlim([15,20])
# ax0.set_ylim([40.5,51])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=1,fontsize=11)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/KE_MurphyComp.pdf', format='pdf',
	bbox_inches='tight')
plt.close()