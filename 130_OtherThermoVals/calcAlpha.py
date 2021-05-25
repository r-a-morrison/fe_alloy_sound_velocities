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

def fitSvibLinear(V,a,b):
	Svib = a*V + b  # in kB/atom
	return Svib

def calc_alphaKT_Linear(V,a,b):
	alphaKT = a*np.ones(len(V))      # in kB/(atom A^3)
	alphaKT = alphaKT*2*1.3806*10**(3)  # Convert to 10^-5 GPa/K
	return alphaKT

def calc_dalphaKT_Linear(V,cov):
	stddevs = np.sqrt(np.diag(cov))
	da = stddevs[0]
	dalphaKT = da*np.ones(len(V))
	dalphaKT = dalphaKT*2*1.3806*10**(3)  # Convert to 10^-5 GPa/K
	return dalphaKT

# def fitSvibQuad(V,a,b,c):
# 	Svib = a*V**2 + b*V + c     # in kB/atom
# 	return Svib

# def calc_alphaKT_Quad(V,a,b,c):
# 	alphaKT = 2*a*V + b               # in kB/(atom A^3)
# 	alphaKT = alphaKT*2*1.3806*10**(3)  # Convert to 10^-5 GPa/K
# 	return alphaKT

# def calc_dalphaKT_Quad(V,cov):
# 	# We only need the covariance between a and b
# 	cov = np.matrix(cov[:2,:2])
# 	# Calculate the Jacobian to propagate errors
# 	J = np.matrix([2*V, 1*np.ones(len(V))]).T
# 	# Calculate the covariance matrix for Pvib
# 	alphaKT_cov = J@(cov@J.T)
# 	# Reduce down to uncertainty for each Pvib
# 	dalphaKT = np.sqrt(np.diag(Pvib_cov))
# 	dalphaKT = dalphaKT*2*1.3806*10**(3)  # Convert to 10^-5 GPa/K
# 	return dalphaKT


# Import data
#############

input_dict = dict()
for study in ['Fe','FeNi','FeNiSi']:
	phox_df = pd.read_csv(phoxpath[study], engine='python')
	EOS_df  = pd.read_csv(EOSpath[study], engine='python')

	input_df = phox_df.merge(EOS_df, on=('Folder','Index'))
	# Correct for Pandas tendency to keep trailing rounding precision vals
	input_df = input_df.round(5)
	input_df.to_csv('Results/'+study+'_Svib_input_values.csv',index=False)
	input_dict[study] = input_df


# Fit Svib with linear and polynomial function accounting for uncertainties
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

	lin_fit[study], lin_cov[study] = curve_fit(fitSvibLinear,input_df['V'],
		input_df['Svib'], sigma=input_df['dSvib'])
	# quad_fit[study], quad_cov[study] = curve_fit(fitSvibQuad,input_df['V'],
	# 	input_df['Svib'], sigma=input_df['dSvib'])

	print(lin_fit[study])
	print(lin_cov[study])
	std_fit = np.sqrt(np.diag(lin_cov[study]))
	print(std_fit)
	corr = np.diag(std_fit**(-1))@lin_cov[study]@np.diag(std_fit**(-1))
	print(corr)

	# print(quad_fit[study])
	# print(quad_cov[study])
	# std_fit = np.sqrt(np.diag(quad_cov[study]))
	# print(std_fit)
	# corr = np.diag(std_fit**(-1))@quad_cov[study]@np.diag(std_fit**(-1))
	# print(corr)


# Calculate alpha
#################
for study in ['Fe','FeNi','FeNiSi']:
	input_df = input_dict[study]

	alphaKT = calc_alphaKT_Linear(input_df['V'],*lin_fit[study])
	dalphaKT = calc_dalphaKT_Linear(input_df['V'],lin_cov['Fe'])
	input_df['alpha'] = alphaKT/input_df['KT']
	input_df['dalpha'] = np.sqrt(( alphaKT*input_df['KT']**(-2)*input_df['dKT'] )**2 +
		( dalphaKT/input_df['KT'] )**2)

	input_dict[study] = input_df
	input_df.to_csv('Results/'+study+'_alphavib.csv',index=False)
 

# Plot Svib
###########

V_array = np.arange(15,21,0.1)
# V_conv = 2*1E24/6.022E23
# MurphyResults = [1.67/V_conv**2,-43.74/V_conv,220.15]

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,6))

input_df = input_dict['Fe']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],bcc_df['Svib'],xerr=bcc_df['dV'],yerr=bcc_df['dSvib'],marker = 'o',
	color = 'gray', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label='Fe')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['Svib'],xerr=hcp_df['dV'],yerr=hcp_df['dSvib'],marker = 'o',
	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label='Fe')
ax0.plot(V_array,fitSvibLinear(V_array,*lin_fit['Fe']),
	'-',color='gray',linewidth=1)
# ax0.plot(V_array,fitSvibQuad(V_array,*quad_fit['Fe']),
# 	'-',color='gray',linewidth=1,label='fit')
# ax0.plot(V_array,fitSvibQuad(V_array,*MurphyResults),
# 	'--',color='gray',linewidth=1,label='fit')

input_df = input_dict['FeNi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],bcc_df['Svib'],xerr=bcc_df['dV'],yerr=bcc_df['dSvib'],marker = 'o',
	color = 'darkorange', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'Fe$_{0.91}$Ni$_{0.09}$')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['Svib'],xerr=hcp_df['dV'],yerr=hcp_df['dSvib'],marker = 'o',
	color = 'darkorange', mfc='#ffc681', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'Fe$_{0.91}$Ni$_{0.09}$')
ax0.plot(V_array,fitSvibLinear(V_array,*lin_fit['FeNi']),
	'-',color='darkorange',linewidth=1)
# ax0.plot(V_array,fitSvibQuad(V_array,*quad_fit['FeNi']),
# 	'-',color='darkorange',linewidth=1,label='fit')

input_df = input_dict['FeNiSi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],bcc_df['Svib'],xerr=bcc_df['dV'],yerr=bcc_df['dSvib'],marker = 'o',
	color = 'deepskyblue', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['Svib'],xerr=hcp_df['dV'],yerr=hcp_df['dSvib'],marker = 'o',
	color = 'deepskyblue', mfc='#86e1ff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$')
ax0.plot(V_array,fitSvibLinear(V_array,*lin_fit['FeNiSi']),
	'-',color='deepskyblue', linewidth=1)
# ax0.plot(V_array,fitSvibQuad(V_array,*quad_fit['FeNiSi']),
# 	'-',color='deepskyblue', linewidth=1,label='fit')

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$S_{vib}$ ($k_B$/atom)',fontsize = 16)
# ax0.set_xlim([15,21])
ax0.set_ylim([1.6,3.3])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=4,fontsize=11)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Svib.pdf', format='pdf',
	bbox_inches='tight')
plt.close()



# Plot Svib vs P for comparison
###############################

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,6))

# Mao 2001
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Mao2001_bccFe.csv')
ax0.plot(lit_df['P'],lit_df['Svib'],'x',color='gray', mfc='white',ms=5,
	label='Fe (Mao et al. 2001)')
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Mao2001_hcpFe.csv')
ax0.plot(lit_df['P'],lit_df['Svib'],'x',color='gray', mfc='white',ms=5,
	label='')

# Gleason 2013
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Gleason2013_hcpFe.csv')
ax0.errorbar(lit_df['P'],lit_df['Svib'],xerr=lit_df['dP'],yerr=lit_df['dSvib'],
	marker='^',color='gray', mfc='white',ms=5,linestyle='None',elinewidth=1,
	label='Fe (Gleason et al. 2013)')

# Lin 2005
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Lin2005_hcpFe.csv')
lit_df = lit_df[lit_df['T']==300] # Only use 300 K data
ax0.errorbar(lit_df['P'],lit_df['Svib'],xerr=lit_df['dP'],yerr=lit_df['dSvib'],
	marker='d',color='gray', mfc='white',ms=5,linestyle='None',elinewidth=1,
	label='Fe (Lin et al. 2005)')

# Lin 2003a
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Lin2003_hcpFeNi.csv')
ax0.errorbar(lit_df['P'],lit_df['Svib'],xerr=lit_df['dP'],yerr=lit_df['dSvib'],
	marker='s',color='darkorange', mfc='white',ms=5,linestyle='None',elinewidth=1,
	label='Fe$_{0.92}$Ni$_{0.08}$ (Lin et al. 2003a)')
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Lin2003_hcpFeSi.csv',
	header=1)
ax0.errorbar(lit_df['P'],lit_df['Svib'],xerr=lit_df['dP'],yerr=lit_df['dSvib'],
	marker='s',color='deepskyblue', mfc='white',ms=5,linestyle='None',elinewidth=1,
	label='Fe$_{0.85}$Si$_{0.15}$ (Lin et al. 2003a)')


input_df = input_dict['Fe']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['P'],bcc_df['Svib'],xerr=bcc_df['dP'],
	yerr=bcc_df['dSvib'],marker = 'o',
	color = 'gray', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label='')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['P'],hcp_df['Svib'],xerr=hcp_df['dP'],
	yerr=hcp_df['dSvib'],marker = 'o',
	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label='')

input_df = input_dict['FeNi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['P'],bcc_df['Svib'],xerr=bcc_df['dP'],
	yerr=bcc_df['dSvib'],marker = 'o',
	color = 'darkorange', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['P'],hcp_df['Svib'],xerr=hcp_df['dP'],
	yerr=hcp_df['dSvib'],marker = 'o',
	color = 'darkorange', mfc='#ffc681', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'')

input_df = input_dict['FeNiSi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['P'],bcc_df['Svib'],xerr=bcc_df['dP'],
	yerr=bcc_df['dSvib'],marker = 'o',
	color = 'deepskyblue', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['P'],hcp_df['Svib'],xerr=hcp_df['dP'],
	yerr=hcp_df['dSvib'],marker = 'o',
	color = 'deepskyblue', mfc='#86e1ff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'')

ax0.set_xlabel(r'$P$ (GPa)',fontsize = 16)
ax0.set_ylabel(r'$S_{vib}$ ($k_B$/atom)',fontsize = 16)
# ax0.set_xlim([15,21])
ax0.set_ylim([1.6,3.3])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=1,fontsize=11)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Svib_vs_P.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot Svib for Murphy comparison
###########

V_array = np.arange(15,21,0.1)
# V_conv = 2*1E24/6.022E23
# MurphyResults = [1.67/V_conv**2,-43.74/V_conv,220.15]

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,6))

# Murphy 2013
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Murphy2013_hcpFe.csv')
ax0.errorbar(lit_df['V'],lit_df['Svib'],xerr=lit_df['dV'],yerr=lit_df['dSvib'],
	marker='o',color='red', mfc='red',ms=7,linestyle='None',elinewidth=1,
	label='Murphy et al. 2013, as published')

input_df = input_dict['Fe']
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['Svib'],xerr=hcp_df['dV'],yerr=hcp_df['dSvib'],marker = 'o',
	color = 'black', mfc='black', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label='Murphy et al. 2013, re-analyzed')
# ax0.plot(V_array,fitSvibLinear(V_array,*lin_fit['Fe']),
# 	'-',color='black',linewidth=1)

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$S_{vib}$ ($k_B$/atom)',fontsize = 16)
# ax0.set_xlim([15,21])
# ax0.set_ylim([1.6,2.7])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=4,fontsize=11)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Svib_MurphyComp.pdf', format='pdf',
	bbox_inches='tight')
plt.close()



# Plot alpha*KT
###############

V_array = np.arange(15,21,0.1)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

input_df = input_dict['Fe']
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*fit['Fe']),
# 	'-',color='gray',linewidth=1)
# ax0.plot(V_array,calcPvibQuad(V_array,*MurphyResults),
# 	'--',color='gray',linewidth=1)
ax0.errorbar(hcp_df['V'],calc_alphaKT_Linear(hcp_df['V'],*lin_fit['Fe']),
	yerr=calc_dalphaKT_Linear(hcp_df['V'],lin_cov['Fe']),marker = 'o',
	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label='hcp Fe')

input_df = input_dict['FeNi']
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*fit['FeNi']),
# 	'-',color='darkorange',linewidth=1)
ax0.errorbar(hcp_df['V'],calc_alphaKT_Linear(hcp_df['V'],*lin_fit['FeNi']),
	yerr=calc_dalphaKT_Linear(hcp_df['V'],lin_cov['FeNi']),marker = 'o',
	color = 'darkorange', mfc='#ffc681', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label=r'hcp Fe$_{0.91}$Ni$_{0.09}$')

input_df = input_dict['FeNiSi']
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibLinear(V_array,*fit['FeNiSi']),
# 	'-',color='deepskyblue',linewidth=1)
ax0.errorbar(hcp_df['V'],calc_alphaKT_Linear(hcp_df['V'],*lin_fit['FeNiSi']),
	yerr=calc_dalphaKT_Linear(hcp_df['V'],lin_cov['FeNiSi']),marker = 'o',
	color = 'deepskyblue', mfc='#86e1ff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label=r'hcp Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$')

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$\alpha K_T$ (10$^{-5}$ GPa/K)',fontsize = 16)
# ax0.set_xlim([15,23])
# ax0.set_ylim([1.5,3.5])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=2,fontsize=13)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/alphaKT.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot alpha
############

V_array = np.arange(15,21,0.1)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

input_df = input_dict['Fe']
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['alpha'],yerr=hcp_df['dalpha'],marker = 'o',
	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label='hcp Fe (Murphy et al. 2013*)')

input_df = input_dict['FeNi']
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['alpha'],yerr=hcp_df['dalpha'],marker = 'o',
	color = 'darkorange', mfc='#ffc681', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label=r'hcp Fe$_{0.91}$Ni$_{0.09}$ (this study)')

input_df = input_dict['FeNiSi']
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['alpha'],yerr=hcp_df['dalpha'],marker = 'o',
	color = 'deepskyblue', mfc='#86e1ff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label=r'hcp Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$ (this study)')

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$\alpha_{vib}$ (10$^{-5}$ K$^{-1}$)',fontsize = 16)
ax0.set_xlim([15,21])
# ax0.set_ylim([0.5,2.4])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=2,fontsize=13)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/alpha.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot alpha for Murphy comparison
###########

V_array = np.arange(15,21,0.1)
# V_conv = 2*1E24/6.022E23
# MurphyResults = [1.67/V_conv**2,-43.74/V_conv,220.15]

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(6,6))

# Murphy 2013
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Murphy2013_hcpFe.csv')
ax0.errorbar(lit_df['V'],lit_df['alphavib'],xerr=lit_df['dV'],yerr=lit_df['dalphavib'],
	marker='o',color='red', mfc='red',ms=7,linestyle='None',elinewidth=1,
	label='Murphy et al. 2013, as published')

input_df = input_dict['Fe']
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['alpha'],xerr=hcp_df['dV'],yerr=hcp_df['dalpha'],marker = 'o',
	color = 'black', mfc='black', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label='Murphy et al. 2013, re-analyzed')
# ax0.plot(V_array,fitSvibLinear(V_array,*lin_fit['Fe']),
# 	'-',color='black',linewidth=1)

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$\alpha_{vib}$ (10$^{-5}$ K$^{-1}$)',fontsize = 16)
# ax0.set_xlim([15,21])
# ax0.set_ylim([1.6,2.7])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=4,fontsize=11)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/alpha_MurphyComp.pdf', format='pdf',
	bbox_inches='tight')
plt.close()