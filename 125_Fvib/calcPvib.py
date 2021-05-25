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
V0_dict = dict()

phoxpath['Fe'] = '../050_phox_Fe_man/Results/phox_valsFromPDOS.csv'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis:
EOSpath['Fe'] = '../081_vD_Fe_PowerLaw/Results/input_values.csv'
V0_dict['Fe'] = 22.428

phoxpath['FeNi'] = '../060_phox_FeNi_man/Results/phox_valsFromPDOS.csv'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis:
EOSpath['FeNi'] = '../091_vD_FeNi_PowerLaw/Results/input_values.csv'
V0_dict['FeNi'] = 22.505

phoxpath['FeNiSi'] = '../070_phox_FeNiSi_man/Results/phox_valsFromPDOS.csv'
# All the relevant info (P, V, etc.) was nicely collected during the sound 
# velocity analysis:
EOSpath['FeNiSi'] = '../101_vD_FeNiSi_PowerLaw/Results/input_values.csv'
V0_dict['FeNiSi'] = 22.952


# Functions
###########

def fitFvibLinear(V,a,b):
	Fvib = a*V + b
	return Fvib

def calcPvibLinear(V,a,b):
	Pvib = -a*np.ones(len(V))  # in mev/(atom A^3)
	unitconv = 2*1.6022*10**(-1)
	Pvib = Pvib*unitconv  # Unit conversion to GPa
	print('P_vib = '+str(-a*unitconv))
	# print('P_vib = '+str(unitconv)+'('+str(-a)+')')
	return Pvib

def calcdPvibLinear(V,cov):
	stddevs = np.sqrt(np.diag(cov))
	da = stddevs[0]
	dPvib = da*np.ones(len(V))
	dPvib = dPvib*2*1.6022*10**(-1)  # Unit conversion to GPa
	return dPvib

def fitFvibQuad(V,a,b,c):
	Fvib = a*V**2 + b*V + c
	return Fvib

def calcPvibQuad(V,a,b,c):
	Pvib = -(2*a*V + b)             # in mev/(atom A^3)
	unitconv = 2*1.6022*10**(-1)
	Pvib = Pvib*unitconv  # Unit conversion to GPa
	print('P_vib = '+str(-2*a*unitconv)+'V + '+str(-b*unitconv))
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

def calcPanh(V,V0,T):
	k = constants.k
	m = 1.87
	a0 = 3.7E-5  # K^-1
	Panh = ((3*k)/(2*V))*m*a0*(V/V0)**m*T**2
	Panh = Panh*2E21 # unit conversion to GPa
	return Panh

def calcPel(V,V0,T):
	k = constants.k
	g = 1.339
	e0 = 1.95E-4  # K^-1
	Pel = ((3*k)/(2*V))*g*e0*(V/V0)**g*T**2
	Pel = Pel*2E21 # unit conversion to GPa
	return Pel


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

fit = dict()
cov = dict()
hcp_dict = dict()

for study in ['Fe','FeNi','FeNiSi']:
	print(study+':')
	input_df = input_dict[study]

	# Only fit hcp data
	hcp_df = input_df[input_df['Phase']=='hcp'].copy()

	fit[study], cov[study] = curve_fit(fitFvibQuad,hcp_df['V'],
		hcp_df['Fvib'])
	hcp_df['Pvib'] = calcPvibQuad(hcp_df['V'],*fit[study])
	hcp_df['dPvib'] = calcdPvibQuad(hcp_df['V'],cov[study])

	# if (study=='Fe') or (study=='FeNi'):
	# 	fit[study], cov[study] = curve_fit(fitFvibQuad,hcp_df['V'],
	# 		hcp_df['Fvib'], sigma=hcp_df['dFvib'])
	# 	hcp_df['Pvib'] = calcPvibQuad(hcp_df['V'],*fit[study])
	# 	hcp_df['dPvib'] = calcdPvibQuad(hcp_df['V'],cov[study])
		
	# elif study=='FeNiSi':
	# 	fit[study], cov[study] = curve_fit(fitFvibLinear,hcp_df['V'],
	# 		hcp_df['Fvib'], sigma=hcp_df['dFvib'])
	# 	hcp_df['Pvib'] = calcPvibLinear(hcp_df['V'],*fit[study])
	# 	hcp_df['dPvib'] = calcdPvibLinear(hcp_df['V'],cov[study])

	hcp_df['Panh_5500K'] = calcPanh(hcp_df['V'],V0_dict[study],5500)
	hcp_df['Pel_5500K']  = calcPel( hcp_df['V'],V0_dict[study],5500)
	
	print('Fvib fit parameters: AV^2 + BV + C')
	print(fit[study])
	print('Fvib fit covariance:')
	print(cov[study])
	print('Fvib fit parameter uncertainty:')
	std_fit = np.sqrt(np.diag(cov[study]))
	print(std_fit)
	print('Fvib fit parameter correlation:')
	corr = np.linalg.inv(np.diag(std_fit))@cov[study]@np.linalg.inv(np.diag(std_fit))
	print(corr)

	hcp_dict[study] = hcp_df
	hcp_df = hcp_df.round(4)
	hcp_df.to_csv('Results/'+study+'_Pth.csv',index=False)


# Plot Fvib
###########

V_array = np.arange(15,21,0.1)
V_conv = 2*1E24/6.022E23
MurphyResults = [1.67/V_conv**2,-43.74/V_conv,220.15]

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

input_df = input_dict['Fe']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],bcc_df['Fvib'],yerr=bcc_df['dFvib'],marker = 'o',
	color = 'gray', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label='Fe')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['Fvib'],yerr=hcp_df['dFvib'],marker = 'o',
	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label='Fe')
# ax0.plot(V_array,fitFvibLinear(V_array,*lin_fit['Fe']),
# 	'--',color='gray',linewidth=1)
ax0.plot(V_array,fitFvibQuad(V_array,*fit['Fe']),
	'-',color='gray',linewidth=1,label='fit')
# ax0.plot(V_array,fitFvibQuad(V_array,*MurphyResults),
# 	'--',color='gray',linewidth=1,label='fit')

input_df = input_dict['FeNi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],bcc_df['Fvib'],yerr=bcc_df['dFvib'],marker = 'o',
	color = 'darkorange', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'Fe$_{0.91}$Ni$_{0.09}$')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['Fvib'],yerr=hcp_df['dFvib'],marker = 'o',
	color = 'darkorange', mfc='#ffc681', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'Fe$_{0.91}$Ni$_{0.09}$')
# ax0.plot(V_array,fitFvibLinear(V_array,*lin_fit['FeNi']),
# 	'--',color='darkorange',linewidth=1)
ax0.plot(V_array,fitFvibQuad(V_array,*fit['FeNi']),
	'-',color='darkorange',linewidth=1,label='fit')

input_df = input_dict['FeNiSi']
bcc_df = input_df[input_df['Phase']=='bcc']
ax0.errorbar(bcc_df['V'],bcc_df['Fvib'],yerr=bcc_df['dFvib'],marker = 'o',
	color = 'deepskyblue', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$')
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['Fvib'],yerr=hcp_df['dFvib'],marker = 'o',
	color = 'deepskyblue', mfc='#86e1ff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label=r'Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$')
# ax0.plot(V_array,fitFvibLinear(V_array,*lin_fit['FeNiSi']),
# 	'--',color='deepskyblue', linewidth=1)
ax0.plot(V_array,fitFvibQuad(V_array,*fit['FeNiSi']),
	'-',color='deepskyblue', linewidth=1,label='fit')

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$F_{vib}$ (meV/atom)',fontsize = 16)
# ax0.set_xlim([15,21])
# ax0.set_ylim([0,75])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=1,fontsize=11)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Fvib.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot Fvib compare Murphy
###########

V_array = np.arange(15,21,0.1)
V_conv = 2*1E24/6.022E23
MurphyResults = [1.67/V_conv**2,-43.74/V_conv,220.15]

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

# Murphy 2013
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Murphy2013_hcpFe.csv')
ax0.errorbar(lit_df['V'],lit_df['Fvib'],xerr=lit_df['dV'],yerr=lit_df['dFvib'],
	marker='o',color='red', mfc='red',ms=7,linestyle='None',elinewidth=1,
	label='Murphy et al. 2011a, as published')

input_df = input_dict['Fe']
hcp_df = input_df[input_df['Phase']=='hcp']
ax0.errorbar(hcp_df['V'],hcp_df['Fvib'],yerr=hcp_df['dFvib'],marker = 'o',
	color = 'black', mfc='black', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	label='Murphy et al. 2011a, re-analyzed')
# ax0.plot(V_array,fitFvibLinear(V_array,*lin_fit['Fe']),
# 	'--',color='gray',linewidth=1)
# ax0.plot(V_array,fitFvibQuad(V_array,*fit['Fe']),
# 	'-',color='gray',linewidth=1,label='fit')
# ax0.plot(V_array,fitFvibQuad(V_array,*MurphyResults),
# 	'--',color='gray',linewidth=1,label='fit')

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$F_{vib}$ (meV/atom)',fontsize = 16)
# ax0.set_xlim([15,21])
# ax0.set_ylim([0,75])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=1,fontsize=11)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Fvib_MurphyComp.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot Pvib
###########

V_array = np.arange(15,21,0.1)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

input_df = input_dict['Fe']
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*fit['Fe']),
# 	'-',color='gray',linewidth=1)
# ax0.plot(V_array,calcPvibQuad(V_array,*MurphyResults),
# 	'--',color='gray',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*fit['Fe']),
	yerr=calcdPvibQuad(hcp_df['V'],cov['Fe']),marker = 'o',
	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label='hcp-Fe')

input_df = input_dict['FeNi']
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*fit['FeNi']),
# 	'-',color='darkorange',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*fit['FeNi']),
	yerr=calcdPvibQuad(hcp_df['V'],cov['FeNi']),marker = 'o',
	color = 'darkorange', mfc='#ffc681', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label=r'hcp-Fe$_{0.91}$Ni$_{0.09}$')

input_df = input_dict['FeNiSi']
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibLinear(V_array,*fit['FeNiSi']),
# 	'-',color='deepskyblue',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*fit['FeNiSi']),
	yerr=calcdPvibQuad(hcp_df['V'],cov['FeNiSi'])*(2/3),marker = 'o',
	color = 'deepskyblue', mfc='#86e1ff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label=r'hcp-Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$')

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$P^h_{vib}$ (GPa)',fontsize = 16)
# ax0.set_xlim([15,23])
# ax0.set_ylim([1.5,3.5])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=1,fontsize=13)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Pvib.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot Pvib Murphy compare
###########

V_array = np.arange(15,21,0.1)

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

# Murphy 2013
lit_df = pd.read_csv('../006_Literature/ExpWithSameUnits/Murphy2013_hcpFe.csv')
ax0.errorbar(lit_df['V'],lit_df['Pvib'],xerr=lit_df['dV'],yerr=lit_df['dPvib'],
	marker='o',color='red', mfc='red',ms=7,linestyle='None',elinewidth=1,
	label='Murphy et al. 2013, as published')

input_df = input_dict['Fe']
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*fit['Fe']),
# 	'-',color='gray',linewidth=1)
# ax0.plot(V_array,calcPvibQuad(V_array,*MurphyResults),
# 	'--',color='gray',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*fit['Fe']),
	yerr=calcdPvibQuad(hcp_df['V'],cov['Fe']),marker = 'o',
	color = 'black', mfc='black', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label='Murphy et al. 2013, re-analyzed')

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$P^h_{vib}$ (GPa)',fontsize = 16)
# ax0.set_xlim([15,23])
# ax0.set_ylim([1.5,3.5])
ax0.tick_params(direction='in',right='on',top='on')

ax0.legend(loc=1,fontsize=13)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Pvib_MurphyComp.pdf', format='pdf',
	bbox_inches='tight')
plt.close()


# Plot Pvib at high T
#####################

V_array = np.arange(15,21,0.1)
T0 = 300
T1 = 2000
T2 = 4000
T3 = 5500
T4 = 6500

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

# Plot Alfe 2001's Fvib harmonic results
Alfe_2000K = pd.read_csv('../006_Literature/ExpWithSameUnits/Alfe2001_2000K.csv')
ax0.plot(Alfe_2000K['V'],Alfe_2000K['Pvib'],color='gray',linewidth=1,linestyle='-')
Alfe_4000K = pd.read_csv('../006_Literature/ExpWithSameUnits/Alfe2001_4000K.csv')
ax0.plot(Alfe_4000K['V'],Alfe_4000K['Pvib'],color='gray',linewidth=1,linestyle='-')
Alfe_6000K = pd.read_csv('../006_Literature/ExpWithSameUnits/Alfe2001_6000K.csv')
ax0.plot(Alfe_6000K['V'],Alfe_6000K['Pvib'],color='gray',linewidth=1,linestyle='-')

input_df = input_dict['Fe']
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*fit['Fe']),
# 	'-',color='gray',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*fit['Fe']),
	yerr=calcdPvibQuad(hcp_df['V'],cov['Fe']),marker = 'o',
	color = 'gray', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label='hcp Fe')
ax0.errorbar(hcp_df['V'],(T1/T0)*calcPvibQuad(hcp_df['V'],*fit['Fe']),
	yerr=(T1/T0)*calcdPvibQuad(hcp_df['V'],cov['Fe']),marker = 'o',
	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T2/T0)*calcPvibQuad(hcp_df['V'],*fit['Fe']),
	yerr=(T2/T0)*calcdPvibQuad(hcp_df['V'],cov['Fe']),marker = 'o',
	color = 'gray', mfc='darkgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T3/T0)*calcPvibQuad(hcp_df['V'],*fit['Fe']),
	yerr=(T3/T0)*calcdPvibQuad(hcp_df['V'],cov['Fe']),marker = 'o',
	color = 'gray', mfc='gray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)


input_df = input_dict['FeNi']
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*fit['FeNi']),
# 	'-',color='darkorange',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*fit['FeNi']),
	yerr=calcdPvibQuad(hcp_df['V'],cov['FeNi']),marker = 'o',
	color = 'darkorange', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label=r'hcp Fe$_{0.91}$Ni$_{0.09}$')
ax0.errorbar(hcp_df['V'],(T1/T0)*calcPvibQuad(hcp_df['V'],*fit['FeNi']),
	yerr=(T1/T0)*calcdPvibQuad(hcp_df['V'],cov['FeNi']),marker = 'o',
	color = 'darkorange', mfc='#ffd6a9', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T2/T0)*calcPvibQuad(hcp_df['V'],*fit['FeNi']),
	yerr=(T2/T0)*calcdPvibQuad(hcp_df['V'],cov['FeNi']),marker = 'o',
	color = 'darkorange', mfc='#ffac52', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T3/T0)*calcPvibQuad(hcp_df['V'],*fit['FeNi']),
	yerr=(T3/T0)*calcdPvibQuad(hcp_df['V'],cov['FeNi']),marker = 'o',
	color = 'darkorange', mfc='darkorange', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)


input_df = input_dict['FeNiSi']
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibLinear(V_array,*fit['FeNiSi']),
# 	'-',color='deepskyblue',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*fit['FeNiSi']),
	yerr=calcdPvibQuad(hcp_df['V'],cov['FeNiSi'])*(2/3),marker = 'o',
	color = 'deepskyblue', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label=r'hcp Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$')
ax0.errorbar(hcp_df['V'],(T1/T0)*calcPvibQuad(hcp_df['V'],*fit['FeNiSi']),
	yerr=(T1/T0)*calcdPvibQuad(hcp_df['V'],cov['FeNiSi'])*(2/3),marker = 'o',
	color = 'deepskyblue', mfc='#aeeeff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T2/T0)*calcPvibQuad(hcp_df['V'],*fit['FeNiSi']),
	yerr=(T2/T0)*calcdPvibQuad(hcp_df['V'],cov['FeNiSi'])*(2/3),marker = 'o',
	color = 'deepskyblue', mfc='#57dbff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T3/T0)*calcPvibQuad(hcp_df['V'],*fit['FeNiSi']),
	yerr=(T3/T0)*calcdPvibQuad(hcp_df['V'],cov['FeNiSi'])*(2/3),marker = 'o',
	color = 'deepskyblue', mfc='deepskyblue', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$P^h_{vib}$ (GPa)',fontsize = 16)
ax0.set_xlim([14.5,21.5])
ax0.set_ylim([0,60])
ax0.tick_params(direction='in',right='on',top='on')
ax0.xaxis.set_minor_locator(AutoMinorLocator(5))
ax0.yaxis.set_minor_locator(AutoMinorLocator(5))
ax0.tick_params(direction='in',which='minor',right='on',top='on')


# ax0.legend(loc=1,fontsize=11)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Pvib_highT.pdf', format='pdf',
	bbox_inches='tight')
plt.close()



# Plot Pth at high T
#####################

V_array = np.arange(15,21,0.1)
T0 = 300
T1 = 2000
T2 = 4000
T3 = 5500
T4 = 6500

fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(7.5,6))

study = 'Fe'
input_df = input_dict[study]
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*fit['Fe']),
# 	'-',color='gray',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T0) + 
	calcPel(hcp_df['V'],V0_dict[study],T0),
	yerr=calcdPvibQuad(hcp_df['V'],cov[study]),marker = 'o',
	color = 'gray', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label='hcp Fe')
ax0.errorbar(hcp_df['V'],(T1/T0)*calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T1) + 
	calcPel(hcp_df['V'],V0_dict[study],T1),
	yerr=(T1/T0)*calcdPvibQuad(hcp_df['V'],cov[study]),marker = 'o',
	color = 'gray', mfc='lightgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T2/T0)*calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T2) + 
	calcPel(hcp_df['V'],V0_dict[study],T2),
	yerr=(T2/T0)*calcdPvibQuad(hcp_df['V'],cov[study]),marker = 'o',
	color = 'gray', mfc='darkgray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T3/T0)*calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T3) + 
	calcPel(hcp_df['V'],V0_dict[study],T3),
	yerr=(T3/T0)*calcdPvibQuad(hcp_df['V'],cov[study]),marker = 'o',
	color = 'gray', mfc='gray', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)

study = 'FeNi'
input_df = input_dict[study]
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibQuad(V_array,*fit['FeNi']),
# 	'-',color='darkorange',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T0) + 
	calcPel(hcp_df['V'],V0_dict[study],T0),
	yerr=calcdPvibQuad(hcp_df['V'],cov[study]),marker = 'o',
	color = 'darkorange', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label=r'hcp Fe$_{0.91}$Ni$_{0.09}$')
ax0.errorbar(hcp_df['V'],(T1/T0)*calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T1) + 
	calcPel(hcp_df['V'],V0_dict[study],T1),
	yerr=(T1/T0)*calcdPvibQuad(hcp_df['V'],cov[study]),marker = 'o',
	color = 'darkorange', mfc='#ffd6a9', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T2/T0)*calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T2) + 
	calcPel(hcp_df['V'],V0_dict[study],T2),
	yerr=(T2/T0)*calcdPvibQuad(hcp_df['V'],cov[study]),marker = 'o',
	color = 'darkorange', mfc='#ffac52', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T3/T0)*calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T3) + 
	calcPel(hcp_df['V'],V0_dict[study],T3),
	yerr=(T3/T0)*calcdPvibQuad(hcp_df['V'],cov[study]),marker = 'o',
	color = 'darkorange', mfc='darkorange', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)

study = 'FeNiSi'
input_df = input_dict[study]
hcp_df = input_df[input_df['Phase']=='hcp']
# ax0.plot(V_array,calcPvibLinear(V_array,*fit['FeNiSi']),
# 	'-',color='deepskyblue',linewidth=1)
ax0.errorbar(hcp_df['V'],calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T0) + 
	calcPel(hcp_df['V'],V0_dict[study],T0),
	yerr=calcdPvibQuad(hcp_df['V'],cov[study])*(2/3),marker = 'o',
	color = 'deepskyblue', mfc='white', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1, label=r'hcp Fe$_{0.8}$Ni$_{0.1}$Si$_{0.1}$')
ax0.errorbar(hcp_df['V'],(T1/T0)*calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T1) + 
	calcPel(hcp_df['V'],V0_dict[study],T1),
	yerr=(T1/T0)*calcdPvibQuad(hcp_df['V'],cov[study])*(2/3),marker = 'o',
	color = 'deepskyblue', mfc='#aeeeff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T2/T0)*calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T2) + 
	calcPel(hcp_df['V'],V0_dict[study],T2),
	yerr=(T2/T0)*calcdPvibQuad(hcp_df['V'],cov[study])*(2/3),marker = 'o',
	color = 'deepskyblue', mfc='#57dbff', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)
ax0.errorbar(hcp_df['V'],(T3/T0)*calcPvibQuad(hcp_df['V'],*fit[study]) +
	calcPanh(hcp_df['V'],V0_dict[study],T3) + 
	calcPel(hcp_df['V'],V0_dict[study],T3),
	yerr=(T3/T0)*calcdPvibQuad(hcp_df['V'],cov[study])*(2/3),marker = 'o',
	color = 'deepskyblue', mfc='deepskyblue', ms=7, markeredgewidth=1,ls='none',elinewidth=1,
	capsize = 1)

ax0.set_xlabel(r'$V$ ($\AA^3$)',fontsize = 16)
ax0.set_ylabel(r'$P_{th}$ (GPa)',fontsize = 16)
ax0.set_xlim([15,21.5])
ax0.set_ylim([0,80])
ax0.tick_params(direction='in',right='on',top='on')
ax0.xaxis.set_minor_locator(AutoMinorLocator(5))
ax0.yaxis.set_minor_locator(AutoMinorLocator(5))
ax0.tick_params(direction='in',which='minor',right='on',top='on')

# ax0.legend(loc=1,fontsize=11)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(wspace=0)

fig.savefig('Plots/Pth_highT.pdf', format='pdf',
	bbox_inches='tight')
plt.close()