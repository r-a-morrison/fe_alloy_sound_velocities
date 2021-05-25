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

V0_dict = dict()
dV0_dict = dict()
M_dict = dict()

# Fe V0 from Dewaele et al. 2006
V0_dict['Fe'] = 22.428
dV0_dict['Fe'] = 0.098
M_dict['Fe'] = 56.942

# Fe0.91Ni0.09 V0 from our EOS study
V0_dict['FeNi'] = 22.505
dV0_dict['FeNi'] = 0.042
M_dict['FeNi'] = 57.100

# Fe0.8Ni0.1Si0.1 V0 from our EOS study
V0_dict['FeNiSi'] = 22.952
dV0_dict['FeNiSi'] = 0.072
M_dict['FeNiSi'] = 54.100

Fe_input_filename = 'Fe/Results/input_values.csv'
FeNi_input_filename = 'FeNi/Results/input_values.csv'
FeNiSi_input_filename = 'FeNiSi/Results/input_values.csv'

Fe_phox_filename = '../050_phox_Fe_man/Results/phox_valsFromPDOS.csv'
FeNi_phox_filename = '../060_phox_FeNi_man/Results/phox_valsFromPDOS.csv'
FeNiSi_phox_filename = '../070_phox_FeNiSi_man/Results/phox_valsFromPDOS.csv'


# Functions
###########

def calcGruneisen(V,dV,V0,dV0,gamma0,dgamma0,q):
	gamma = gamma0*(V/V0)**q
	dgamma_dgamma0 = (V/V0)**q
	dgamma_dV = (gamma0*q/V)*(V/V0)**q
	dgamma_dV0 = -(gamma0*q/V0)*(V/V0)**q
	dgamma = np.sqrt( (dgamma_dgamma0*dgamma0)**2 + (dgamma_dV*dV)**2 + (dgamma_dV*dV0)**2 )
	return gamma, dgamma

def calcPvib(V,dV,gammavib,dgammavib,Uvib,dUvib,Cvib,dCvib,Cel):
	Pvib = (Cvib*gammavib/(Cvib+Cel))*(Uvib/V)
	dP_dCvib = Cel*gammavib*Uvib/((Cvib+Cel)**2*V)
	dP_dgamma = (Cvib/(Cvib+Cel))*(Uvib/V)
	dP_dU = (Cvib*gammavib/(Cvib+Cel))*(1/V)
	dP_dV = -(Cvib*gammavib/(Cvib+Cel))*(Uvib/V**2)
	dPvib = np.sqrt( (dP_dCvib*dCvib)**2 + (dP_dgamma*dgammavib)**2 +
		(dP_dU*dUvib)**2 + (dP_dV*dV)**2)
	return Pvib, dPvib

def calcCel(V,V0,T,M):
	# From Fei et al. 2016
	beta0 = 0.07 # in J/(kg K^2)
	k = 1.34
	NA = 6.022141*10**23 # mol^(-1)
	kB = 1.38065*10**(-23) # J/K
	beta0 = beta0*M/(10**3*kB*NA) # in kB/(atom K)
	Cel = beta0*(V/V0)**k*T # in kB/atom
	return Cel


# Combine input and phox datasets
#################################

input_dict = dict()

input_df = pd.read_csv(Fe_input_filename)
phox_df = pd.read_csv(Fe_phox_filename)
# Only use hcp phases
input_df = input_df[input_df['Phase']=='hcp']
input_df = input_df[['Folder','Index','Phase','V','dV','P','dP']]
# Combine input and phox dataframes
input_df = input_df.merge(phox_df,on='Index')
input_dict['Fe'] = input_df

input_df = pd.read_csv(FeNi_input_filename)
phox_df = pd.read_csv(FeNi_phox_filename)
# Only use hcp phases
input_df = input_df[input_df['Phase']=='hcp']
input_df = input_df[['Folder','Index','Phase','V','dV','P','dP']]
# Combine input and phox dataframes
input_df = input_df.merge(phox_df,on='Index')
input_dict['FeNi'] = input_df

input_df = pd.read_csv(FeNiSi_input_filename)
phox_df = pd.read_csv(FeNiSi_phox_filename)
# Only use hcp phases
input_df = input_df[input_df['Phase']=='hcp']
input_df = input_df[['Folder','Index','Phase','V','dV','P','dP']]
# Combine input and phox dataframes
input_df = input_df.merge(phox_df,on='Index')
input_dict['FeNiSi'] = input_df


# Calculate Gruneisen parameter at each volume: Fe
##################################################

T = 300

input_df = input_dict['Fe']
V0 = V0_dict['Fe']
dV0 = dV0_dict['Fe']
M = M_dict['Fe']
q = 1
gammavib0 = 2.04
dgammavib0 = 0.1
gammaD0 = 1.74
dgammaD0 = 0.1

results_df = input_df.copy()
results_df = results_df[['Folder','Index','Phase','V','dV','P','dP','KE','dKE',
	'Cvib','dCvib']]

results_df['gamma_vib'], results_df['dgamma_vib'] = calcGruneisen(
	results_df['V'],results_df['dV'],V0,dV0,gammavib0,dgammavib0,q)
results_df['gamma_D'], results_df['dgamma_D'] = calcGruneisen(
	results_df['V'],results_df['dV'],V0,dV0,gammaD0,dgammaD0,q)

# Calculate Pvib at each volume

Cel = calcCel(results_df['V'],V0,T,M)
results_df['Pvib'],results_df['dPvib'] = calcPvib(results_df['V'],results_df['dV'],
	results_df['gamma_vib'], results_df['dgamma_vib'],
	2*results_df['KE'],2*results_df['dKE'],results_df['Cvib'],results_df['dCvib'],
	Cel)
print(results_df)

results_df = results_df.round({'gamma_vib':2,'dgamma_vib':2,'gamma_D':2,'dgamma_D':2,
	'Pvib':2,'dPvib':2})
results_df.to_csv('Fe/Results/GruneisenResults.csv',index=False)



# Calculate Gruneisen parameter at each volume: FeNi
####################################################

input_df = input_dict['FeNi']
V0 = V0_dict['FeNi']
dV0 = dV0_dict['FeNi']
M = M_dict['FeNi']
q = 1
gammavib0 = 2.07
dgammavib0 = 0.1
gammaD0 = 1.69
dgammaD0 = 0.1

results_df = input_df.copy()
results_df = results_df[['Folder','Index','Phase','V','dV','P','dP','KE','dKE',
	'Cvib','dCvib']]

results_df['gamma_vib'], results_df['dgamma_vib'] = calcGruneisen(
	results_df['V'],results_df['dV'],V0,dV0,gammavib0,dgammavib0,q)
results_df['gamma_D'], results_df['dgamma_D'] = calcGruneisen(
	results_df['V'],results_df['dV'],V0,dV0,gammaD0,dgammaD0,q)

# Calculate Pvib at each volume

Cel = calcCel(results_df['V'],V0,T,M)
results_df['Pvib'],results_df['dPvib'] = calcPvib(results_df['V'],results_df['dV'],
	results_df['gamma_vib'], results_df['dgamma_vib'],
	2*results_df['KE'],2*results_df['dKE'],results_df['Cvib'],results_df['dCvib'],
	Cel)
print(results_df)

results_df = results_df.round({'gamma_vib':2,'dgamma_vib':2,'gamma_D':2,'dgamma_D':2,
	'Pvib':2,'dPvib':2})
results_df.to_csv('FeNi/Results/GruneisenResults.csv',index=False)



# Calculate Gruneisen parameter at each volume: FeNiSi
######################################################

input_df = input_dict['FeNiSi']
V0 = V0_dict['FeNiSi']
dV0 = dV0_dict['FeNiSi']
M = M_dict['FeNiSi']
q = 1
gammavib0 = 2.03
dgammavib0 = 0.1
gammaD0 = 1.82
dgammaD0 = 0.1

results_df = input_df.copy()
results_df = results_df[['Folder','Index','Phase','V','dV','P','dP','KE','dKE',
	'Cvib','dCvib']]

results_df['gamma_vib'], results_df['dgamma_vib'] = calcGruneisen(
	results_df['V'],results_df['dV'],V0,dV0,gammavib0,dgammavib0,q)
results_df['gamma_D'], results_df['dgamma_D'] = calcGruneisen(
	results_df['V'],results_df['dV'],V0,dV0,gammaD0,dgammaD0,q)

# Calculate Pvib at each volume

Cel = calcCel(results_df['V'],V0,T,M)
results_df['Pvib'],results_df['dPvib'] = calcPvib(results_df['V'],results_df['dV'],
	results_df['gamma_vib'], results_df['dgamma_vib'],
	2*results_df['KE'],2*results_df['dKE'],results_df['Cvib'],results_df['dCvib'],
	Cel)
print(results_df)

results_df = results_df.round({'gamma_vib':2,'dgamma_vib':2,'gamma_D':2,'dgamma_D':2,
	'Pvib':2,'dPvib':2})
results_df.to_csv('FeNiSi/Results/GruneisenResults.csv',index=False)


