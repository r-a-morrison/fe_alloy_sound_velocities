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


# Input file paths and info
###########################

m_res_isotope = 57
fit_type = 'Constrained_Power'

PDOSpath = dict()
Emin_Debye = dict()
Emax_Debye = dict()
Emin_PowerLaw = dict()
Emax_PowerLaw = dict()

PDOSpath['Fe'] = '../081_vD_Fe_PowerLaw/2009Oct_53GPa/PDOS_Velocity.csv'
PDOSpath['FeNi'] = '../091_vD_FeNi_PowerLaw/2012Apr_DAC11_P5/PDOS_Velocity.csv'
PDOSpath['FeNiSi'] = '../101_vD_FeNiSi_PowerLaw/2015Mar_DAC13_P3/PDOS_Velocity.csv'

# P = 30 GPa (P1)
Emin_Debye['Fe'] = 0      # 3
Emax_Debye['Fe'] = 16.0        # 14.5
Emin_PowerLaw['Fe'] = 0   # 3
Emax_PowerLaw['Fe'] = 22   # 23

# P = 41 GPa (P5)
Emin_Debye['FeNi'] = 0      # 3
Emax_Debye['FeNi'] = 14.5        # 14.5
Emin_PowerLaw['FeNi'] = 0   # 3
Emax_PowerLaw['FeNi'] = 20   # 23

# P = 41 GPa (P3)
Emin_Debye['FeNiSi'] = 0     # 3
Emax_Debye['FeNiSi'] = 12.5  # 13
Emin_PowerLaw['FeNiSi'] = 0    # 3
Emax_PowerLaw['FeNiSi'] = 20    # 22


# PDOSpath['FeNi'] = '../091_vD_FeNi_PowerLaw/2012Apr_DAC11_P9/PDOS_Velocity.csv'
# PDOSpath['FeNiSi'] = '../101_vD_FeNiSi_PowerLaw/2015Mar_DAC13_P6/PDOS_Velocity.csv'

# # P = 83 GPa (P9)
# Emin_Debye['FeNi'] = 3.5
# Emax_Debye['FeNi'] = 17.5
# Emin_PowerLaw['FeNi'] = 3.5
# Emax_PowerLaw['FeNi'] = 27.5

# # P = 86 GPa (P6)
# Emin_Debye['FeNiSi'] = 3.5
# Emax_Debye['FeNiSi'] = 16.5
# Emin_PowerLaw['FeNiSi'] = 3.5
# Emax_PowerLaw['FeNiSi'] = 27


# PDOSpath['FeNi'] = '../091_vD_FeNi_PowerLaw/2012Apr_DAC11_P7/PDOS_Velocity.csv'
# PDOSpath['FeNiSi'] = '../101_vD_FeNiSi_PowerLaw/2015Mar_DAC13_P5/PDOS_Velocity.csv'

# # P = 83 GPa (P7)
# Emin_Debye['FeNi'] = 3.5
# Emax_Debye['FeNi'] = 16.5
# Emin_PowerLaw['FeNi'] = 3.5
# Emax_PowerLaw['FeNi'] = 25

# # P = 86 GPa (P5)
# Emin_Debye['FeNiSi'] = 3.5
# Emax_Debye['FeNiSi'] = 16.5
# Emin_PowerLaw['FeNiSi'] = 3.5
# Emax_PowerLaw['FeNiSi'] = 25



# Functions
###########

def set_E_range(dos_df, Emin, Emax):
	dos_crop_df = dos_df[dos_df['E'] >= Emin]
	dos_crop_df = dos_crop_df[dos_crop_df['E'] <= Emax]

	E_data = dos_crop_df['E']
	vD_data = dos_crop_df['Velocity']
	dvD_data = dos_crop_df['dVelocity']

	return E_data, vD_data, dvD_data


def set_start_params(fit_type,dos_df):
	dos_beg_df = dos_df[dos_df['E'] <= 5]
	A0 = dos_beg_df['Velocity'].mean()
	if fit_type == 'Debye':
		return A0
	elif fit_type == 'Constrained_Power':
		A1 = 10
		return A0, A1
	elif fit_type == 'Unconstrained_Power':
		A1 = 10
		A2 = 2
		return A0, A1, A2
	else:
		print('Error: Invalid model selected.\n')


def Debye_fit(x,A0):
	return A0


def Constrained_Power_fit(x,A0,A1):
	return A0 - (x/A1)**4


# Fit a curve to data
def fit_curve(E_data,vD_data,dvD_data,fit_type):
	
	no_conv = False

	# Define start parameters for fit
	start_params = set_start_params(fit_type,dos_df)

	# Fit curve to data depending on selected model
	if fit_type == 'Debye':
		try:
			par_opt,par_cov = curve_fit(Debye_fit,E_data,vD_data,
				p0=[start_params], sigma = dvD_data)
			vD_fit = [Debye_fit(E,*par_opt) for E in E_data.values]
		except:
			no_conv = True
	elif fit_type == 'Constrained_Power':
		try:
			par_opt,par_cov = curve_fit(Constrained_Power_fit,E_data,vD_data,
				p0=[*start_params], sigma = dvD_data)
			vD_fit = [Constrained_Power_fit(E,*par_opt) for E in E_data.values]
		except:
			no_conv = True
	else:
		print('Error: Model not found.')

	if not no_conv:
		return vD_fit
	else:
		print('No convergence')
		return None


# Import data
#############

dos_df_dict = dict()
dos_df_dict['Fe'] = pd.read_csv(PDOSpath['Fe'])
dos_df_dict['FeNi'] = pd.read_csv(PDOSpath['FeNi'])
dos_df_dict['FeNiSi'] = pd.read_csv(PDOSpath['FeNiSi'])


# Make plots
############

# fig, (ax0) = plt.subplots(nrows = 1, ncols=1, sharex=True, figsize=(8, 6))

# dos_df = dos_df_dict['Fe']

# # Plot Debye fit
# fit_type = 'Debye'
# E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_Debye['Fe'], Emax_Debye['Fe'])
# # Fit curve to data
# vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
# ax0.plot(E_data, vD_fit, '-', color='gray')

# # Plot power law fit
# fit_type = 'Constrained_Power'
# E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_PowerLaw['Fe'], Emax_PowerLaw['Fe'])
# # Fit curve to data
# vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
# ax0.plot(E_data, vD_fit, '--', color='gray')

# # ax0.errorbar(dos_df['E'], dos_df['Velocity'], yerr=dos_df['dVelocity'],
# # 	marker = '.', linestyle='', color='black', mfc='White', elinewidth = 1)
# ax0.plot(dos_df['E'], dos_df['Velocity'],
# 	marker = '.', linestyle='', color='gray', mfc='White')
# ax0.fill_between(dos_df['E'],dos_df['Velocity']-dos_df['dVelocity'],
# 	dos_df['Velocity']+dos_df['dVelocity'],facecolor='gray',alpha=0.3)


# dos_df = dos_df_dict['FeNi']

# # Plot Debye fit
# fit_type = 'Debye'
# E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_Debye['FeNi'], Emax_Debye['FeNi'])
# # Fit curve to data
# vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
# ax0.plot(E_data, vD_fit, '-', color='darkorange')

# # # Plot Debye fit
# # fit_type = 'Debye'
# # E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_Debye['FeNi'], 16)
# # # Fit curve to data
# # vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
# # ax0.plot(E_data, vD_fit, '-', color='darkorange')

# # Plot power law fit
# fit_type = 'Constrained_Power'
# E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_PowerLaw['FeNi'], Emax_PowerLaw['FeNi'])
# # Fit curve to data
# vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
# ax0.plot(E_data, vD_fit, '--', color='darkorange')

# # ax0.errorbar(dos_df['E'], dos_df['Velocity'], yerr=dos_df['dVelocity'],
# # 	marker = '.', linestyle='', color='black', mfc='White', elinewidth = 1)
# ax0.plot(dos_df['E'], dos_df['Velocity'],
# 	marker = '.', linestyle='', color='darkorange', mfc='White')
# ax0.fill_between(dos_df['E'],dos_df['Velocity']-dos_df['dVelocity'],
# 	dos_df['Velocity']+dos_df['dVelocity'],facecolor='darkorange',alpha=0.3)


# dos_df = dos_df_dict['FeNiSi']

# # Plot Debye fit
# fit_type = 'Debye'
# E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_Debye['FeNiSi'], Emax_Debye['FeNiSi'])
# # Fit curve to data
# vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
# ax0.plot(E_data, vD_fit, '-', color='deepskyblue')
# print(vD_fit[0])

# # # Plot Debye fit
# # fit_type = 'Debye'
# # E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_Debye['FeNiSi'], 16)
# # # Fit curve to data
# # vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
# # ax0.plot(E_data, vD_fit, '-', color='deepskyblue')
# # print(vD_fit[0])

# # Plot power law fit
# fit_type = 'Constrained_Power'
# E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_PowerLaw['FeNiSi'], Emax_PowerLaw['FeNiSi'])
# # Fit curve to data
# vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
# ax0.plot(E_data, vD_fit, '--', color='deepskyblue')

# # ax0.errorbar(dos_df['E'], dos_df['Velocity'], yerr=dos_df['dVelocity'],
# # 	marker = '.', linestyle='', color='blue', mfc='White', elinewidth = 1)
# ax0.plot(dos_df['E'], dos_df['Velocity'],
# 	marker = '.', linestyle='', color='deepskyblue', mfc='White')
# ax0.fill_between(dos_df['E'],dos_df['Velocity']-dos_df['dVelocity'],
# 	dos_df['Velocity']+dos_df['dVelocity'],facecolor='deepskyblue',alpha=0.3)

# # if fit != None:
# # 	[E_fit,vD_fit] = fit
# # 	ax0.plot(E_fit, vD_fit, '-', color='red')

# ax0.set_xlabel(r'Energy (meV)', fontsize=14)
# ax0.set_ylabel(r'$(E^2/D(E,V))^{1/3}$ (m/s)', fontsize=14)
# ax0.set_xlim([0,30])
# ax0.set_ylim([2500,6000])

# plt.tight_layout()
# fig = plt.gcf()
# fig.savefig('PDOS_velocity_fit.pdf', format='pdf')
# plt.close()


# Old plot with FeNi and FeNiSi
###############################

fig, (ax0,ax1) = plt.subplots(nrows = 2, ncols=1, sharex=True, figsize=(6, 7))

dos_df = dos_df_dict['FeNi']

# Plot Debye fit
fit_type = 'Debye'
E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_Debye['FeNi'], Emax_Debye['FeNi'])
# Fit curve to data
vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
ax0.plot(E_data, vD_fit*np.ones(len(vD_fit))*10**-3, '-', color='black')

# Plot power law fit
fit_type = 'Constrained_Power'
E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_PowerLaw['FeNi'], Emax_PowerLaw['FeNi'])
# Fit curve to data
vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
ax0.plot(E_data, vD_fit*np.ones(len(vD_fit))*10**-3, '--', color='black')

ax0.errorbar(dos_df['E'], dos_df['Velocity']/1E3, yerr=dos_df['dVelocity']/1E3,
	marker = '.', linestyle='', color='darkorange', mfc='White', elinewidth = 1)
# ax0.plot(dos_df['E'], dos_df['Velocity'],
# 	marker = '.', linestyle='', color='darkorange', mfc='White')
# ax0.fill_between(dos_df['E'],dos_df['Velocity']-dos_df['dVelocity'],
# 	dos_df['Velocity']+dos_df['dVelocity'],facecolor='darkorange',alpha=0.3)

ax0.set_ylim([2.7,4.7])
ax0.set_ylabel(r'$v(E)$ (km/s)', fontsize=14)


dos_df = dos_df_dict['FeNiSi']

# Plot Debye fit
fit_type = 'Debye'
E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_Debye['FeNiSi'], Emax_Debye['FeNiSi'])
# Fit curve to data
vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
ax1.plot(E_data, vD_fit*np.ones(len(vD_fit))*10**-3, '-', color='black')
print(vD_fit[0])

# Plot power law fit
fit_type = 'Constrained_Power'
E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_PowerLaw['FeNiSi'], Emax_PowerLaw['FeNiSi'])
# Fit curve to data
vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
ax1.plot(E_data, vD_fit*np.ones(len(vD_fit))*10**-3, '--', color='black')

ax1.errorbar(dos_df['E'], dos_df['Velocity']/1E3, yerr=dos_df['dVelocity']/1E3,
	marker = '.', linestyle='', color='deepskyblue', mfc='White', elinewidth = 1)
# ax1.plot(dos_df['E'], dos_df['Velocity'],
# 	marker = '.', linestyle='', color='deepskyblue', mfc='White')
# ax1.fill_between(dos_df['E'],dos_df['Velocity']-dos_df['dVelocity'],
# 	dos_df['Velocity']+dos_df['dVelocity'],facecolor='deepskyblue',alpha=0.3)

ax1.set_xlabel(r'Energy (meV)', fontsize=14)
ax1.set_ylabel(r'$v(E)$ (km/s)', fontsize=14)
ax1.set_xlim([0,28])
ax1.set_ylim([3,4.7])

plt.tight_layout()
fig = plt.gcf()
fig.savefig('PDOS_velocity_fit.pdf', format='pdf')
plt.close()


# Old plot with FeNi and FeNiSi
###############################

fig, (ax0,ax1,ax2) = plt.subplots(nrows = 3, ncols=1, sharex=True, figsize=(7, 10.4))

dos_df = dos_df_dict['Fe']

# Plot Debye fit
fit_type = 'Debye'
E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_Debye['Fe'], Emax_Debye['Fe'])
# Fit curve to data
vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
ax0.plot(E_data, vD_fit*np.ones(len(vD_fit))*10**-3, '-', color='black')

# Plot power law fit
fit_type = 'Constrained_Power'
E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_PowerLaw['Fe'], Emax_PowerLaw['Fe'])
# Fit curve to data
vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
ax0.plot(E_data, vD_fit*np.ones(len(vD_fit))*10**-3, '--', color='black')

ax0.errorbar(dos_df['E'], dos_df['Velocity']/1E3, yerr=dos_df['dVelocity']/1E3,
	marker = '.', linestyle='', color='gray', mfc='White', elinewidth = 1)

ax0.set_ylim([3.1,5.1])
ax0.set_ylabel(r'$v(E)$ (km/s)', fontsize=18)


dos_df = dos_df_dict['FeNi']

# Plot Debye fit
fit_type = 'Debye'
E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_Debye['FeNi'], Emax_Debye['FeNi'])
# Fit curve to data
vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
ax1.plot(E_data, vD_fit*np.ones(len(vD_fit))*10**-3, '-', color='black')

# Plot power law fit
fit_type = 'Constrained_Power'
E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_PowerLaw['FeNi'], Emax_PowerLaw['FeNi'])
# Fit curve to data
vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
ax1.plot(E_data, vD_fit*np.ones(len(vD_fit))*10**-3, '--', color='black')

ax1.errorbar(dos_df['E'], dos_df['Velocity']/1E3, yerr=dos_df['dVelocity']/1E3,
	marker = '.', linestyle='', color='darkorange', mfc='White', elinewidth = 1)
# ax1.plot(dos_df['E'], dos_df['Velocity'],
# 	marker = '.', linestyle='', color='darkorange', mfc='White')
# ax1.fill_between(dos_df['E'],dos_df['Velocity']-dos_df['dVelocity'],
# 	dos_df['Velocity']+dos_df['dVelocity'],facecolor='darkorange',alpha=0.3)

ax1.set_ylim([2.7,4.7])
ax1.set_ylabel(r'$v(E)$ (km/s)', fontsize=18)


dos_df = dos_df_dict['FeNiSi']

# Plot Debye fit
fit_type = 'Debye'
E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_Debye['FeNiSi'], Emax_Debye['FeNiSi'])
# Fit curve to data
vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
ax2.plot(E_data, vD_fit*np.ones(len(vD_fit))*10**-3, '-', color='black')

# Plot power law fit
fit_type = 'Constrained_Power'
E_data, vD_data, dvD_data = set_E_range(dos_df, Emin_PowerLaw['FeNiSi'], Emax_PowerLaw['FeNiSi'])
# Fit curve to data
vD_fit = fit_curve(E_data,vD_data,dvD_data,fit_type)
ax2.plot(E_data, vD_fit*np.ones(len(vD_fit))*10**-3, '--', color='black')

ax2.errorbar(dos_df['E'], dos_df['Velocity']/1E3, yerr=dos_df['dVelocity']/1E3,
	marker = '.', linestyle='', color='deepskyblue', mfc='White', elinewidth = 1)
# ax2.plot(dos_df['E'], dos_df['Velocity'],
# 	marker = '.', linestyle='', color='deepskyblue', mfc='White')
# ax2.fill_between(dos_df['E'],dos_df['Velocity']-dos_df['dVelocity'],
# 	dos_df['Velocity']+dos_df['dVelocity'],facecolor='deepskyblue',alpha=0.3)

ax2.set_xlabel(r'Energy (meV)', fontsize=18)
ax2.set_ylabel(r'$v(E)$ (km/s)', fontsize=18)
ax2.set_xlim([0,28])
ax2.set_ylim([3,4.7])

plt.tight_layout()
fig = plt.gcf()
fig.savefig('PDOS_velocity_fit2.pdf', format='pdf')
plt.close()