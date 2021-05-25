# Substitute for psvl

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
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
from scipy.interpolate import spline

# Seaborn, useful for graphics
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
fit_type = 'Debye'

phoxpath = '../060_phox_FeNi_man/'
K_bcc_path = '../015_BulkModulusDet/Results/bccFeNi.csv'
K_hcp_path = '../015_BulkModulusDet/Results/hcpFeNi.csv'

start_fit_path = 'StartFitParams.csv'
choices_file_path = 'ChoiceFitParams.csv'


# Functions
###########

def transform_DOS(dos_df,m_res_isotope,rho,drho):
	# f(E) = [(m/(2*pi^2*hbar^3*rho))*E^2/PDOS]^(1/3)
	# units:      (kg/atom)        (meV*eV/1000meV)^2
	#        ------------------- * ------------------
	#        (eV*s)^3 * (kg/m^3)       1/eV
	m_res_isotope_SI = m_res_isotope/(constants.N_A*1000) # g/mol to kg/atom
	hbar = 6.58212E-16  # eV s
	rho_SI = rho*1000   # g/cm^3 to kg/m^3
	drho_SI = drho*1000 # g/cm^3 to kg/m^3
	coeff = m_res_isotope_SI/(2*np.pi**2*hbar**3)
	# PSVL has an extra 3 factor that I don't understand.
	# Including this factor gives identical results to psvl.
	coeff = 3*coeff
	# energy in meV, but we need eV for calculation
	# DOS in 1/eV

	# If DOS = 0, divide by zero warning will appear later. We can't use this first
	# data point anyway.
	dos_df = dos_df[dos_df['DOS'] != 0.0].copy()
		
	# Transform to velocity units
	# Pandas can't handle the cubic root of a negative number or of zero.
	# This corrects for that
	cubic_root = lambda x: 0 if x==0 else np.sign(x) * np.power(abs(x), 1./3)
	dos_df['Velocity'] = [cubic_root((coeff/rho_SI)*E**2/DOS)
		for E, DOS in zip(dos_df['E']/1000,dos_df['DOS'])]

	# Calculate uncertainty in velocity units
	# df(E) = (1/3)*(E^2/DOS^4)^(1/3)*dDOS
	dos_df['dVelocity'] = [np.sqrt( ((1/3)*cubic_root((coeff/rho_SI)*E**2/(DOS**4))*dDOS)**2 + 
		((1/3)*cubic_root((coeff/rho_SI**4)*E**2/(DOS))*drho_SI)**2)
		for E, DOS, dDOS in zip(dos_df['E']/1000,dos_df['DOS'],dos_df['dDOS'])]
		
	return(dos_df) # Units of m/s


# Create a plot of transformed PDOS
def plot_transform_DOS(folder,index,dos_df,fit):
	fig, (ax0) = plt.subplots(nrows = 1, ncols=1, sharex=True, figsize=(8, 6))

	ax0.plot(dos_df['E'], dos_df['Velocity'],
		marker = '', linestyle='', color='black', mfc='White')
	ax0.errorbar(dos_df['E'], dos_df['Velocity'], yerr=dos_df['dVelocity'],
		marker = '', linestyle='', color='black', mfc='White', elinewidth = 1)

	if fit != None:
		[E_fit,V_fit] = fit
		ax0.plot(E_fit, V_fit, '-', color='red')
	
	ax0.set_xlabel(r'Energy (meV)', fontsize=14)
	ax0.set_ylabel(r'$(E^2/D(E,V))^{1/3}$ (m/s)', fontsize=14)
	ax0.set_xlim([0,40])
	ax0.set_ylim([2000,8000])
	
	plt.tight_layout()
	fig = plt.gcf()
	fig.savefig(folder+'/PDOS_Velocity_'+index+'.pdf', format='pdf')
	plt.close()


def set_E_range(dos_df, Emin, Emax):
	dos_crop_df = dos_df[dos_df['E'] >= Emin]
	dos_crop_df = dos_crop_df[dos_crop_df['E'] <= Emax]

	E_data = dos_crop_df['E']
	V_data = dos_crop_df['Velocity']
	dV_data = dos_crop_df['dVelocity']

	return E_data, V_data, dV_data


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
	return A0*(1 - (x/A1)**2)

def Unconstrained_Power_fit(x,A0,A1,A2):
	return A0*(1 - (x/A1)**A2)


# Fit a curve to data
def fit_curve(E_data,V_data,dV_data,fit_type):
	
	no_conv = False

	# Define start parameters for fit
	start_params = set_start_params(fit_type,dos_df)

	# Fit curve to data depending on selected model
	if fit_type == 'Debye':
		try:
			par_opt,par_cov = curve_fit(Debye_fit,E_data,V_data,
				p0=start_params, sigma = dV_data)
			V_fit = [Debye_fit(E,*par_opt) for E in E_data.values]
		except:
			no_conv = True
	elif fit_type == 'Constrained_Power':
		try:
			par_opt,par_cov = curve_fit(Constrained_Power_fit,E_data,V_data,
				p0=[*start_params], sigma = dV_data)
			V_fit = [Constrained_Power_fit(E,*par_opt) for E in E_data.values]
		except:
			no_conv = True
	elif fit_type == 'Unconstrained_Power':
		try:
			par_opt,par_cov = curve_fit(Unconstrained_Power_fit,E_data,V_data,
				p0=[*start_params], sigma = dV_data)
			V_fit = [Unconstrained_Power_fit(E,*par_opt) for E in E_data.values]
		except:
			no_conv = True
	else:
		print('Error: Model not found.')

	if not no_conv:
		# Get number of data points
		N = len(E_data)

		# Calculate reduced chi^2
		chi2 = calc_chi2(V_data,V_fit,dV_data,N,M)
		
		# Calculate AICc
		AICc = calc_AICc(chi2, M, N)

		# Calclulate "probability" from AICc
		prob = calc_pdfAICc(AICc)

		return par_opt, par_cov, V_fit, N, chi2, AICc, prob, no_conv
	else:
		return None, None, None, None, None, None, None, no_conv


# Calculate reduced chi^2
def calc_chi2(data,fit,sigma,N,M):
	chi2 = np.sum((data-fit)**2/sigma**2)/(N-M)
	return chi2


# Calculate AICc - corrected Aikake information criterion
# Useful for comparing fits with different numbers of model parameters or data points
# Correction accounts for small number of data points
# Most other information criteria assume N is large, which wouldn't work here
# M = number of model parameters
# N = number of data points
def calc_AICc(misfit, M, N):
	AIC = misfit + 2*M
	correction = 2*M*(M+1)/(N-M-1) # Corrects for cases with small N
	AICc = AIC + correction
	return AICc


# AICc is analogous to a log posterior. Transform back to posterior space to create
# a pdf
def calc_pdfAICc(AICc):
	return np.exp(-(1/2)*AICc)


def V_D_from_fit(par_opt, par_cov):
	V_D = par_opt[0]
	dV_D = np.sqrt(par_cov[0,0])
	return V_D, dV_D


# Gaussian function for fitting
def Gaussian(x,amplitude,mean,sigma):
	return amplitude*np.exp(-(x-mean)**2/(2*sigma**2))


# Bin data by V_D to get pdf of V_D
def bin_V_D(psvl_all_df,bin_width):
	min_V_D = psvl_all_df['V_D'].min()
	max_V_D = psvl_all_df['V_D'].max()
	V_D_start = min_V_D - (min_V_D%bin_width) - bin_width
	V_D_end = max_V_D - (max_V_D%bin_width) + 2*bin_width
	V_D_bins = np.arange(V_D_start,V_D_end,bin_width)
	psvl_all_df['Bin'] = pd.cut(psvl_all_df['V_D'],V_D_bins)
	psvl_all_df['Bin Midpoint'] = [psvl_all_df['Bin'].iloc[i].mid for i in psvl_all_df.index]
	grouped = psvl_all_df.groupby(['Bin Midpoint'])
	binned_df = grouped.sum()[['Prob']].reset_index()
	return binned_df


# Pivot data for plotting
def df_2_pivotdata(psvl_all_df,param):
	psvl_all_df = psvl_all_df.drop_duplicates(subset=['Emin', 'Emax'], keep='first')
	# NaN automatically assigned when missing
	pivotplot = psvl_all_df.pivot('Emin', 'Emax', param)
	Emin_array = pivotplot.index.values
	Emax_array = pivotplot.columns.values
	param_pivot = pivotplot.values
	return Emin_array, Emax_array, param_pivot


def plot_AICc_analysis(psvl_all_df,binned_df,Gauss_opt,folder,index,fit_type,
	E_choice,title):
	# Pivot data so we can make a plot
	Emin_array, Emax_array, chi2_pivot = df_2_pivotdata(psvl_all_df,'chi2')
	Emin_array, Emax_array, prob_pivot = df_2_pivotdata(psvl_all_df,'Prob')
	Emin_array, Emax_array, V_D_pivot  = df_2_pivotdata(psvl_all_df,'V_D')
	Emin_array, Emax_array, dV_D_pivot = df_2_pivotdata(psvl_all_df,'dV_D')

	fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows = 4, ncols=1, figsize=(16,16))
	cmap_choice = 'magma'
	big_dot_size = 10

	if E_choice != None:
		Emin_choice = E_choice[0]
		Emax_choice = E_choice[1]

	pdflim_min = max([min(binned_df['Bin Midpoint']),2000])
	pdflim_max = min([max(binned_df['Bin Midpoint']),10000])

	# Plot chi2 data
	cf0 = ax0.pcolormesh(Emax_array,Emin_array,chi2_pivot,edgecolor='face',
		cmap=cmap_choice+'_r')
	cbar0 = plt.colorbar(cf0, ax = ax0, aspect=10)
	cbar0.ax.set_ylabel(r'$\chi^2$',fontsize=18)
	if E_choice != None:
		ax0.plot(Emax_choice,Emin_choice,'o',color='white',ms=big_dot_size)
	ax0.set_xlim([min(Emax_array),max(Emax_array)])
	ax0.set_ylim([min(Emin_array),max(Emin_array)])
	ax0.tick_params(direction='out')
	ax0.set_ylabel('Emin',fontsize=18)
	ax0.set_title(title,fontsize=18)

	cf1 = ax1.pcolormesh(Emax_array,Emin_array,prob_pivot,edgecolor='face',cmap=cmap_choice)
	cbar1 = plt.colorbar(cf1, ax = ax1, aspect=10)
	cbar1.ax.set_ylabel(r'Posterior (unnormalized)',fontsize=18)
	if E_choice != None:
		ax1.plot(Emax_choice,Emin_choice,'o',color='white',ms=big_dot_size)
	ax1.set_xlim([min(Emax_array),max(Emax_array)])
	ax1.set_ylim([min(Emin_array),max(Emin_array)])
	ax1.tick_params(direction='out')
	ax1.set_ylabel('Emin',fontsize=18)

	# Plot V_D data
	cf2 = ax2.pcolormesh(Emax_array,Emin_array,V_D_pivot,edgecolor='face',cmap=cmap_choice,
		vmax=pdflim_max)
	cbar2 = plt.colorbar(cf2, ax = ax2, aspect=10)
	cbar2.ax.set_ylabel(r'$v_D$',fontsize=20)
	if E_choice != None:
		ax2.plot(Emax_choice,Emin_choice,'o',color='white',ms=big_dot_size)
	ax2.set_xlim([min(Emax_array),max(Emax_array)])
	ax2.set_ylim([min(Emin_array),max(Emin_array)])
	ax2.tick_params(direction='out')
	ax2.set_xlabel('Emax',fontsize=18)
	ax2.set_ylabel('Emin',fontsize=18)

	# Plot V_D pdf data
	ax3.plot(binned_df['Bin Midpoint'],binned_df['Prob'],'.')
	ax3.plot(binned_df['Bin Midpoint'],Gaussian(binned_df['Bin Midpoint'],*Gauss_opt),
		color='red')
	ax3.set_xlim([pdflim_min,pdflim_max])
	ax3.tick_params(direction='out')
	ax3.set_xlabel(r'$v_D$',fontsize=18)
	ax3.set_ylabel(r'pdf of $v_D$',fontsize=18)

	# Fine-tune figure; make subplots close to each other and hide x ticks for
	# all but bottom plot.
	fig.subplots_adjust(hspace=0.05)
	plt.setp([ax0.get_xticklabels() for a in fig.axes[:1]], visible=False);
	plt.setp([ax1.get_xticklabels() for a in fig.axes[:1]], visible=False);

	plt.tight_layout()
	fig.savefig(folder+'/AICc_plot_'+index+'_'+fit_type+'.pdf', format='pdf')
	plt.close()


# Load datasets
###############

K_bcc_df = pd.read_csv(K_bcc_path, engine='python')
K_hcp_df = pd.read_csv(K_hcp_path, engine='python')


# Get results from XRD and phox analysis
############################################

# phox directory structure needs to be the same as psvl directory structure

# Find the filepath of all .res NRIXS files in phox directories
respath_list = [filepath for filepath in glob.glob(phoxpath+'*/*.res')]

# Prep lists, dictionaries, and df to store input data in
folder_list = []
index_dict = dict()
dospath_dict = dict()
in_psvlpath_dict = dict()

input_df = pd.DataFrame(columns = ['Folder','Index','Phase','P','dP', #'V','dV',
									'rho','drho','K_S','dK_S'])

# Collect folders, indices, paths, and input values
for respath in respath_list:

	# Determine filepaths for dos and in_psvl
	folder = re.findall('([A-Za-z0-9_]+)/[A-Za-z0-9_]+.res',respath)[0]
	index = re.findall('/([A-Za-z0-9_]+).res',respath)[0]
	dospath = phoxpath+folder+'/Output/'+index+'_dos.dat'
	in_psvlpath = folder+'/in_psvl'

	# Store to use paths later
	folder_list.append(folder)
	index_dict[folder] = index
	dospath_dict[folder] = dospath
	in_psvlpath_dict[folder] = in_psvlpath

	# Determine if measurement is in bcc or hcp phase (stored in XRD and bulk modulus data)
	if index in K_bcc_df['Index'].values:
		phase = 'bcc'
		missingXRD = False
	elif index in K_hcp_df['Index'].values:
		phase = 'hcp'
		missingXRD = False
	else:
		print('Warning: Phase not found for '+index+
			'. The in_psvl file will not be updated for '+index+'.')
		missingXRD = True
	
	# Only procede for the measurements that are good.
	# XRD analysis was not done for some of the bad NRIXS data.
	if not missingXRD:

		# We'll assume K_S = K_T for bcc phases
		# For hcp phases, use calculated K_S and dK_S.

		# Get values for enriched density, KT (bcc), KS (hcp) and uncertainties
		# Bulk modulus code saves enriched value of density from XRD calcs
		if phase == 'bcc':
			# First select the row that matches the filename.
			# Then pick the value in the last (only) row in the 'rho' column
			P = K_bcc_df[K_bcc_df['Index']==index].iloc[-1]['P']
			dP = K_bcc_df[K_bcc_df['Index']==index].iloc[-1]['dP']
			# V = K_bcc_df[K_bcc_df['Index']==index].iloc[-1]['V']
			# dV = K_bcc_df[K_bcc_df['Index']==index].iloc[-1]['dV']
			rho = K_bcc_df[K_bcc_df['Index']==index].iloc[-1]['rho']
			drho = K_bcc_df[K_bcc_df['Index']==index].iloc[-1]['drho']
			K = K_bcc_df[K_bcc_df['Index']==index].iloc[-1]['KT']
			dK = K_bcc_df[K_bcc_df['Index']==index].iloc[-1]['dKT']
		if phase == 'hcp':
			P = K_hcp_df[K_hcp_df['Index']==index].iloc[-1]['P']
			dP = K_hcp_df[K_hcp_df['Index']==index].iloc[-1]['dP']
			# V = K_hcp_df[K_hcp_df['Index']==index].iloc[-1]['V']
			# dV = K_hcp_df[K_hcp_df['Index']==index].iloc[-1]['dV']
			rho = K_hcp_df[K_hcp_df['Index']==index].iloc[-1]['rho']
			drho = K_hcp_df[K_hcp_df['Index']==index].iloc[-1]['drho']
			K = K_hcp_df[K_hcp_df['Index']==index].iloc[-1]['KS']
			dK = K_hcp_df[K_hcp_df['Index']==index].iloc[-1]['dKS']

		# Store in a dataframe to use later
		input_df = input_df.append(pd.DataFrame([[folder,index,phase,P,dP, #V,dV,
			rho,drho,K,dK]], columns = input_df.columns))

if not os.path.isdir('Results'):
	os.mkdir('Results')

# Save input values for reference
input_df.to_csv('Results/input_values.csv',index=False)

# Update folder_list to exclude files not found in XRD
folder_list = input_df['Folder'].values


# Create directories if they don't exist
########################################

for folder in folder_list:
	if not os.path.isdir(folder):
		os.mkdir(folder)


# Import PDOS, transform, and plot
##################################

for folder in folder_list:
	index = index_dict[folder]
	rho = input_df[input_df['Folder']==folder]['rho'].iloc[0]
	drho = input_df[input_df['Folder']==folder]['drho'].iloc[0]
	
	# Import PDOS
	dos_df = pd.read_csv(dospath_dict[folder], sep='\s+', comment='@', header=None,
		names = ['E','DOS','dDOS'])

	dos_df = transform_DOS(dos_df,m_res_isotope,rho,drho)

	# Create a plot of transformed PDOS
	plot_transform_DOS(folder,index,dos_df,None)

	# Save results to csv
	dos_df.to_csv(folder+'/PDOS_Velocity.csv',index=False)


# # Fit model to transformed PDOS for a range of energies
# #######################################################

# # Import csv file with fit ranges
# fit_range_df = pd.read_csv(start_fit_path, engine='python')

# for folder in folder_list:
# 	index = index_dict[folder]

# 	print('Starting fit for '+index)
# 	start_fit_time = time.time()

# 	# Import transformed PDOS
# 	dos_df = pd.read_csv(folder+'/PDOS_Velocity.csv')
# 	dos_df = dos_df.dropna()

# 	# Get fit ranges for this pressure point
# 	Emin_start = fit_range_df[fit_range_df['Folder']==folder].iloc[-1]['Emin start']
# 	Emin_end = fit_range_df[fit_range_df['Folder']==folder].iloc[-1]['Emin end']
# 	Emax_start = fit_range_df[fit_range_df['Folder']==folder].iloc[-1]['Emax start']
# 	Emax_end = fit_range_df[fit_range_df['Folder']==folder].iloc[-1]['Emax end']
# 	# If exact stepsize is chosen, occasionally points may be missed.
# 	E_stepsize = fit_range_df[fit_range_df['Folder']==folder].iloc[-1]['Step size']

# 	# Get an array of Emin values
# 	Emin_steps = (Emin_end - Emin_start)/E_stepsize + 1
# 	Emin_array = np.linspace(Emin_start,Emin_end,Emin_steps)

# 	vel_results_df = pd.DataFrame(columns=['Emin', 'Emax', 'M', 'N', 'chi2',
# 		'AICc', 'Prob', 'V_D', 'dV_D'])

# 	M = np.array(set_start_params(fit_type,dos_df)).size
	
# 	for Emin in Emin_array[:-(M+1)]:
# 		print('Emin: '+str(Emin))
# 		# N-M-1 must be greater than zero or AICc will fail
# 		Emax_safe_start = Emin + E_stepsize*(M+2)
# 		Emax_steps = (Emax_end - Emax_safe_start)/E_stepsize + 1
# 		Emax_array = np.linspace(Emax_safe_start,Emax_end,Emax_steps)
		
# 		for Emax in Emax_array:
# 			E_data, V_data, dV_data = set_E_range(dos_df, Emin, Emax)

# 			# Fit curve to data
# 			par_opt, par_cov, V_fit, N, chi2, AICc, prob, no_conv = fit_curve(E_data,V_data,
# 																		dV_data,fit_type)
# 			if not no_conv:
# 				# Get Debye velocity and uncertainty from this fit
# 				V_D, dV_D = V_D_from_fit(par_opt, par_cov)

# 				temp_df = pd.DataFrame([[Emin, Emax, M, N, chi2,
# 										AICc, prob, V_D, dV_D]],
# 										 columns = ['Emin', 'Emax', 'M', 'N', 'chi2',
# 													'AICc', 'Prob', 'V_D', 'dV_D'])
# 				vel_results_df = vel_results_df.append(temp_df)

# 	# Save the results
# 	vel_results_df.to_csv(folder+'/fit_search_results_'+fit_type+'.csv',index=False)

# 	# How long did looping over Emin and Emax take?
# 	print('Fit time for '+index+': '+str(int(time.time() - start_fit_time))+' s')

# 	# # Create a plot of transformed PDOS with fit curve
# 	# plot_transform_DOS(folder,index,dos_df,[E_data,V_fit])


# Collect fit ranges selected by previous analyses for comparison
#############################################################################

# Import our choices if available
if os.path.exists(choices_file_path):
	choices_df = pd.read_csv(choices_file_path)
	choices_avail = True
else:
	choices_avail = False

E_choice = dict()

for folder in folder_list:
	index = index_dict[folder]

	if choices_avail:
		choices_Emin = choices_df[choices_df['Folder']==folder].iloc[-1]['Emin']
		choices_Emax = choices_df[choices_df['Folder']==folder].iloc[-1]['Emax']
		E_choice[folder] = [choices_Emin, choices_Emax]
	else:
		E_choice[folder] = None


# Make plots of fits, calculate V_D pdf, best V_D, and uncertainty
##################################################################

V_D_est_df = pd.DataFrame(columns = ['Folder','Index','V_D','dV_D','Prob Amp'])

for folder in folder_list:
	index = index_dict[folder]
	
	vel_results_df = pd.read_csv(folder+'/fit_search_results_'+fit_type+'.csv')

	# Calculate pdf of V_D using AICc
	print(index+': Calculating Debye velocity pdf')
	bin_width = 1
	binned_df = bin_V_D(vel_results_df,bin_width)

	# Save V_D pdf results
	binned_df['Bin Midpoint'] = binned_df['Bin Midpoint'].astype(float).round(1)
	binned_df['Prob'] = binned_df['Prob'].astype(float).round(4)
	binned_df.to_csv(folder+'/pdf_results_'+fit_type+'.csv',index=False)

	# Fit Gaussian to V_D pdf
	print(index+': Fitting Gaussian to Debye velocity pdf')
	xlim_min = max([min(binned_df['Bin Midpoint']),2000])
	xlim_max = min([max(binned_df['Bin Midpoint']),10000])
	V_D_guess = xlim_min + (xlim_max - xlim_min)/2
	Gauss_opt,Gauss_cov = curve_fit(Gaussian,binned_df['Bin Midpoint'],binned_df['Prob'],
		p0=[1,V_D_guess,50])
	# Save results of V_D to df
	V_D_norm = 1/Gauss_opt[0]
	prob_amp = Gauss_opt[0]
	V_D = Gauss_opt[1]
	dV_D = Gauss_opt[2]
	temp_df = pd.DataFrame([[folder,index,V_D,dV_D,prob_amp]],
		columns = ['Folder','Index','V_D','dV_D','Prob Amp'])
	V_D_est_df = V_D_est_df.append(temp_df)

	# Plot AICc, V_D, pdf of V_D
	print(index+': Plotting AICc and Debye velocity pdf')
	P = input_df[input_df['Folder']==folder]['P'].iloc[0]
	plot_title = index+'    '+fit_type+' Fit    P = '+str(P)+' GPa'
	plot_AICc_analysis(vel_results_df,binned_df,Gauss_opt,folder,index,fit_type,
		E_choice[folder],plot_title)

# Save V_D pdf fit results
V_D_est_df['V_D'] = V_D_est_df['V_D'].astype(int)
V_D_est_df['dV_D'] = V_D_est_df['dV_D'].astype(int)
V_D_est_df.to_csv('Results/V_D_results_'+fit_type+'.csv',index=False)

print('Total time elapsed: '+str(int(time.time() - start_time))+' s')


# Create a plot using an example best fit
#########################################

for folder in folder_list:
	index = index_dict[folder]
	
	# Import transformed PDOS
	dos_df = pd.read_csv(folder+'/PDOS_Velocity.csv')
	dos_df = dos_df.dropna()

	Emin_choice = E_choice[folder][0]
	Emax_choice = E_choice[folder][1]

	M = np.array(set_start_params(fit_type,dos_df)).size

	E_data, V_data, dV_data = set_E_range(dos_df, Emin_choice, Emax_choice)

	# Fit curve to data
	par_opt, par_cov, V_fit, N, chi2, AICc, prob, no_conv = fit_curve(E_data,V_data,
																dV_data,fit_type)

	# Create a plot of transformed PDOS
	plot_transform_DOS(folder,index,dos_df,[E_data, V_fit])
