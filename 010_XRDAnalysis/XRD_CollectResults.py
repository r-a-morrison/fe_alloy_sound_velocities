import warnings
warnings.simplefilter(action = 'ignore', category = UserWarning)

# Front matter
import os
import glob
import re
import pandas as pd
import numpy as np
import scipy.constants as constants
import sympy as sp
from sympy import Matrix, Symbol
from sympy.utilities.lambdify import lambdify
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec

# Seaborn, useful for graphics
import seaborn as sns

matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

rc = {'lines.linewidth': 1, 
      'axes.labelsize': 20, 
      'axes.titlesize': 20,
      'legend.fontsize': 26,
      'xtick.direction': u'in',
      'ytick.direction': u'in'}
sns.set_style('ticks', rc=rc)

# Functions
def calc_V_bcc(a):
	return a**3

def calc_V_hcp(a,c):
	return (np.sqrt(3)/2)*a**2*c

def calc_dV_bcc(a,da):
	return 3*a**2*da

def calc_dV_hcp(a,c,da,dc):
	return np.sqrt( (np.sqrt(3)*a*c*da)**2 + ((np.sqrt(3)/2)*a**2*dc)**2 )

# Numeric Vinet EOS, used for everything except calculating dP
def VinetEOS(V,V0,K0,Kprime0):
    A = V/V0
    P = 3*K0*A**(-2/3) * (1-A**(1/3)) * np.exp((3/2)*(Kprime0-1)*(1-A**(1/3)))
    return P

# Symbolic Vinet EOS, needed to calculate dP
def VinetEOS_sym(V,V0,K0,Kprime0):
    A = V/V0
    P = 3*K0*A**(-2/3) * (1-A**(1/3)) * sp.exp((3/2)*(Kprime0-1)*(1-A**(1/3)))
    return P

# Create a covariance matrix from EOS_df with V0, K0, and K0prime; used to get dP
def getCov3(EOS_df, phase):
	dV0 = np.float(EOS_df[EOS_df['Phase'] == phase]['dV0'])
	dK0 = np.float(EOS_df[EOS_df['Phase'] == phase]['dK0'])
	dKprime0 = np.float(EOS_df[EOS_df['Phase'] == phase]['dKprime0'])
	V0K0_corr = np.float(EOS_df[EOS_df['Phase'] == phase]['V0K0 corr'])
	V0Kprime0_corr = np.float(EOS_df[EOS_df['Phase'] == phase]['V0Kprime0 corr'])
	K0Kprime0_corr = np.float(EOS_df[EOS_df['Phase'] == phase]['K0Kprime0 corr'])
	corr_matrix = np.eye(3)
	corr_matrix[0,1] = V0K0_corr
	corr_matrix[1,0] = V0K0_corr
	corr_matrix[0,2] = V0Kprime0_corr
	corr_matrix[2,0] = V0Kprime0_corr
	corr_matrix[1,2] = K0Kprime0_corr
	corr_matrix[2,1] = K0Kprime0_corr
	sigmas = np.array([[dV0,dK0,dKprime0]])
	cov = (sigmas.T@sigmas)*corr_matrix
	return cov

# Create a covariance matrix with V, V0, K0, and K0prime; used to get dP
def getVinetCov(dV, EOS_df, phase):
	cov3 = getCov3(EOS_df, phase)
	cov = np.eye(4)
	cov[1:4,1:4] = cov3
	cov[0,0] = dV**2
	return cov

def calc_dP_VinetEOS(V, dV, EOS_df, phase):
	# Create function for Jacobian of Vinet EOS
	a,b,c,d = Symbol('a'),Symbol('b'),Symbol('c'),Symbol('d') # Symbolic variables V, V0, K0, K'0
	Vinet_matrix = Matrix([VinetEOS_sym(a,b,c,d)])            # Create a symbolic Vinet EOS matrix
	param_matrix = Matrix([a,b,c,d])                          # Create a matrix of symbolic variables
	# Symbolically take the Jacobian of the Vinet EOS and turn into a column matrix
	J_sym = Vinet_matrix.jacobian(param_matrix).T
	# Create a numpy function for the above expression
	# (easier to work with numerically)
	J_Vinet = lambdify((a,b,c,d), J_sym, 'numpy')
	J = J_Vinet(V,*getEOSparams(EOS_df, phase)) # Calculate Jacobian
	cov = getVinetCov(dV, EOS_df, phase) # Calculate covariance matrix
	dP = (J.T@cov@J).item() # Calculate uncertainty and convert to a scalar
	return dP

def getEOSparams(EOS_df, phase):
	V0 = np.float(EOS_df[EOS_df['Phase'] == phase]['V0'])
	K0 = np.float(EOS_df[EOS_df['Phase'] == phase]['K0'])
	Kprime0 = np.float(EOS_df[EOS_df['Phase'] == phase]['Kprime0'])
	return V0, K0, Kprime0

def calc_rho(V,dV,M):
	# Convert from cubic angstroms to cm^3/mol
	V_ccpermol = (V/2)*constants.N_A/(10**24)
	rho = M/V_ccpermol
	drho = (M*2*10**24/constants.N_A)*(dV/(V**2))
	return rho, drho

# Import EOS information
EOS_df = pd.read_csv('FeAlloyEOS.csv')

# Find the filepath of all .xy XRD pattern files
patternfilepath_list = [filepath for filepath in glob.glob('*/*/*.xy')]

allresults_df = pd.DataFrame()
bccFe_df = pd.DataFrame()
bccFeNi_df = pd.DataFrame()
bccFeNiSi_df = pd.DataFrame()
hcpFeNi_df = pd.DataFrame()
hcpFeNiSi_df = pd.DataFrame()

for patternfilepath in patternfilepath_list:
	filepath = re.findall('([A-Za-z0-9/_.]+/)[A-Za-z0-9_]+.xy',patternfilepath)[0]
	filename = re.findall('/([A-Za-z0-9_]+).xy',patternfilepath)[0]

	patternresults_df = pd.read_csv(filepath+filename+'_results.csv')
	paramresults_df = pd.read_csv(filepath+filename+'_params.csv')

	# We don't need the rogue peak info in the collection of results
	if 'extrapeak a' in patternresults_df.columns:
		patternresults_df = patternresults_df.drop('extrapeak a', 1)

	if 'extrapeak Vol' in patternresults_df.columns:
		patternresults_df = patternresults_df.drop('extrapeak Vol', 1)

	# Determine NRIXS experiment designation
	if np.all(patternresults_df['File name']=='FeNiSi_DAC13_P4a_13BMC_014'):
		patternresults_df['NRIXS exp'] = ['FeNiSi_DAC13_P4']
	elif np.all(patternresults_df['File name']=='FeNiSi_DAC14_P4b_028'):
		patternresults_df['NRIXS exp'] = ['FeNiSi_Ambient']
	else:
		patternresults_df['NRIXS exp']=patternresults_df['File name'].str.replace(r'[abd]*_[0-9]+', '')

	# Gather the uncertainties of parameters into the pattern results df
	# and collect results by phase
	if not paramresults_df[paramresults_df['Parameters']=='bcc Fe a'].empty:
		std_a = paramresults_df[paramresults_df['Parameters']=='bcc Fe a']['Std'].iloc[0]
		patternresults_df['bcc Fe da'] = std_a
		bccFe_df = bccFe_df.append(patternresults_df)
	
	if not paramresults_df[paramresults_df['Parameters']=='bcc FeNi a'].empty:
		std_a = paramresults_df[paramresults_df['Parameters']=='bcc FeNi a']['Std'].iloc[0]
		patternresults_df['bcc FeNi da'] = std_a
		bccFeNi_df = bccFeNi_df.append(patternresults_df)

	if not paramresults_df[paramresults_df['Parameters']=='bcc FeNiSi a'].empty:
		std_a = paramresults_df[paramresults_df['Parameters']=='bcc FeNiSi a']['Std'].iloc[0]
		patternresults_df['bcc FeNiSi da'] = std_a
		bccFeNiSi_df = bccFeNiSi_df.append(patternresults_df)
	
	# if not paramresults_df[paramresults_df['Parameters']=='hcp Fe a'].empty:
	# 	std_a = paramresults_df[paramresults_df['Parameters']=='hcp Fe a']['Std'].iloc[0]
	# 	patternresults_df['hcp Fe da'] = std_a
	# 	std_c = paramresults_df[paramresults_df['Parameters']=='hcp Fe c']['Std'].iloc[0]
	# 	patternresults_df['hcp Fe dc'] = std_c

	if not paramresults_df[paramresults_df['Parameters']=='hcp FeNi a'].empty:
		std_a = paramresults_df[paramresults_df['Parameters']=='hcp FeNi a']['Std'].iloc[0]
		patternresults_df['hcp FeNi da'] = std_a
		std_c = paramresults_df[paramresults_df['Parameters']=='hcp FeNi c']['Std'].iloc[0]
		patternresults_df['hcp FeNi dc'] = std_c
		hcpFeNi_df = hcpFeNi_df.append(patternresults_df)

	if not paramresults_df[paramresults_df['Parameters']=='hcp FeNiSi a'].empty:
		std_a = paramresults_df[paramresults_df['Parameters']=='hcp FeNiSi a']['Std'].iloc[0]
		patternresults_df['hcp FeNiSi da'] = std_a
		std_c = paramresults_df[paramresults_df['Parameters']=='hcp FeNiSi c']['Std'].iloc[0]
		patternresults_df['hcp FeNiSi dc'] = std_c
		hcpFeNiSi_df = hcpFeNiSi_df.append(patternresults_df)

	allresults_df = allresults_df.append(patternresults_df)

# Rename columns to ease later calculations and fix index so it's not all zeros
bccFe_df = bccFe_df.rename(columns={'bcc Fe a': 'a', 'bcc Fe da': 'da',
									'bcc Fe Vol': 'V'}).reset_index(drop=True)
bccFeNi_df = bccFeNi_df.rename(columns={'bcc FeNi a': 'a', 'bcc FeNi da': 'da',
									'bcc FeNi Vol': 'V'}).reset_index(drop=True)
bccFeNiSi_df = bccFeNiSi_df.rename(columns={'bcc FeNiSi a': 'a', 'bcc FeNiSi da': 'da',
									'bcc FeNiSi Vol': 'V'}).reset_index(drop=True)
hcpFeNi_df = hcpFeNi_df.rename(columns={'hcp FeNi a': 'a', 'hcp FeNi da': 'da',
	'hcp FeNi c':'c', 'hcp FeNi dc': 'dc',
	'hcp FeNi Vol':'V', 'hcp FeNi c/a':'c/a'}).reset_index(drop=True)
hcpFeNiSi_df = hcpFeNiSi_df.rename(columns={'hcp FeNiSi a': 'a', 'hcp FeNiSi da': 'da',
	'hcp FeNiSi c':'c', 'hcp FeNiSi dc': 'dc',
	'hcp FeNiSi Vol':'V', 'hcp FeNiSi c/a':'c/a'}).reset_index(drop=True)
allresults_df = allresults_df.reset_index(drop=True)

# Calculate P using equations of state in FeAlloyEOS.csv
bccFe_df['P']     = VinetEOS(bccFe_df['V'],    *getEOSparams(EOS_df, 'bcc Fe'))
bccFeNi_df['P']   = VinetEOS(bccFeNi_df['V'],  *getEOSparams(EOS_df, 'bcc FeNi'))
bccFeNiSi_df['P'] = VinetEOS(bccFeNiSi_df['V'],*getEOSparams(EOS_df, 'bcc FeNiSi'))
hcpFeNi_df['P']   = VinetEOS(hcpFeNi_df['V'],  *getEOSparams(EOS_df, 'hcp FeNi'))
hcpFeNiSi_df['P'] = VinetEOS(hcpFeNiSi_df['V'],*getEOSparams(EOS_df, 'hcp FeNiSi'))

# Save results
allresults_df.to_csv('Results/XRD_results_all.csv',index=False)
bccFe_df.to_csv('Results/XRD_results_all_bccFe.csv',index=False)
bccFeNi_df.to_csv('Results/XRD_results_all_bccFeNi.csv',index=False)
bccFeNiSi_df.to_csv('Results/XRD_results_all_bccFeNiSi.csv',index=False)
hcpFeNi_df.to_csv('Results/XRD_results_all_hcpFeNi.csv',index=False)
hcpFeNiSi_df.to_csv('Results/XRD_results_all_hcpFeNiSi.csv',index=False)

# Reduce to NRIXS experiment space (average lat params)
# Also calculate range of V and P
bccFe_reduced_df = bccFe_df[['NRIXS exp']].copy().drop_duplicates().reset_index(drop=True)
bccFeNi_reduced_df = bccFeNi_df[['NRIXS exp']].copy().drop_duplicates().reset_index(drop=True)
bccFeNiSi_reduced_df = bccFeNiSi_df[['NRIXS exp']].copy().drop_duplicates().reset_index(drop=True)
hcpFeNi_reduced_df = hcpFeNi_df[['NRIXS exp']].copy().drop_duplicates().reset_index(drop=True)
hcpFeNiSi_reduced_df = hcpFeNiSi_df[['NRIXS exp']].copy().drop_duplicates().reset_index(drop=True)

def reduce_latpar(new_df, old_df):
	# Initialize variables
	a = np.zeros(len(new_df))
	da = np.zeros(len(new_df))
	if 'c' in old_df.columns:
		c = np.zeros(len(new_df))
		dc = np.zeros(len(new_df))
	Vrange = np.zeros(len(new_df))
	Prange = np.zeros(len(new_df))
	for i, NRIXS_exp in enumerate(new_df['NRIXS exp']):
		# average lattice parameter a and propagate uncertainties to da
		mult_a = old_df[old_df['NRIXS exp']==NRIXS_exp]['a']
		mult_da = old_df[old_df['NRIXS exp']==NRIXS_exp]['da']
		a[i] = np.mean(mult_a)
		da[i] = np.linalg.norm(mult_da)/len(mult_da)
		if 'c' in old_df.columns:
				# average lattice parameter c and propagate uncertainties to dc
			mult_c = old_df[old_df['NRIXS exp']==NRIXS_exp]['c']
			mult_dc = old_df[old_df['NRIXS exp']==NRIXS_exp]['dc']
			c[i] = np.mean(mult_c)
			dc[i] = np.linalg.norm(mult_dc)/len(mult_da)
		mult_V = old_df[old_df['NRIXS exp']==NRIXS_exp]['V']
		mult_P = old_df[old_df['NRIXS exp']==NRIXS_exp]['P']
		Vrange[i] = max(mult_V) - min(mult_V)
		Prange[i] = max(mult_P) - min(mult_P)
	new_df['a'] = a
	new_df['da'] = da
	if 'c' in old_df.columns:
		new_df['c'] = c
		new_df['dc'] = dc
	new_df['V range'] = Vrange
	new_df['P range'] = Prange
	return new_df

bccFe_reduced_df     = reduce_latpar(bccFe_reduced_df, bccFe_df)
bccFeNi_reduced_df   = reduce_latpar(bccFeNi_reduced_df, bccFeNi_df)
bccFeNiSi_reduced_df = reduce_latpar(bccFeNiSi_reduced_df, bccFeNiSi_df)
hcpFeNi_reduced_df   = reduce_latpar(hcpFeNi_reduced_df, hcpFeNi_df)
hcpFeNiSi_reduced_df = reduce_latpar(hcpFeNiSi_reduced_df, hcpFeNiSi_df)

# The XRD fitting program underestimates lattice parameter uncertainty
# See lab book pg. 25-27 (1/16/18) for a detailed explanation for why these values were chosen
Fe_scaleparam = 100
bccFeNi_scaleparam = 40
hcpFeNi_scaleparam = 8
bccFeNiSi_scaleparam = 20
hcpFeNiSi_scaleparam = 4
 
bccFe_reduced_df['da']     = Fe_scaleparam * bccFe_reduced_df['da']    
bccFeNi_reduced_df['da']   = bccFeNi_scaleparam * bccFeNi_reduced_df['da']  
bccFeNiSi_reduced_df['da'] = bccFeNiSi_scaleparam * bccFeNiSi_reduced_df['da']
hcpFeNi_reduced_df['da']   = hcpFeNi_scaleparam * hcpFeNi_reduced_df['da']  
hcpFeNiSi_reduced_df['da'] = hcpFeNiSi_scaleparam * hcpFeNiSi_reduced_df['da']
hcpFeNi_reduced_df['dc']   = hcpFeNi_scaleparam * hcpFeNi_reduced_df['dc']  
hcpFeNiSi_reduced_df['dc'] = hcpFeNiSi_scaleparam * hcpFeNiSi_reduced_df['dc']

# Calculate V
bccFe_reduced_df['V'] = calc_V_bcc(bccFe_reduced_df['a'])
bccFeNi_reduced_df['V'] = calc_V_bcc(bccFeNi_reduced_df['a'])
bccFeNiSi_reduced_df['V'] = calc_V_bcc(bccFeNiSi_reduced_df['a'])
hcpFeNi_reduced_df['V'] = calc_V_hcp(hcpFeNi_reduced_df['a'],hcpFeNi_reduced_df['c'])
hcpFeNiSi_reduced_df['V'] = calc_V_hcp(hcpFeNiSi_reduced_df['a'],hcpFeNiSi_reduced_df['c'])

# Calculate V uncertainty
bccFe_reduced_df['dV'] = calc_dV_bcc(bccFe_reduced_df['a'],bccFe_reduced_df['da'])
bccFeNi_reduced_df['dV'] = calc_dV_bcc(bccFeNi_reduced_df['a'],bccFeNi_reduced_df['da'])
bccFeNiSi_reduced_df['dV'] = calc_dV_bcc(bccFeNiSi_reduced_df['a'],bccFeNiSi_reduced_df['da'])
hcpFeNi_reduced_df['dV'] = calc_dV_hcp(hcpFeNi_reduced_df['a'],hcpFeNi_reduced_df['c'],
			hcpFeNi_reduced_df['da'],hcpFeNi_reduced_df['dc'])
hcpFeNiSi_reduced_df['dV'] = calc_dV_hcp(hcpFeNiSi_reduced_df['a'],hcpFeNiSi_reduced_df['c'],
			hcpFeNiSi_reduced_df['da'],hcpFeNiSi_reduced_df['dc'])

# Calculate density
Fe_M = EOS_df[EOS_df['Phase']=='bcc Fe']['Molecular mass'].iloc[0]
FeNi_M = EOS_df[EOS_df['Phase']=='bcc FeNi']['Molecular mass'].iloc[0]
FeNiSi_M = EOS_df[EOS_df['Phase']=='bcc FeNiSi']['Molecular mass'].iloc[0]

bccFe_reduced_df['rho'],    bccFe_reduced_df['drho']     = calc_rho(bccFe_reduced_df['V'],
															bccFe_reduced_df['dV'],Fe_M)
bccFeNi_reduced_df['rho'],  bccFeNi_reduced_df['drho']   = calc_rho(bccFeNi_reduced_df['V'],
															bccFeNi_reduced_df['dV'],FeNi_M)
bccFeNiSi_reduced_df['rho'],bccFeNiSi_reduced_df['drho'] = calc_rho(bccFeNiSi_reduced_df['V'],
															bccFeNiSi_reduced_df['dV'],FeNiSi_M)
hcpFeNi_reduced_df['rho'],  hcpFeNi_reduced_df['drho']   = calc_rho(hcpFeNi_reduced_df['V'],
															hcpFeNi_reduced_df['dV'],FeNi_M)
hcpFeNiSi_reduced_df['rho'],hcpFeNiSi_reduced_df['drho'] = calc_rho(hcpFeNiSi_reduced_df['V'],
															hcpFeNiSi_reduced_df['dV'],FeNiSi_M)


# Calculate P. Calculating P uncertainty requires EOS param
# correlations.
bccFe_reduced_df['P']     = VinetEOS(bccFe_reduced_df['V'],    *getEOSparams(EOS_df, 'bcc Fe'))
bccFeNi_reduced_df['P']   = VinetEOS(bccFeNi_reduced_df['V'],  *getEOSparams(EOS_df, 'bcc FeNi'))
bccFeNiSi_reduced_df['P'] = VinetEOS(bccFeNiSi_reduced_df['V'],*getEOSparams(EOS_df, 'bcc FeNiSi'))
hcpFeNi_reduced_df['P']   = VinetEOS(hcpFeNi_reduced_df['V'],  *getEOSparams(EOS_df, 'hcp FeNi'))
hcpFeNiSi_reduced_df['P'] = VinetEOS(hcpFeNiSi_reduced_df['V'],*getEOSparams(EOS_df, 'hcp FeNiSi'))

bccFe_reduced_df['dP']       = [calc_dP_VinetEOS(V,dV,EOS_df, 'bcc Fe') for V, dV 
							in zip(bccFe_reduced_df['V'].values,bccFe_reduced_df['dV'].values)]
bccFeNi_reduced_df['dP']     = [calc_dP_VinetEOS(V,dV,EOS_df, 'bcc FeNi') for V, dV 
							in zip(bccFeNi_reduced_df['V'].values,bccFeNi_reduced_df['dV'].values)]
bccFeNiSi_reduced_df['dP']     = [calc_dP_VinetEOS(V,dV,EOS_df, 'bcc FeNiSi') for V, dV 
							in zip(bccFeNiSi_reduced_df['V'].values,bccFeNiSi_reduced_df['dV'].values)]
hcpFeNi_reduced_df['dP']     = [calc_dP_VinetEOS(V,dV,EOS_df, 'hcp FeNi') for V, dV 
							in zip(hcpFeNi_reduced_df['V'].values,hcpFeNi_reduced_df['dV'].values)]
hcpFeNiSi_reduced_df['dP']     = [calc_dP_VinetEOS(V,dV,EOS_df, 'hcp FeNiSi') for V, dV 
							in zip(hcpFeNiSi_reduced_df['V'].values,hcpFeNiSi_reduced_df['dV'].values)]

# bccFe_reduced_df = bccFe_reduced_df.drop(['V range','P range'],1)
# bccFeNi_reduced_df = bccFeNi_reduced_df.drop(['V range','P range'],1)
# hcpFeNi_reduced_df = hcpFeNi_reduced_df.drop(['V range','P range'],1)
# bccFeNiSi_reduced_df = bccFeNiSi_reduced_df.drop(['V range','P range'],1)
# hcpFeNiSi_reduced_df = hcpFeNiSi_reduced_df.drop(['V range','P range'],1)

# We know the pressure for ambient measurements is actually 0 GPa, so let's fix that
bccFe_reduced_df.loc[bccFe_reduced_df['NRIXS exp']=='Fe_Ambient', 'P'] = 0
bccFe_reduced_df.loc[bccFe_reduced_df['NRIXS exp']=='Fe_Ambient', 'dP'] = 0
bccFeNi_reduced_df.loc[bccFeNi_reduced_df['NRIXS exp']=='FeNi_Ambient', 'P'] = 0
bccFeNi_reduced_df.loc[bccFeNi_reduced_df['NRIXS exp']=='FeNi_Ambient', 'dP'] = 0
bccFeNiSi_reduced_df.loc[bccFeNiSi_reduced_df['NRIXS exp']=='FeNiSi_Ambient', 'P'] = 0
bccFeNiSi_reduced_df.loc[bccFeNiSi_reduced_df['NRIXS exp']=='FeNiSi_Ambient', 'dP'] = 0

# # Import bulk modulus values (manually calculated with MINUTI)
# K_T_df = pd.read_csv('../015_BulkModulusDet/BulkMod.csv')
# bccFe_reduced_df['K_T'] = np.array(K_T_df[K_T_df['Phase']=='bcc Fe']['K_T'])
# bccFeNi_reduced_df['K_T'] = np.array(K_T_df[K_T_df['Phase']=='bcc FeNi']['K_T'])
# bccFeNiSi_reduced_df['K_T'] = np.array(K_T_df[K_T_df['Phase']=='bcc FeNiSi']['K_T'])
# hcpFeNi_reduced_df['K_T'] = np.array(K_T_df[K_T_df['Phase']=='hcp FeNi']['K_T'])
# hcpFeNiSi_reduced_df['K_T'] = np.array(K_T_df[K_T_df['Phase']=='hcp FeNiSi']['K_T'])
# # print(np.array(K_T_df[K_T_df['Phase']=='bcc FeNi']['K_T']))

# Round results so we're not saving a ridiculous number of sig figs
latdec = 5
Vdec = 3
Pdec = 2
# Kdec = 0
rhodec = 3
bccFe_reduced_df = bccFe_reduced_df.round({'a': latdec, 'da': latdec,'V': Vdec, 'dV': Vdec,
												'P': Pdec, 'dP': Pdec, #'K_T':Kdec,
												'V range': Vdec, 'P range': Pdec})
# For some reason, 'rho' and 'drho' aren't rounding with the above method...
bccFe_reduced_df[['rho','drho']] = np.round(bccFe_reduced_df[['rho','drho']].astype(float),decimals=rhodec)
bccFeNi_reduced_df = bccFeNi_reduced_df.round({'a': latdec, 'da': latdec,'V': Vdec, 'dV': Vdec,
												'P': Pdec, 'dP': Pdec, #'K_T':Kdec,
												'V range': Vdec, 'P range': Pdec})
bccFeNi_reduced_df[['rho','drho']] = np.round(bccFeNi_reduced_df[['rho','drho']].astype(float),decimals=rhodec)
bccFeNiSi_reduced_df = bccFeNiSi_reduced_df.round({'a': latdec, 'da': latdec,'V': Vdec, 'dV': Vdec,
												'P': Pdec, 'dP': Pdec, #'K_T':Kdec,
												'V range': Vdec, 'P range': Pdec})
bccFeNiSi_reduced_df[['rho','drho']] = np.round(bccFeNiSi_reduced_df[['rho','drho']].astype(float),decimals=rhodec)
hcpFeNi_reduced_df = hcpFeNi_reduced_df.round({'a': latdec, 'da': latdec,'c': latdec, 'dc': latdec,
												'V': Vdec, 'dV': Vdec,'P': Pdec, 'dP': Pdec, #'K_T':Kdec, 
												'V range': Vdec, 'P range': Pdec})
hcpFeNi_reduced_df[['rho','drho']] = np.round(hcpFeNi_reduced_df[['rho','drho']].astype(float),decimals=rhodec)
hcpFeNiSi_reduced_df = hcpFeNiSi_reduced_df.round({'a': latdec, 'da': latdec,'c': latdec, 'dc': latdec,
												'V': Vdec, 'dV': Vdec,'P': Pdec, 'dP': Pdec, #'K_T':Kdec, 
												'V range': Vdec, 'P range': Pdec})
hcpFeNiSi_reduced_df[['rho','drho']] = np.round(hcpFeNiSi_reduced_df[['rho','drho']].astype(float),decimals=rhodec)

print(bccFe_reduced_df)
print(bccFeNi_reduced_df)
print(bccFeNiSi_reduced_df)
print(hcpFeNi_reduced_df)
print(hcpFeNiSi_reduced_df)

# Save results
bccFe_reduced_df.to_csv('Results/XRD_results_bccFe.csv',index=False)
bccFeNi_reduced_df.to_csv('Results/XRD_results_bccFeNi.csv',index=False)
bccFeNiSi_reduced_df.to_csv('Results/XRD_results_bccFeNiSi.csv',index=False)
hcpFeNi_reduced_df.to_csv('Results/XRD_results_hcpFeNi.csv',index=False)
hcpFeNiSi_reduced_df.to_csv('Results/XRD_results_hcpFeNiSi.csv',index=False)