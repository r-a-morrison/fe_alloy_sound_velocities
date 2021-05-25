# Front matter
import os
import glob
import re
import pandas as pd
import numpy as np
import scipy.constants as constants

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec

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


# Functions
###########

# Calculate Ks
def calc_Ks(rho,vphi):
    return vphi**2*rho

# Calculate dKs
def calc_dKs(rho,drho,vphi,dvphi):
    return np.sqrt( (2*vphi*rho*dvphi)**2 + (vphi**2*drho)**2 )

# Given 3 curves each with uncertainty (where which curve lies in the middle is known),
# determine the overall range spanned by these uncertainties. Return the overall
# uncertainty.
def calc_uncertainty(result_df,a_df,b_df,quantity):
	result_df = result_df[['P',quantity]]
	a_df = a_df[['P',quantity,'d'+quantity]]
	a_df = a_df.rename(columns={quantity: quantity+'_a','d'+quantity: 'd'+quantity+'_a'})
	b_df = b_df[['P',quantity,'d'+quantity]]
	b_df = b_df.rename(columns={quantity: quantity+'_b','d'+quantity: 'd'+quantity+'_b'})
	comb_df = pd.merge(a_df,b_df,on='P')
	comb_df = pd.merge(result_df,comb_df,on='P')

	comb_df[quantity+'_a_upper'] = comb_df[quantity+'_a'] + comb_df['d'+quantity+'_a']
	comb_df[quantity+'_a_lower'] = comb_df[quantity+'_a'] - comb_df['d'+quantity+'_a']
	comb_df[quantity+'_b_upper'] = comb_df[quantity+'_b'] + comb_df['d'+quantity+'_b']
	comb_df[quantity+'_b_lower'] = comb_df[quantity+'_b'] - comb_df['d'+quantity+'_b']
	comb_df[quantity+'_upper'] = comb_df[[quantity+'_a_upper',quantity+'_b_upper']].max(axis=1)
	comb_df[quantity+'_lower'] = comb_df[[quantity+'_a_lower',quantity+'_b_lower']].min(axis=1)
	comb_df['d'+quantity] = ((comb_df[quantity+'_upper'] - comb_df[quantity]) + 
		(comb_df[quantity] - comb_df[quantity+'_lower']))/2
	return comb_df['d'+quantity]


# Define data sets and input parameters
#######################################
studylist = []
phase = dict()
MINUTIpath = dict()
XRDpath = dict()

# bcc Fe
study = 'bccFe' # Our NRIXS data
studylist.append(study)
phase[study] = 'bcc'
MINUTIpath[study] = 'bccFe_Dewaele/Output/bccFe_Dewaele' # Dewaele's EOS data
XRDpath[study] = '../010_XRDAnalysis/Results/XRD_results_bccFe.csv'

# bcc FeNi
study = 'bccFeNi'
studylist.append(study)
phase[study] = 'bcc'
MINUTIpath[study] = 'bccFeNi/Output/bccFeNi'
XRDpath[study] = '../010_XRDAnalysis/Results/XRD_results_bccFeNi.csv'


# bcc FeNiSi
study = 'bccFeNiSi'
studylist.append(study)
phase[study] = 'bcc'
MINUTIpath[study] = 'bccFeNiSi/Output/bccFeNiSi'
XRDpath[study] = '../010_XRDAnalysis/Results/XRD_results_bccFeNiSi.csv'


# hcp Fe
study = 'hcpFe_Murphy' # Caitlin's NRIXS data
studylist.append(study)
phase[study] = 'hcp'
MINUTIpath[study] = 'hcpFe_Dewaele/fromData/Output/hcpFe_Dewaele' # Dewaele's EOS data
XRDpath[study] = '../005_PubNRIXSVals/hcpFe_Murphy.csv'

# hcp Fe: gammahigh
study = 'hcpFe_Murphy_gammahigh'
studylist.append(study)
phase[study] = 'hcp'
MINUTIpath[study] = 'hcpFe_Dewaele/fromData_gammahigh/Output/hcpFe_Dewaele'
XRDpath[study] = '../005_PubNRIXSVals/hcpFe_Murphy.csv'

# hcp Fe: gammalow
study = 'hcpFe_Murphy_gammalow'
studylist.append(study)
phase[study] = 'hcp'
MINUTIpath[study] = 'hcpFe_Dewaele/fromData_gammalow/Output/hcpFe_Dewaele'
XRDpath[study] = '../005_PubNRIXSVals/hcpFe_Murphy.csv'

# hcp FeNi
study = 'hcpFeNi'
studylist.append(study)
phase[study] = 'hcp'
MINUTIpath[study] = 'hcpFeNi/fromData/Output/hcpFeNi'
XRDpath[study] = '../010_XRDAnalysis/Results/XRD_results_hcpFeNi.csv'

# hcp FeNi: gammahigh
study = 'hcpFeNi_gammahigh'
studylist.append(study)
phase[study] = 'hcp'
MINUTIpath[study] = 'hcpFeNi/fromData_gammahigh/Output/hcpFeNi'
XRDpath[study] = '../010_XRDAnalysis/Results/XRD_results_hcpFeNi.csv'

# hcp FeNi: gammalow
study = 'hcpFeNi_gammalow'
studylist.append(study)
phase[study] = 'hcp'
MINUTIpath[study] = 'hcpFeNi/fromData_gammalow/Output/hcpFeNi'
XRDpath[study] = '../010_XRDAnalysis/Results/XRD_results_hcpFeNi.csv'

# hcp FeNiSi
study = 'hcpFeNiSi'
studylist.append(study)
phase[study] = 'hcp'
MINUTIpath[study] = 'hcpFeNiSi/fromData/Output/hcpFeNiSi'
XRDpath[study] = '../010_XRDAnalysis/Results/XRD_results_hcpFeNiSi.csv'

# hcp FeNiSi: gammahigh
study = 'hcpFeNiSi_gammahigh'
studylist.append(study)
phase[study] = 'hcp'
MINUTIpath[study] = 'hcpFeNiSi/fromData_gammahigh/Output/hcpFeNiSi'
XRDpath[study] = '../010_XRDAnalysis/Results/XRD_results_hcpFeNiSi.csv'

# hcp FeNiSi: gammalow
study = 'hcpFeNiSi_gammalow'
studylist.append(study)
phase[study] = 'hcp'
MINUTIpath[study] = 'hcpFeNiSi/fromData_gammalow/Output/hcpFeNiSi'
XRDpath[study] = '../010_XRDAnalysis/Results/XRD_results_hcpFeNiSi.csv'


# Load data and calculate Ks
############################
XRD_dfdict = dict()
MINUTI_dfdict = dict()

for study in studylist:
	# Import XRD data
	init_XRD_df = pd.read_csv(XRDpath[study],engine='python')
	# Change index name to match other data sets
	if 'NRIXS exp' in init_XRD_df:
	    init_XRD_df = init_XRD_df.rename(columns={'NRIXS exp': 'Index'})
	# We only need some of the info here
	if phase[study] == 'bcc':
	    XRD_df = init_XRD_df[['Index','a','da','V','dV','rho','drho','P','dP']]
	else:
	    XRD_df = init_XRD_df[['Index','a','da','c','dc','V','dV','rho','drho','P','dP']]
	# Be careful here: Density here is enriched. We want enriched density for sound 
	# velocity calcs, but we need natural density for K_S calculations.

	# Import MINUTI data
	KT_df      = pd.read_csv(MINUTIpath[study]+'_bms.dat',header=None,
				 names = ['P','KT','dKT'],sep='\s+',engine='python')
	KTprime_df = pd.read_csv(MINUTIpath[study]+'_bdp.dat',header=None,
				 names = ['P','KTprime','dKTprime'],sep='\s+',engine='python')
	# Density is natural abundance; need to get back to K_S
	rho_df = pd.read_csv(MINUTIpath[study]+'_dns.dat',header=None,
				 names = ['P','rho','drho'],sep='\s+',engine='python')
	vphi_df    = pd.read_csv(MINUTIpath[study]+'_smv.dat',header=None,
				 names = ['P','vphi','dvphi'],sep='\s+',engine='python')
	# Combine dfs into one df
	init_MINUTI_df = pd.merge(KT_df,KTprime_df,on='P')
	init_MINUTI_df = pd.merge(init_MINUTI_df,rho_df,on='P')
	init_MINUTI_df = pd.merge(init_MINUTI_df,vphi_df,on='P')

	# Select pressures in MINUTI_df that correspond with XRD pressures and give them an index
	# What pressures correspond to our measurements?
	P_df = XRD_df[['Index','P']]
	# Round P so we can look them up in MINUTI_df
	P_df = P_df.round({'P': 1})
	# Get rows in MINUTI_df that correspond to our measurements, which are listed in P_df
	# (We're not interested in any of the other pressures.)
	MINUTI_df = init_MINUTI_df[init_MINUTI_df['P'].isin(P_df['P'].values)]
	# Add the index designation to each line in MINUTI_df, so we can easily link it back to
	# our dataset
	MINUTI_df = pd.merge(MINUTI_df,P_df,on='P')

	# Calculate K_S
	MINUTI_df['KS'] = calc_Ks(MINUTI_df['rho'],MINUTI_df['vphi'])
	MINUTI_df['dKS'] = calc_dKs(MINUTI_df['rho'],MINUTI_df['drho'],
						MINUTI_df['vphi'],MINUTI_df['dvphi'])

	# Add results to dictionaries so we can easily access them again
	XRD_dfdict[study] = XRD_df
	MINUTI_dfdict[study] = MINUTI_df


# Calculate K_S uncertainty for hcp phases
##########################################

# We'll assume K_T = K_S for bcc phases. For hcp phases, we need to account for the effect
# of parameters gamma0, q, DebyeT on the uncertainty of K_S

results_dfdict = dict() # To save final results in before saving to csv

for study in ['hcpFe_Murphy','hcpFeNi','hcpFeNiSi']:
	XRD_df = XRD_dfdict[study]
	MINUTI_df = MINUTI_dfdict[study]
	gammahigh_df = MINUTI_dfdict[study+'_gammahigh']
	gammalow_df = MINUTI_dfdict[study+'_gammalow']

	# Calculate dvphi and dKS uncertainty accounting for gamma0, q, DebyeT uncertainties
	# Update dvphi and dKS in df
	MINUTI_df['dvphi'] = calc_uncertainty(MINUTI_df,gammahigh_df,gammalow_df,'vphi')
	MINUTI_df['dKS'] = calc_uncertainty(MINUTI_df,gammahigh_df,gammalow_df,'KS')

	# These are the new results we want to use
	# NOTE: vphi is NOT stored because it assumes natural enrichment, but we want
	# the enriched vphi
	results_df = MINUTI_df[['Index','KT','dKT','KTprime','dKTprime','KS','dKS']]
	# Combine with existing XRD data
	results_df = pd.merge(XRD_df,results_df,on='Index')

	# Note that we're using the density from XRD (enriched) in our final results. This is
	# what we want to use in our next step to calculate sound velocities.

	# Save results to dfdict so we can save results as csv
	results_dfdict[study] = results_df


# Add final results for bcc phases to results_df
################################################

# We'll assume K_S = K_T for bcc phases. The K_S calculated from MINUTI has the wrong
# gamma0, q, and DebyeT entered, so the K_S in MINUTI_df is wrong. That means vphi
# is wrong too.

for study in ['bccFe','bccFeNi','bccFeNiSi']:
	XRD_df = XRD_dfdict[study]
	MINUTI_df = MINUTI_dfdict[study]

	# These are the new results we want to use
	results_df = MINUTI_df[['Index','KT','dKT','KTprime','dKTprime']]
	# Combine with existing XRD data
	results_df = pd.merge(XRD_df,results_df,on='Index')

	# Save results to dfdict so we can save results as csv
	results_dfdict[study] = results_df


# Save results to csv
#####################

for study in ['bccFe','bccFeNi','bccFeNiSi','hcpFe_Murphy','hcpFeNi','hcpFeNiSi']:
	results_df = results_dfdict[study]

	# Round results so we're not saving a ridiculous number of sig figs
	latdec = 5
	Vdec = 3
	Pdec = 2
	rhodec = 3
	Kdec = 1
	Kprimedec = 2
	vdec = 3
	if phase[study]=='bcc':
		results_df = results_df.round({'a': latdec, 'da': latdec,
			'V': Vdec, 'dV': Vdec, 'rho': rhodec, 'drho': rhodec,
			'P': Pdec, 'dP': Pdec,'KT': Kdec, 'dKT': Kdec,'KTprime': Kprimedec,
			'dKTprime': Kprimedec})
		# For some obnoxious reason the above method doesn't work when data is saved to a csv.
		# The method below doesn't work either. Trailing zeros are lost.
		# results_df[['a','da']] = np.round(results_df[['a','da']].astype(float),decimals=latdec)
		# results_df[['V','dV']] = np.round(results_df[['V','dV']].astype(float),decimals=Vdec)
		# results_df[['P','dP']] = np.round(results_df[['P','dP']].astype(float),decimals=Pdec)
		# results_df[['rho','drho']] = np.round(results_df[['rho','drho']].astype(float),decimals=rhodec)
		# results_df[['KT','dKT']] = np.round(results_df[['KT','dKT']].astype(float),decimals=Kdec)
		# results_df[['KTprime','dKTprime']] = np.round(results_df[['KTprime','dKTprime']].astype(float),decimals=Kprimedec)
	if phase[study]=='hcp':
		results_df = results_df.round({'a': latdec, 'da': latdec, 'c': latdec, 'dc': latdec,
			'V': Vdec, 'dV': Vdec, 'rho': rhodec, 'drho': rhodec,
			'P': Pdec, 'dP': Pdec,'KT': Kdec, 'dKT': Kdec,'KTprime': Kprimedec,
			'dKTprime': Kprimedec,'KS': Kdec, 'dKS': Kdec})
		# results_df[['a','da','c','dc']] = np.round(results_df[['a','da','c','dc']].astype(float),decimals=latdec)
		# results_df[['V','dV']] = np.round(results_df[['V','dV']].astype(float),decimals=Vdec)
		# results_df[['P','dP']] = np.round(results_df[['P','dP']].astype(float),decimals=Pdec)
		# results_df[['rho','drho']] = np.round(results_df[['rho','drho']].astype(float),decimals=rhodec)
		# results_df[['KT','dKT','KS','dKS']] = np.round(results_df[['KT','dKT','KS','dKS']].astype(float),decimals=Kdec)
		# results_df[['KTprime','dKTprime']] = np.round(results_df[['KTprime','dKTprime']].astype(float),decimals=Kprimedec)
		# results_df[['vphi','dvphi']] = np.round(results_df[['vphi','dvphi']].astype(float),decimals=vdec)

	results_df.to_csv('Results/'+study+'.csv',index=False)
