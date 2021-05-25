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

pdf_filename = 'pdf_results_Constrained_Power.csv'


# Functions
###########

# Gaussian function with an amplitude of 1 (not normalized)
def Gaussian(x,amplitude,mean,sigma):
	return amplitude*np.exp(-(x-mean)**2/(2*sigma**2))

# Lorentzian function with an amplitude of 1 (not normalized)
def Lorentzian(x,amplitude,mean,sigma):
    return amplitude*(sigma/2)**2/((x-mean)**2+(sigma/2)**2)

# Gaussian function with an amplitude of 1 (not normalized)
def Gaussian1(x,mean,sigma):
	return np.exp(-(x-mean)**2/(2*sigma**2))

# Lorentzian function with an amplitude of 1 (not normalized)
def Lorentzian1(x,mean,sigma):
    return (sigma/2)**2/((x-mean)**2+(sigma/2)**2)

# Pseudo-Voigt function
def pseudoVoigt(x,amplitude,mean,sigma,shape):
    return amplitude*( shape*Lorentzian1(x,mean,sigma) +
         (1-shape)*Gaussian1(x,mean,sigma) )

# Function for an asymmetric Gaussian peak
def asymGaussian(x,amplitude,mean,sigmaLeft,sigmaRight):
	# Piecewise function
	y = np.zeros(len(x)) # Framework for adding pieces of function
	y += Gaussian(x,amplitude,mean,sigmaLeft) * (x <= mean) # Left side
	y += Gaussian(x,amplitude,mean,sigmaRight) * (x > mean) # Right side
	return y

# Function for an asymmetric Lorentzian peak
def asymLorentzian(x,amplitude,mean,sigmaLeft,sigmaRight):
	# Piecewise function
	y = np.zeros(len(x)) # Framework for adding pieces of function
	y += Lorentzian(x,amplitude,mean,sigmaLeft) * (x <= mean) # Left side
	y += Lorentzian(x,amplitude,mean,sigmaRight) * (x > mean) # Right side
	return y

# Function for an asymmetric pseudo-Voigt peak
def asymPseudoVoigt(x,amplitude,mean,sigmaLeft,shapeLeft,sigmaRight,shapeRight):
	# Piecewise function
	y = np.zeros(len(x)) # Framework for adding pieces of function
	y += pseudoVoigt(x,amplitude,mean,sigmaLeft,shapeLeft) * (x <= mean) # Left side
	y += pseudoVoigt(x,amplitude,mean,sigmaRight,shapeRight) * (x > mean) # Right side
	return y

# Wolfgang's function
def WSfunction(x,amplitude,mean,sigmaLeft,shapeLeft,sigmaRight,shapeRight):
	# Piecewise function
	y = np.zeros(len(x)) # Framework for adding pieces of function
	y += amplitude*2**(-np.abs((x-mean)/sigmaLeft)**shapeLeft) * (x <= mean) # Left side
	y += amplitude*2**(-np.abs((x-mean)/sigmaRight)**shapeRight) * (x > mean) # Right side
	return y


# Plot Debye pdf
def plot_vD_pdf(binned_df,Loren_opt,aLoren_opt,folder):

	fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(16,4))

	# Plot vD pdf data
	ax0.plot(binned_df['Bin Midpoint'],binned_df['Prob'],'.')
	# ax0.plot(binned_df['Bin Midpoint'],Gaussian(binned_df['Bin Midpoint'],*Gauss_opt),
	# 	color='red')
	ax0.plot(binned_df['Bin Midpoint'],Lorentzian(binned_df['Bin Midpoint'],*Loren_opt),
		color='orange')
 	# ax0.plot(binned_df['Bin Midpoint'],pseudoVoigt(binned_df['Bin Midpoint'],*pVoigt_opt),
	# 	color='green')
	# ax0.plot(binned_df['Bin Midpoint'],asymGaussian(binned_df['Bin Midpoint'],*aGauss_opt),
	# 	color='blue')
	ax0.plot(binned_df['Bin Midpoint'],asymLorentzian(binned_df['Bin Midpoint'],*aLoren_opt),
		color='purple')
	# ax0.plot(binned_df['Bin Midpoint'],asymPseudoVoigt(binned_df['Bin Midpoint'],*apVoigt_opt),
	# 	color='gray')
	# ax0.plot(binned_df['Bin Midpoint'],WSfunction(binned_df['Bin Midpoint'],*WS_opt),
	# 	color='red')
	
	plt.tight_layout()
	fig.savefig('test_plots/vD_pdf_'+folder+'.pdf', format='pdf')
	plt.close()



# Fit asymmetric peak to pdf data
#################################

# Find the filepath of all .res NRIXS files in phox directories
pdf_path_list = [filepath for filepath in glob.glob('*/'+pdf_filename)]

# Prep lists, dictionaries
folder_list = []
vD_pdf_dict = dict()

# Collect folders, indices, paths, and input values
for pdf_path in pdf_path_list:
	folder = re.findall('([A-Za-z0-9_]+)/'+pdf_filename,pdf_path)[0]

	binned_df = pd.read_csv(pdf_path)

	# Parameter bounds
	amp_min = 0.1
	amp_max = 30
	mean_min = 100
	mean_max = 20000
	sigma_min = 1
	sigma_max = 1000
	shape_min = 0
	shape_max = 1

	# Initial parameter guess
	amp_guess = 1
	mean_guess = np.average(binned_df['Bin Midpoint'],weights = binned_df['Prob'])
	sigma_guess = 50
	shape_guess = 0.5

	print('Fitting peak to '+folder)

	# Fit Gaussian peak to pdf
	Gauss_opt,cov = curve_fit(Gaussian,binned_df['Bin Midpoint'],binned_df['Prob'],
		p0=[amp_guess,mean_guess,sigma_guess],
		bounds = ([amp_min,mean_min,sigma_min],
			[amp_max,mean_max,sigma_max]))
	print(Gauss_opt)

	# Fit Lorentzian peak to pdf
	Loren_opt,cov = curve_fit(Lorentzian,binned_df['Bin Midpoint'],binned_df['Prob'],
		p0=[amp_guess,mean_guess,sigma_guess],
		bounds = ([amp_min,mean_min,sigma_min],
			[amp_max,mean_max,sigma_max]))
	print(Loren_opt)

	# # Fit pseudo-Voigt peak to pdf
	# pVoigt_opt,cov = curve_fit(pseudoVoigt,binned_df['Bin Midpoint'],binned_df['Prob'],
	# 	p0=[amp_guess,mean_guess,sigma_guess,shape_guess],
	# 	bounds = ([amp_min,mean_min,sigma_min,shape_min],
	# 		[amp_max,mean_max,sigma_max,shape_max]))

	# # Fit asymmetric Gaussian peak to pdf
	# aGauss_opt,cov = curve_fit(asymGaussian,binned_df['Bin Midpoint'],binned_df['Prob'],
	# 	p0=[amp_guess,mean_guess,sigma_guess,sigma_guess],
	# 	bounds = ([amp_min,mean_min,sigma_min,sigma_min],
	# 		[amp_max,mean_max,sigma_max,sigma_max]))

	# Fit asymmetric Lorentzian peak to pdf
	aLoren_opt,cov = curve_fit(asymLorentzian,binned_df['Bin Midpoint'],binned_df['Prob'],
		p0=[amp_guess,mean_guess,sigma_guess,sigma_guess],
		bounds = ([amp_min,mean_min,sigma_min,sigma_min],
			[amp_max,mean_max,sigma_max,sigma_max]))
	print(aLoren_opt)

	# # Fit asymmetric pseudo-Voigt peak to pdf
	# apVoigt_opt,cov = curve_fit(asymPseudoVoigt,binned_df['Bin Midpoint'],binned_df['Prob'],
	# 	p0=[amp_guess,mean_guess,sigma_guess,shape_guess,sigma_guess,shape_guess],
	# 	bounds = ([amp_min,mean_min,sigma_min,shape_min,sigma_min,shape_min],
	# 		[amp_max,mean_max,sigma_max,shape_max,sigma_max,shape_max]))

	# # Fit Wolfgang's function to pdf
	# WS_opt,cov = curve_fit(WSfunction,binned_df['Bin Midpoint'],binned_df['Prob'],
	# 	p0=[amp_guess,mean_guess,sigma_guess,shape_guess,sigma_guess,shape_guess],
	# 	bounds = ([amp_min,mean_min,sigma_min,shape_min,sigma_min,shape_min],
	# 		[amp_max,mean_max,sigma_max,shape_max,sigma_max,shape_max]))

	# Create plot of pdf
	plot_vD_pdf(binned_df,Loren_opt,aLoren_opt,folder)


