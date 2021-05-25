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


# Input data
############

pdf_path = dict() # Key is tuple of composition and pressure
vD = dict()
dvD = dict()

# Fe
key = ('Fe',30)
pdf_path[key] = '../../../psvl_binwidth_test/vD_Fe_PowerLaw_bin10/2009Oct_30GPa/pdf_results_Constrained_Power.csv'
vD[key] = 4.421
dvD[key] = 0.063

key = ('Fe',53)
pdf_path[key] = '../../../psvl_binwidth_test/vD_Fe_PowerLaw_bin10/2009Oct_53GPa/pdf_results_Constrained_Power.csv'
vD[key] = 4.609
dvD[key] = 0.044

key = ('Fe',90)
pdf_path[key] = '../../../psvl_binwidth_test/vD_Fe_PowerLaw_bin10/2009Oct_90GPa/pdf_results_Constrained_Power.csv'
vD[key] = 5.132
dvD[key] = 0.029


# FeNi
key = ('FeNi',23)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNi_PowerLaw_bin10/2012Apr_DAC11_P4/pdf_results_Constrained_Power.csv'
vD[key] = 3.936
dvD[key] = 0.034

key = ('FeNi',41)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNi_PowerLaw_bin10/2012Apr_DAC11_P5/pdf_results_Constrained_Power.csv'
vD[key] = 4.243
dvD[key] = 0.013

key = ('FeNi',48)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNi_PowerLaw_bin10/2012Apr_DAC11_P6/pdf_results_Constrained_Power.csv'
vD[key] = 4.359
dvD[key] = 0.035

key = ('FeNi',63)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNi_PowerLaw_bin10/2012Apr_DAC11_P7/pdf_results_Constrained_Power.csv'
vD[key] = 4.519
dvD[key] = 0.027

key = ('FeNi',83)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNi_PowerLaw_bin10/2012Apr_DAC11_P9/pdf_results_Constrained_Power.csv'
vD[key] = 4.672
dvD[key] = 0.037

key = ('FeNi',104)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNi_PowerLaw_bin10/2012Oct_DAC11_P14/pdf_results_Constrained_Power.csv'
vD[key] = 4.972
dvD[key] = 0.036


# FeNiSi
key = ('FeNiSi',28)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNiSi_PowerLaw_bin10/2014Feb_DAC13_P1/pdf_results_Constrained_Power.csv'
vD[key] = 3.972
dvD[key] = 0.045

key = ('FeNiSi',37)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNiSi_PowerLaw_bin10/2014Feb_DAC13_P2/pdf_results_Constrained_Power.csv'
vD[key] = 4.158
dvD[key] = 0.035

key = ('FeNiSi',41)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNiSi_PowerLaw_bin10/2015Mar_DAC13_P3/pdf_results_Constrained_Power.csv'
vD[key] = 4.302
dvD[key] = 0.054

key = ('FeNiSi',55)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNiSi_PowerLaw_bin10/2015Mar_DAC13_P4/pdf_results_Constrained_Power.csv'
vD[key] = 4.373
dvD[key] = 0.041

key = ('FeNiSi',69)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNiSi_PowerLaw_bin10/2015Mar_DAC13_P5/pdf_results_Constrained_Power.csv'
vD[key] = 4.546
dvD[key] = 0.056

key = ('FeNiSi',86)
pdf_path[key] = '../../../psvl_binwidth_test/vD_FeNiSi_PowerLaw_bin10/2015Mar_DAC13_P6/pdf_results_Constrained_Power.csv'
vD[key] = 4.683
dvD[key] = 0.043


# Functions
###########

def make_subplot(study,P,ax,xmin,xmax,xaxislabelon,yaxislabelon):
	key = (study,P)
	pdf_df = pdf_dict[key]
	ymin = min(pdf_df['Prob']) - 0.05*max(pdf_df['Prob'])
	ymax = 1.05*max(pdf_df['Prob'])
	ax.plot(pdf_df['Bin Midpoint']/1000,pdf_df['Prob'],'.',markersize=7)
	ax.plot(vD[key]*np.ones(2), [ymin,ymax], color='red',
		linewidth=2)
	ax.fill_between([vD[key]-dvD[key],vD[key]+dvD[key]], ymin*np.ones(2), ymax*np.ones(2),
		facecolor='red', alpha=0.35)

	ax.set_xlim([xmin,xmax])
	ax.set_ylim([ymin,ymax])
	ax.tick_params(direction='in',left='off',top='on')
	if xaxislabelon:
		ax.set_xlabel(r'$v_D$ (km/s)',fontsize=18)
	else:
		ax.xaxis.set_ticklabels([])
	if yaxislabelon:
		ax.set_ylabel(r'Probability distribution',fontsize=18)
	ax.yaxis.set_ticklabels([])


# Load data
###########

pdf_dict = dict()

for key in pdf_path.keys():
	pdf_dict[key] = pd.read_csv(pdf_path[key])


# Plot Fe, FeNi, FeNiSi pdfs at various pressures
#####################################

fig, ((ax1, ax4, ax7),
	(ax2, ax5, ax8),
	(ax3, ax6, ax9)) = plt.subplots(nrows = 3, ncols=3, figsize=(14,8))

minrange = 4.0
maxrange = 5.3
study = 'Fe'
P = 30 # GPa
make_subplot(study,P,ax1,minrange,maxrange,False,False)

study = 'Fe'
P = 53 # GPa
make_subplot(study,P,ax2,minrange,maxrange,False,True)

study = 'Fe'
P = 90 # GPa
make_subplot(study,P,ax3,minrange,maxrange,True,False)

minrange = 3.7
maxrange = 5
study = 'FeNi'
P = 23 # GPa
make_subplot(study,P,ax4,minrange,maxrange,False,False)

study = 'FeNi'
P = 48 # GPa
make_subplot(study,P,ax5,minrange,maxrange,False,False)

study = 'FeNi'
P = 83 # GPa
make_subplot(study,P,ax6,minrange,maxrange,True,False)

minrange = 3.7
maxrange = 5
study = 'FeNiSi'
P = 28 # GPa
make_subplot(study,P,ax7,minrange,maxrange,False,False)

study = 'FeNiSi'
P = 55 # GPa
make_subplot(study,P,ax8,minrange,maxrange,False,False)

study = 'FeNiSi'
P = 86 # GPa
make_subplot(study,P,ax9,minrange,maxrange,True,False)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0.1)

#plt.tight_layout()
fig.savefig('vD_PDF_FeAlloycompare.pdf', format='pdf', transparent=True,
	bbox_inches='tight')
plt.close()


# # Plot Fe, FeNi, FeNiSi pdf comparison
# ######################################

# fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, ncols=1, figsize=(6,10))

# study = 'Fe'
# P = 53 # GPa
# make_subplot(study,P,ax1,4.05,4.8)

# study = 'FeNi'
# P = 41 # GPa
# make_subplot(study,P,ax2,4.05,4.8)

# study = 'FeNiSi'
# P = 41 # GPa
# make_subplot(study,P,ax3,4.05,4.8)

# plt.tight_layout()
# fig.savefig('vD_PDF_compcompare.pdf', format='pdf', transparent=True)
# plt.close()


# # Plot FeNi
# ######################################

# fig, (ax1) = plt.subplots(nrows = 1, ncols=1, figsize=(9.7,3.45))

# study = 'FeNi'
# P = 41 # GPa
# make_subplot(study,P,ax1,4.1,4.5)
# # ax1.xaxis.set_ticklabels(['4.0','4.1','4.2','4.3','4.4','4.5'])
# ax1.xaxis.set_ticks([4,4.1,4.2,4.3,4.4,4.5])
# ax1.tick_params(direction='out',left='off')

# plt.tight_layout()
# fig.savefig('vD_PDF_FeNi.pdf', format='pdf', transparent=True)
# plt.close()


# # Plot FeNi pdfs at various pressures
# #####################################

# fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, ncols=1, figsize=(6,10))

# # study = 'FeNi'
# # P = 23 # GPa
# # make_subplot(study,P,ax1,3.7,5.25)

# study = 'FeNi'
# P = 41 # GPa
# make_subplot(study,P,ax1,4.1,5.3)

# study = 'FeNi'
# P = 63 # GPa
# make_subplot(study,P,ax2,4.1,5.3)

# # study = 'FeNi'
# # P = 83 # GPa
# # make_subplot(study,P,ax3,4.1,4.85)

# study = 'FeNi'
# P = 104 # GPa
# make_subplot(study,P,ax3,4.1,5.3)

# plt.tight_layout()
# fig.savefig('vD_PDF_FeNicompare.pdf', format='pdf', transparent=True)
# plt.close()

# # Plot FeNiSi pdfs at various pressures
# #######################################

# fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, ncols=1, figsize=(6,10))

# study = 'FeNiSi'
# P = 28 # GPa
# make_subplot(study,P,ax1,3.8,4.95)

# study = 'FeNiSi'
# P = 55 # GPa
# make_subplot(study,P,ax2,3.8,4.95)

# study = 'FeNiSi'
# P = 86 # GPa
# make_subplot(study,P,ax3,3.8,4.95)

# plt.tight_layout()
# fig.savefig('vD_PDF_FeNiSicompare.pdf', format='pdf', transparent=True)
# plt.close()


# # Plot Fe, FeNi, FeNiSi pdf comparison for presentation
# ######################################

# def make_subplot(study,P,ax,xmin,xmax):
# 	key = (study,P)
# 	pdf_df = pdf_dict[key]
# 	ymin = min(pdf_df['Prob']) - 0.05*max(pdf_df['Prob'])
# 	ymax = 1.05*max(pdf_df['Prob'])
# 	ax.plot(pdf_df['Bin Midpoint']/1000,pdf_df['Prob'],'.',markersize=8)
# 	ax.plot(vD[key]*np.ones(2), [ymin,ymax], color='red',
# 		linewidth=2)
# 	ax.fill_between([vD[key]-dvD[key],vD[key]+dvD[key]], ymin*np.ones(2), ymax*np.ones(2),
# 		facecolor='red', alpha=0.35)

# 	ax.set_xlim([xmin,xmax])
# 	ax.set_ylim([ymin,ymax])
# 	ax.tick_params(direction='in',left='off')
# 	ax.set_xlabel(r'$v_D$ (km/s)',fontsize=18)
# 	ax.set_ylabel(r'Probability distribution',fontsize=18)
# 	ax.yaxis.set_ticklabels([])

# fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, ncols=1, figsize=(9,9))

# study = 'Fe'
# P = 53 # GPa
# make_subplot(study,P,ax1,4.05,4.8)

# study = 'FeNi'
# P = 41 # GPa
# make_subplot(study,P,ax2,4.05,4.8)

# study = 'FeNiSi'
# P = 41 # GPa
# make_subplot(study,P,ax3,4.05,4.8)

# plt.tight_layout()
# fig.savefig('vD_PDF_compcompare_Pres.pdf', format='pdf', transparent=True)
# plt.close()