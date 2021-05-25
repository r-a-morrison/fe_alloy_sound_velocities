# Are there multiple solutions to vP and vP for a given vD and vphi?

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

# Functions
###########

def vP_eqn(vP):
	return (3/(2*vD**3) - 1/(2*vP**3))**(-2) - (3*vP**2/4 - 3*vphi**2/4)**3

def vS_eqn(vS):
	return (3/(vD**3) - 2/(vS**3))**(-2) - (vphi**2 + 4*vS**2/3)**3

def vPvS_vDeqn(vP,vS,vD):
	return 1/vP**3 + 2/vS**3 - 3/vD**3

def vPvS_vphieqn(vP,vS,vphi):
	return vP**2 - (4/3)*vS**2 - vphi**2

def vP_vS_eqns(initialguess,vD,vphi):
	vP, vS = initialguess
	eqn1 = vPvS_vDeqn(vP,vS,vD)
	eqn2 = vPvS_vphieqn(vP,vS,vphi)
	return eqn1, eqn2



# Set parameter and vP, vS space to explore
###########################################

vD = 4330.0
dvD = 53.0
vphi = 5698.0
dvphi = 25.0
minvP = 6000
maxvP = 8000
minvS = 3000
maxvS = 5000
stepsize = 20

# Build a function that calculates vP, dvP, vS, dvS from vD, dvD, vphi, dvphi

def get_vP_vS(vD,dvD,vphi,dvphi):
	num_samples = 10000
	# Create random normal distributions for vD and vphi
	vD_dist = np.random.normal(vD, dvD, num_samples)
	vphi_dist = np.random.normal(vphi, dvphi, num_samples)
	initialguess = (vD+3000, vD-500)
	# Solve system of equations defined in vP_vS_eqns()
	temp = [fsolve(vP_vS_eqns, initialguess, args = (vD_i,vphi_i))
		for vD_i,vphi_i in zip(vD_dist,vphi_dist)]
	vP_dist = np.array(temp)[:,0]
	vS_dist = np.array(temp)[:,1]

	vP = vP_dist.mean()
	dvP = vP_dist.std()
	vS = vS_dist.mean()
	dvS = vS_dist.std()
	return vP, dvP, vS, dvS

vP, dvP, vS, dvS = get_vP_vS(vD,dvD,vphi,dvphi)
print(vP, dvP, vS, dvS)

# sns.distplot(vP_dist, rug=True)
# plt.savefig('vP_dist.pdf', format='pdf')
# plt.close()

# sns.distplot(vS_dist, rug=True)
# plt.savefig('vS_dist.pdf', format='pdf')
# plt.close()

# sns.jointplot(vP_dist, vS_dist, kind='hex')
# plt.savefig('test_corr.pdf', format='pdf')
# plt.close()

# # Test roots of vP and vS: Solve for vP and vS
# ##############################################
# print('1D plots of vP and vS')
# vParray = np.arange(minvP,maxvP+1,stepsize)
# vSarray = np.arange(minvS,maxvS+1,stepsize)

# # Plot solutions as funtion of vP and vS. Where function crosses zero is a solutions
# fig, (ax0,ax1,ax2,ax3) = plt.subplots(nrows = 4, ncols=1, figsize=(16,16))

# ax0.plot(vParray,vP_eqn(vParray))
# ax0.plot(vParray,np.zeros(len(vParray)))
# ax0.set_xlim([minvP,maxvP])
# #ax0.set_ylim([min(Emin_array),max(Emin_array)])
# ax0.tick_params(direction='in')
# ax0.set_xlabel('$v_P$',fontsize=18)
# #ax0.set_title(title,fontsize=18)

# ax1.plot(vParray,vP_eqn(vParray))
# ax1.plot(vParray,np.zeros(len(vParray)))
# ax1.set_xlim([minvP,maxvP])
# ax1.set_ylim([-1e21,1e21])
# ax1.tick_params(direction='in')
# ax1.set_xlabel('$v_P$',fontsize=18)

# ax2.plot(vSarray,vS_eqn(vSarray))
# ax2.plot(vSarray,np.zeros(len(vSarray)))
# ax2.set_xlim([minvS,maxvS])
# #ax2.set_ylim([min(Emin_array),max(Emin_array)])
# ax2.tick_params(direction='in')
# ax2.set_xlabel('$v_S$',fontsize=18)

# ax3.plot(vSarray,vS_eqn(vSarray))
# ax3.plot(vSarray,np.zeros(len(vSarray)))
# ax3.set_xlim([minvS,maxvS])
# ax3.set_ylim([-1e24,1e24])
# ax3.tick_params(direction='in')
# ax3.set_xlabel('$v_S$',fontsize=18)

# plt.tight_layout()
# fig.savefig('test_vPvS.pdf', format='pdf')
# plt.close()

# # Test roots of vP and vS: Contour plots
# ########################################
# print('Contour plots of vP and vS solutions')

# print('Create contour meshes')
# vPmesh, vSmesh = np.meshgrid(vParray, vSarray)
# vDeqn_contours = [0]
# vphieqn_contours = [0]

# # Plot solutions as funtion of vP and vS. Where function crosses zero is a solutions
# fig, (ax0) = plt.subplots(nrows = 1, ncols=1, figsize=(8,8))

# print('Create contour plot of vD constraint')
# ax0.contourf(vPmesh, vSmesh, vPvS_vDeqn(vPmesh,vSmesh),alpha=0.6,cmap='Blues')
# ax0.set_xlim([minvP,maxvP])
# ax0.set_ylim([minvS,maxvS])
# ax0.set_xlabel('$v_P$',fontsize=18)
# ax0.set_ylabel('$v_S$',fontsize=18)

# hc1 = ax0.contour(vPmesh, vSmesh, vPvS_vDeqn(vPmesh,vSmesh),vDeqn_contours,colors='b')
# ax0.clabel(hc1, inline=1, fontsize=10)
# ax0.set_xlim([minvP,maxvP])
# ax0.set_ylim([minvS,maxvS])
# ax0.set_xlabel('$v_P$',fontsize=18)
# ax0.set_ylabel('$v_S$',fontsize=18)

# print('Create contour plot of vphi constraint')
# ax0.contourf(vPmesh, vSmesh, vPvS_vphieqn(vPmesh,vSmesh),alpha=0.4,cmap='Reds')
# ax0.set_xlim([minvP,maxvP])
# ax0.set_ylim([minvS,maxvS])
# ax0.set_xlabel('$v_P$',fontsize=18)
# ax0.set_ylabel('$v_S$',fontsize=18)

# hc2 = ax0.contour(vPmesh, vSmesh, vPvS_vphieqn(vPmesh,vSmesh),vphieqn_contours,colors='r')
# ax0.clabel(hc2, inline=1, fontsize=10)
# ax0.set_xlim([minvP,maxvP])
# ax0.set_ylim([minvS,maxvS])
# ax0.set_xlabel('$v_P$',fontsize=18)
# ax0.set_ylabel('$v_S$',fontsize=18)

# plt.tight_layout()
# fig.savefig('test_vPvScontour.pdf', format='pdf')
# plt.close()