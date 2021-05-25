
# Written by Rachel Morrison on July 1, 2016
# Updated October 2017

# ## Input Data and Parameters

from Fit_XRD_Input import *

# ## Outline

# - Import data
# - Specify 2theta range
# - Identify phases
# - Set starting parameter values for a (and c)
# - Identify peaks present in 2theta range for given a (and c)
# - Get starting parameters for background
# - Create model
# - Fit peak intensitites, peak widths, peak shapes, lattice parameters, and background

# ## Front matter

import warnings
warnings.simplefilter(action = 'ignore', category = FutureWarning)
warnings.simplefilter(action = 'ignore', category = UserWarning)

# Front matter
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import constants
from scipy.fftpack import rfft, irfft, fftfreq
import scipy
import time
import itertools
import os
import sys


# Seaborn, useful for graphics
import seaborn as sns

# Magic function to make matplotlib inline; other style specs must come AFTER
# %matplotlib inline

matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20) 

rc = {'lines.linewidth': 1, 
      'axes.labelsize': 20, 
      'axes.titlesize': 20,
      'legend.fontsize': 26,
      'xtick.direction': u'in',
      'ytick.direction': u'in'}
sns.set_style('ticks', rc=rc)

pd.options.mode.chained_assignment = None  # default='warn'


# Constants
N_A = constants.value('Avogadro constant')


# ## Functions

# Import pattern
def import_pattern(file_path):
    pattern = pd.read_csv(file_path, delimiter='\s+', comment='#',
                      header = None, engine='python')
    pattern.columns = ('2theta', 'Intensity')
    return pattern

# Estimate the baseline
def get_baseline(y, lam, p, niter=10):
    L = len(y)
    D = scipy.sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = scipy.sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = scipy.sparse.linalg.spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

# Set pattern range
def set_pattern_range(pattern,pattern_range):
    cropped_pattern = pattern[(pattern['2theta']>=pattern_range[0]) &
                              (pattern['2theta']<=pattern_range[1])]
    return cropped_pattern

# Polynomial function with an arbitrary number of terms
def bg_func(x,*params):
    return sum([params*(x**i) for i, params in enumerate(params)])

# Determine bcc hkl positions
def get_bcc_hkl(max_hkl,max_peaks):
    possible_hkl = pd.DataFrame(list(itertools.product(list(range(0,max_hkl)),repeat = 3)),
                               columns=['h','k','l'])
    possible_hkl['h+k+l'] = possible_hkl['h']+possible_hkl['k']+possible_hkl['l']
    possible_hkl = possible_hkl[possible_hkl['h+k+l']%2 == 0]
    possible_hkl = possible_hkl.drop(0)
    possible_hkl['h2+k2+l2'] = (possible_hkl['h']**2+possible_hkl['k']**2+
                                possible_hkl['l']**2)
    possible_hkl = possible_hkl.sort(['h2+k2+l2','h','k','l'],
                                         ascending=[1, 0, 0, 0]).reset_index(drop=True)
    possible_hkl = possible_hkl.drop_duplicates(subset=['h2+k2+l2'],
                                                    keep='first').reset_index(drop=True)
    possible_hkl = possible_hkl.ix[:max_peaks-1]
    return np.array(possible_hkl[['h','k','l']])

# Determine hcp hkl positions
def get_hcp_hkl(max_hkl,max_peaks):
    possible_hkl = pd.DataFrame(list(itertools.product(list(range(0,max_hkl)),repeat = 3)),
                               columns=['h','k','l'])
    possible_hkl['h+2k'] = possible_hkl['h']+2*possible_hkl['k']
    # Drop the systematically absent peaks
    possible_hkl = possible_hkl.drop(possible_hkl[(possible_hkl['h+2k']%3==0)
                                                  &(possible_hkl['l']%2==1)].index)
    # Drop (0,0,0)
    possible_hkl = possible_hkl.drop(0)
    # Assumes a c/a ratio of 1.6 for choosing peak order (approx. correct)
    possible_hkl['a2/d2'] = (4/3)*(possible_hkl['h']**2+possible_hkl['h']*possible_hkl['k']+
                                          possible_hkl['k']**2)+possible_hkl['l']**2/(1.6)**2
    possible_hkl = possible_hkl.sort(['a2/d2','h','k','l'],
                                         ascending=[1, 0, 0, 0]).reset_index(drop=True)
    possible_hkl = possible_hkl.drop_duplicates(subset=['a2/d2'],
                                                    keep='first').reset_index(drop=True)
    possible_hkl = possible_hkl.ix[:max_peaks-1]
    return np.array(possible_hkl[['h','k','l']])

# Function for a Lorentzian peak with an amplitude of 1 (not normalized)
def lorentzian(x,peakwidth,peakloc):
    m = (1/2)*peakwidth/((x-peakloc)**2+((1/2)*peakwidth)**2)
    return m

# Function for a Gaussian peak with an amplitude of 1 (not normalized)
def gaussian(x,peakwidth,peakloc):
    m = np.exp(-(x-peakloc)**2/(2*peakwidth**2))
    return m

# Function for a pseudovoigt peak
def pseudovoigt(x,peakshape,peakintensity,peakwidth,peakloc):
    m = peakintensity*(peakshape*lorentzian(x,peakwidth,peakloc) +
         (1-peakshape)*gaussian(x,peakwidth,peakloc))
    return m

# Convert a d-spacing value to 2theta location
def d_to_2theta(d,wavelength):
    theta_in_rad = np.arcsin(wavelength/(2*d))
    theta_in_deg = theta_in_rad*180/np.pi
    return 2*theta_in_deg

# Function for XRD peaks for a single phase
def peaks(x,struct,peaknum,peakshape,peakintensities,peakwidth,latpar,wavelength):
    if struct == 'bcc':
        a = latpar
        all_hkl = pd.DataFrame(get_bcc_hkl(15,peaknum),columns=['h','k','l'])
        all_d = a/np.sqrt(all_hkl['h']**2+all_hkl['k']**2+all_hkl['l']**2)
    if struct == 'hcp':
        a = latpar[0]
        c = latpar[1]
        all_hkl = pd.DataFrame(get_hcp_hkl(15,peaknum),columns=['h','k','l'])
        all_d = ((4/3)*(all_hkl['h']**2+all_hkl['h']*all_hkl['k']+all_hkl['k']**2)/a**2 + 
                 all_hkl['l']**2/c**2)**(-1/2)       
    m = 0
    for i,d in enumerate(all_d):
        peakintensity = peakintensities[i]
        peakloc = d_to_2theta(d,wavelength)
        peak = pseudovoigt(x,peakshape,peakintensity,peakwidth,peakloc)
        m = m + peak
    return m

# Function for a combination of XRD peaks
def allpeaks(x,phasedict,peaknumdict,allpeakintensities,peakwidths,peakshapes,alllatpar,
            wavelength):
    peakcount = 0
    latparcount = 0
    phasepeaklist = []
    for i,phase in enumerate(list(phasedict.keys())):
        struct = phasedict[phase]
        peaknum = peaknumdict[phase]
        peakintensities = allpeakintensities[peakcount:peakcount+peaknum]
        peakwidth = peakwidths[i]
        peakshape = peakshapes[i]
        if struct == 'bcc':
            latpar = alllatpar[latparcount:latparcount+1]
            peaks_in_phase = peaks(x,struct,peaknum,peakshape,peakintensities,peakwidth,
                                   latpar,wavelength)
            latparcount = latparcount + 1
        if struct == 'hcp':
            latpar = alllatpar[latparcount:latparcount+2]
            peaks_in_phase = peaks(x,struct,peaknum,peakshape,peakintensities,peakwidth,
                                   latpar,wavelength)
            latparcount = latparcount + 2
        phasepeaklist.append(peaks_in_phase)
        
        peakcount = peakcount + peaknum
        
    m = sum(phasepeaklist)
    return m


# Make the model of XRD patterns + background
def make_model(phasedict,peaknumdict,
              peakintensities0, peakwidths0, peakshapes0, latpar0, bgpar0, wavelength):
    len1 = len(peakintensities0)
    len2 = len(peakwidths0)
    len3 = len(peakshapes0)
    len4 = len(latpar0)
    
    def model(x,*params):
        peakintensitypar = params[:len1]
        peakwidthpar = params[len1:(len1+len2)]
        peakshapepar = params[(len1+len2):(len1+len2+len3)]
        latpar = params[(len1+len2+len3):(len1+len2+len3+len4)]
        bgpar = params[(len1+len2+len3+len4):]
        
        # Multicomponent parameters MUST have a star in front of them
        m = bg_func(x,*bgpar) + allpeaks(x,phasedict,peaknumdict,peakintensitypar,
                                         peakwidthpar,peakshapepar,latpar,wavelength)
        return m
    return model

# Determine the initial parameters needed
def get_par0(peakintensities0,peakwidths0,peakshapes0,latpar0,bgpar0,
               peakintensitynames, peakwidthnames, peakshapenames, latparnames, bgparnames):
    intstartbnds = list(np.zeros(len(peakintensities0)))
    intendbnds = list(np.ones(len(peakintensities0))*99999999)
    wstartbnds = list(np.zeros(len(peakwidths0)))
    wendbnds = list(np.ones(len(peakwidths0))*2) # in degrees
    shapestartbnds = list(np.zeros(len(peakshapes0)))
    shapeendbnds = list(np.ones(len(peakshapes0)))
    latstartbnds = list(np.ones(len(latpar0))*0.1)
    latendbnds = list(np.ones(len(latpar0))*10)
    bgstartbnds = list(np.ones(len(bgpar0))*(-99999999))
    bgendbdns = list(np.ones(len(bgpar0))*99999999)
    
    par0 = list(itertools.chain(peakintensities0,peakwidths0,peakshapes0,latpar0,bgpar0))
    parnames = list(itertools.chain(peakintensitynames,peakwidthnames,peakshapenames,
                                    latparnames,bgparnames))
    startbounds = list(itertools.chain(intstartbnds,wstartbnds,shapestartbnds,
                                       latstartbnds,bgstartbnds))
    endbounds = list(itertools.chain(intendbnds,wendbnds,shapeendbnds,latendbnds,
                                         bgendbdns))
    parbounds = [startbounds,endbounds]
    return par0,parnames,parbounds


# ## Body

# Determine how many and what parameters we need to fit

# What phases are present?
phases = list(phasedict.keys())

# We need our starting parameters in list form. Create a list of names too to keep
# things organized.

# We can start the peak intensities at 100
peakintensities0 = []
peakintensitynames = []
for phase in phases:
    peakintensities0.extend(peakintdict[phase]*np.ones(peaknumdict[phase]))
    for i in range(0,peaknumdict[phase]):
        peakintensitynames.append(phase+' Peak '+str(i)+' Intensity')

# The peak widths are user specified
peakwidths0 = []
peakwidthnames = []
for phase in phases:
    peakwidths0.append(peakwidthdict[phase])
    peakwidthnames.append(phase + ' Peak Width')

# The peak shapes are user specified
peakshapes0 = []
peakshapenames = []
for phase in phases:
    peakshapes0.append(peakshapedict[phase])
    peakshapenames.append(phase + ' Peak Shape')

# The lattice parameters are user specified
latpar0 = []
latparnames = []
for phase in phases:
    latpars = latpardict[phase]
    if phasedict[phase] == 'bcc':
        latpar0.append(latpars['a'])
        latparnames.append(phase+' a')
    if phasedict[phase] == 'hcp':
        latpar0.append(latpars['a'])
        latpar0.append(latpars['c'])
        latparnames.append(phase+' a')
        latparnames.append(phase+' c')

# The background starting parameters will be found
bgparnames = []
for i in range(bg_par_num):
    bgparnames.append('Background '+str(i))


# Start Fitting
print('Fitting '+file_name+'.xy')
sys.stdout.flush()
start = time.time()

file_path = file_name+'.xy'

# Import data in a Pandas data frame
pattern = import_pattern(file_path)

# Estimate the baseline of the pattern and save to data frame
baseline = get_baseline(pattern['Intensity'],10**6,0.01)
pattern['Baseline Est'] = baseline

# Set range and save to new dataframe "data"
data = set_pattern_range(pattern,range_2theta)

# Fit the background with a polynomial function and save to data frame

# Default bg starting parameters
bg_par_start = np.ones(bg_par_num)

# The initial guess p0 here is necessary to indicate the number of parameters to use
bgpar0,bgcov0 = scipy.optimize.curve_fit(bg_func, data['2theta'], data['Baseline Est'],
         p0=bg_par_start)

data['Background'] = bg_func(data['2theta'],*bgpar0)



# Fit the data with the model
par0,parnames,parbounds = get_par0(peakintensities0,peakwidths0,peakshapes0,latpar0,bgpar0,
               peakintensitynames, peakwidthnames, peakshapenames, latparnames, bgparnames)

m = make_model(phasedict,peaknumdict,
              peakintensities0, peakwidths0, peakshapes0, latpar0, bgpar0, wavelength)
par,cov = scipy.optimize.curve_fit(m, data['2theta'], data['Intensity'],
                                   sigma = np.sqrt(data['Intensity']),
                                   p0=par0, bounds=parbounds)




# Break up parameters into individual groups
len1 = len(peakintensities0)
len2 = len(peakwidths0)
len3 = len(peakshapes0)
len4 = len(latpar0)

peakintensitypar = par[:len1]
peakwidthpar = par[len1:(len1+len2)]
peakshapepar = par[(len1+len2):(len1+len2+len3)]
latpar = par[(len1+len2+len3):(len1+len2+len3+len4)]
bgpar = par[(len1+len2+len3+len4):]



# Create a data frame with data and fit information for plotting
results_curves = data[['2theta','Intensity',]]
results_curves['Model'] = m(data['2theta'],*par)
results_curves['Residual'] = results_curves['Intensity'] - results_curves['Model']
results_curves['Background'] = bg_func(data['2theta'],*bgpar)



# Calculate fit for each phase
peakcount = 0
latparcount = 0
phasepeaklist = [] 
for i,phase in enumerate(phases):
    struct = phasedict[phase]
    peaknum = peaknumdict[phase]
    peakintensities = peakintensitypar[peakcount:peakcount+peaknum]
    peakwidth = peakwidthpar[i]
    peakshape = peakshapepar[i]
    if struct == 'bcc':
        lp = latpar[latparcount:latparcount+1]
        latparcount = latparcount + 1
    if struct == 'hcp':
        lp = latpar[latparcount:latparcount+2]
        latparcount = latparcount + 2
    peakcount = peakcount + peaknum

    results_curves[phase] = peaks(data['2theta'],struct,peaknum,peakshape,
                   peakintensities,peakwidth,lp,wavelength)

    
# Create a dataframe with covariance information
cov_df = pd.DataFrame(cov,columns = parnames,index = parnames)

# Create a dictionary with parameter information
par_vals = dict(zip(parnames,par))
par_std = dict(zip(parnames,np.sqrt(cov.diagonal())))

# Create a dataframe with parameter information
par_df = pd.DataFrame(parnames, columns=['Parameters'])
par_df['Values'] = par
par_df['Std'] = np.sqrt(cov.diagonal())

# Determine weighted profile R-factor
wR = np.sqrt(sum((1/data['Intensity'])*(data['Intensity'] -
        m(data['2theta'],*par))**2))/np.sqrt(sum(data['Intensity']))
# Determine the expected weighted profile R-factor
expR = np.sqrt(len(data['Intensity'])-len(par))/np.sqrt(sum(data['Intensity']))
# Determine the reduced chi^2 (GOF)
chi2 = np.sqrt(sum((1/data['Intensity'])*(data['Intensity'] -
        m(data['2theta'],*par))**2))/np.sqrt(len(data['Intensity'])-len(par))

# Create a dataframe with results
results = pd.DataFrame([[file_name,wR,chi2]],columns = ['File name','wR','Red. chi^2'])
results['Start'] = range_2theta[0]
results['End'] = range_2theta[1]
results['Wavelength'] = wavelength
results['Bg num'] = bg_par_num
for latparname in latparnames:
    results[latparname] = par_vals[latparname]
#     results[latparname + ' std'] = par_std[latparname]
for phase in phases:
    if phasedict[phase] == 'hcp':
        a = par_vals[phase+' a']
        c = par_vals[phase+' c']
        results[phase + ' Vol'] = (np.sqrt(3)/2)*a**2*c
        results[phase + ' c/a'] = c/a
#         a_std = par_std[phase+' a']
#         c_std = par_std[phase+' c']
#         results[phase + ' Vol Std'] = (np.sqrt(3)/2)*np.sqrt((2*a*c*a_std)**2+(a**2*c_std)**2)
    if phasedict[phase] == 'bcc':
        a = par_vals[phase+' a']
        results[phase + ' Vol'] = a**3

# Save data
dir_name = '../'+file_name
if not os.path.isdir(dir_name):
    # Save stuff
    os.makedirs(dir_name)
    results.to_csv(dir_name+'/'+file_name+'_results.csv',index=False)
    par_df.to_csv(dir_name+'/'+file_name+'_params.csv',index=False)
    results_curves.to_csv(dir_name+'/'+file_name+'_curves.csv',index=False)
    cov_df.to_csv(dir_name+'/'+file_name+'_cov.csv')
else:
#         ext = 2
#         while os.path.isdir(dir_name+'_'+str(ext)):
#             ext += 1
#         # Save stuff
#         dir_name = dir_name+'_'+str(ext)
#         os.makedirs(dir_name)
    results.to_csv(dir_name+'/'+file_name+'_results.csv',index=False)
    par_df.to_csv(dir_name+'/'+file_name+'_params.csv',index=False)
    results_curves.to_csv(dir_name+'/'+file_name+'_curves.csv',index=False)
    cov_df.to_csv(dir_name+'/'+file_name+'_cov.csv')


# Calculate the correlation
D = np.diag(np.diag(cov)**(-1/2))
corr = D.dot(cov.dot(D))
# Plot and save the correlation matrix
fig = plt.figure(figsize=(11, 9))
plt.pcolor(corr,cmap='seismic', vmin=-1, vmax=1)
plt.colorbar()
plt.yticks(np.arange(0.5,len(corr)+0.5),parnames)
plt.xticks(np.arange(0.5,len(corr)+0.5),parnames,rotation='vertical');
plt.xlim([0,len(corr)])
plt.ylim([0,len(corr)])
plt.tick_params(axis='both', which='major', labelsize=20*23/len(par))
plt.savefig(dir_name+'/'+file_name+'_coor.pdf',format = 'pdf',bbox_inches='tight')
plt.close(fig)


# Plot and save the data and model comparison
fig = plt.figure(figsize=(10, 8))
plt.plot(results_curves['2theta'], results_curves['Intensity'],'black',label = 'Data')
plt.plot(results_curves['2theta'], results_curves['Background'],"#e46a2e",
         label = 'Background')
phasecolors = ["#00cdff","#e63a85","#93d147","#b45bdc"]
for i,phase in enumerate(phases):
    plt.plot(results_curves['2theta'], results_curves[phase]+results_curves['Background'],
             color = phasecolors[i],label = phase)
plt.plot(results_curves['2theta'], results_curves['Residual'],'darkgray',label = 'Residual')
plt.xlabel(r'2$\theta$',fontsize=20)
plt.ylabel('Counts',fontsize=20)
plt.title(file_name,fontsize=20)
plt.legend(fontsize=18)
plt.xlim(range_2theta);
plt.savefig(dir_name+'/'+file_name+'_plot.pdf',format = 'pdf',bbox_inches='tight')
plt.close(fig)


end = time.time()
print(file_name+': chi^2 = %.2f' % chi2)
print('time elapsed: '+str(end-start)+' s')
sys.stdout.flush()

print('Completed')