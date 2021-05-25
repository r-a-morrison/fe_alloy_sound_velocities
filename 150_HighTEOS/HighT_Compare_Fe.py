# Front matter
import datetime
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import time

# Seaborn, useful for graphics
import seaborn as sns

matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

rc = {'lines.linewidth': 1, 
      'axes.labelsize': 18, 
      'axes.titlesize': 18,
      'legend.fontsize': 22,
      'xtick.direction': u'in',
      'ytick.direction': u'in'}
sns.set_style('ticks', rc=rc)


# Functions
###########

def Ks_from_vphi_rho(vphi,rho):
    return vphi**2*rho

def vphi_from_vp_vs(vp,vs):
    return np.sqrt(vp**2-(4/3)*vs**2)


# Define input data
###################

studylist = []
labelchoice = dict()
EOSpath = dict()
colorchoice = dict()
linestylechoice = dict()

# hcp Fe Dewaele 2006 re-analyzed
study = 'hcpFe_Refit'
studylist.append(study)
labelchoice[study] = r'$\gamma$ from Murphy et al. 2011a'
EOSpath[study] = 'MINUTI/hcpFe_Dewaele_Refit_HighT'
colorchoice[study] = 'gray'
linestylechoice[study] = '-'

# hcp Fe Dewaele 2006
study = 'hcpFe_Dewaele'
studylist.append(study)
labelchoice[study] = r'$\gamma$ from this study'
EOSpath[study] = 'MINUTI/hcpFe_Dewaele_HighT'
colorchoice[study] = 'black'
linestylechoice[study] = '--'

# hcp Fe Fei 2016
study = 'hcpFe_Fei'
studylist.append(study)
labelchoice[study] = r'$\gamma$ from Fei et al. 2016'
EOSpath[study] = 'MINUTI/hcpFe_Fei_HighT'
colorchoice[study] = 'green'
linestylechoice[study] = '-.'


# Import data
#############

EOS_dfdict = dict()
temperatures = [5500]

for study in studylist:
    print('Now importing '+study)
    for T in temperatures:
        
        density_path = EOSpath[study]+'/'+str(T)+'K_q_rep/'+study+'_dns.dat'
        rho_df = pd.read_csv(density_path,header=None,sep='\s+',comment='#',engine='python')
        rho_df.columns = ['P','rho']

        vphi_path = EOSpath[study]+'/'+str(T)+'K_q_rep/'+study+'_smv.dat'
        vphi_df = pd.read_csv(vphi_path,header=None,sep='\s+',comment='#',engine='python')
        vphi_df.columns = ['P','vphi']

        EOS_df = pd.merge(rho_df,vphi_df,on='P')

        EOS_df['Ks'] = Ks_from_vphi_rho(EOS_df['vphi'],EOS_df['rho'])

        EOS_dfdict[(study,T)] = EOS_df


# Load Seismic Observations
###########################

# Load Earth values and do some basic calcs so we're working with complete data sets
PREM_innercore = pd.read_csv('../../FeAlloyEOS/PREM/PREM_innercore.csv')
PREM_innercore['Ks'] = Ks_from_vphi_rho(PREM_innercore['vphi'],PREM_innercore['rho'])

AK135_innercore = pd.read_csv('../../FeAlloyEOS/AK135/AK135F_innercore.csv')
AK135_innercore['Ks'] = Ks_from_vphi_rho(AK135_innercore['vphi'],AK135_innercore['rho'])
# AK135 doesn't report pressure, so we'll assume the same pressure/depth relation as in PREM
PressureEst = interp1d(PREM_innercore['depth'], PREM_innercore['P'])
AK135_innercore['P'] = PressureEst(AK135_innercore['depth'])

# Attanayaka et al. 2014 measured Vp as a function of depth and latitude
# Use AK135-F Vs and Attanayaka Vp from bin 5 (Western hemisphere)
Attanayaka = pd.read_csv('../../FeAlloyEOS/Attanayaka/Attanayaka_Vp.csv')
# Attanayaka doesn't report Vs, so we'll assume the same Vs/depth relation as in AK135
VsEst = interp1d(AK135_innercore['depth'], AK135_innercore['vs'])
Attanayaka['vs'] = VsEst(Attanayaka['depth'])
# Attanayaka doesn't report pressure, so we'll assume the same pressure/depth relation as in PREM
Attanayaka['P'] = PressureEst(Attanayaka['depth'])
Attanayaka['BIN 6 vphi'] = vphi_from_vp_vs(Attanayaka['BIN 6 vp'],Attanayaka['vs'])
Attanayaka['BIN 3 vphi'] = vphi_from_vp_vs(Attanayaka['BIN 3 vp'],Attanayaka['vs'])
# Attanayaka doesn't report density, so we'll assume the same density/depth relation as in AK135
rhoEst = interp1d(AK135_innercore['depth'], AK135_innercore['rho'])
Attanayaka['rho'] = rhoEst(Attanayaka['depth'])
Attanayaka['BIN 6 Ks'] = Attanayaka['BIN 6 vphi']**2*Attanayaka['rho']
Attanayaka['BIN 3 Ks'] = Attanayaka['BIN 3 vphi']**2*Attanayaka['rho']

# Define uncertainty on inner core density, Ks, and vphi
rho_sigma = 0.02
Ks_sigma = 0.023
vphi_sigma = 0.0056

# Get the high and low bounds on PREM
PREM_innercore['rho max']  = PREM_innercore['rho']*(1+rho_sigma)
PREM_innercore['Ks max']   = PREM_innercore['Ks']*(1+Ks_sigma)
PREM_innercore['vphi max'] = PREM_innercore['vphi']*(1+vphi_sigma)
PREM_innercore['rho min']  = PREM_innercore['rho']*(1-rho_sigma)
PREM_innercore['Ks min']   = PREM_innercore['Ks']*(1-Ks_sigma)
PREM_innercore['vphi min'] = PREM_innercore['vphi']*(1-vphi_sigma)

# Get the high and low bounds on AK135
AK135_innercore['rho max']  = AK135_innercore['rho']*(1+rho_sigma)
AK135_innercore['Ks max']   = AK135_innercore['Ks']*(1+Ks_sigma)
AK135_innercore['vphi max'] = AK135_innercore['vphi']*(1+vphi_sigma)
AK135_innercore['rho min']  = AK135_innercore['rho']*(1-rho_sigma)
AK135_innercore['Ks min']   = AK135_innercore['Ks']*(1-Ks_sigma)
AK135_innercore['vphi min'] = AK135_innercore['vphi']*(1-vphi_sigma)


# Create Plots
##############

fig, (ax0, ax1, ax2) = plt.subplots(nrows = 3, ncols=1, figsize=(4, 8), sharex=True)

# Density
htemp, = ax0.plot(AK135_innercore['P'],AK135_innercore['rho'],'-',
                   label='AK135-F',color='red',lw=1.5)
ax0.fill_between(AK135_innercore['P'], AK135_innercore['rho max'],
                     AK135_innercore['rho min'],facecolor='red', alpha=0.3)
h2 = [htemp]

for study in studylist:
    for T in temperatures:
	    EOS_df      = EOS_dfdict[(study,T)]
	    
	    htemp, = ax0.plot(EOS_df['P'],EOS_df['rho'],'-',
	             label=labelchoice[study],color=colorchoice[study],lw=1,
	             linestyle=linestylechoice[study])
	    # ax0.fill_between(EOS_df['P'], EOS_df['rho']+0.032,
	    # 	     EOS_df['rho']-0.032, facecolor=colorchoice[study], alpha=0.3)
    h2.append(htemp)

ax0.set_ylabel(r'$\rho$ (g/cm$^3$)', fontsize=14)
ax0.set_ylim(12.4,14.3)
ax0.tick_params(direction='in',right='on',top='on')


# Isentropic Bulk Modulus
ax1.plot(AK135_innercore['P'],AK135_innercore['Ks'],'-',
                   label='AK-135F',color='red',lw=1.5)
ax1.fill_between(AK135_innercore['P'], AK135_innercore['Ks max'],
                     AK135_innercore['Ks min'],facecolor='red', alpha=0.3)
# htemp, = ax1.plot(Attanayaka['P'],Attanayaka['BIN 3 Ks'],':',dash_capstyle='round',
#                    label='AK135 + Attanayake 2014 East (Bin 3)',color='black',lw=1.5)
# h2.append(htemp)
# htemp, = ax1.plot(Attanayaka['P'],Attanayaka['BIN 6 Ks'],'-',
#                    label='AK135 + Attanayake 2014 West (Bin 6)',color='black',lw=1.5)
# h2.append(htemp)

for study in studylist:
    for T in temperatures:
	    EOS_df      = EOS_dfdict[(study,T)]
	    
	    ax1.plot(EOS_df['P'],EOS_df['Ks'],'-',
	             label=labelchoice[study],color=colorchoice[study],lw=1,
	             linestyle=linestylechoice[study])
	    # ax1.fill_between(EOS_df['P'], EOS_df['Ks']+17,
	    # 	     EOS_df['Ks']-17, facecolor=colorchoice[study], alpha=0.3)

ax1.set_ylabel(r'$K_S$ (GPa)', fontsize=14)
ax1.set_ylim(1270,1500)
ax1.tick_params(direction='in',right='on',top='on')


# Bulk Sound Speed
ax2.plot(AK135_innercore['P'],AK135_innercore['vphi'],'-',
                   label='AK-135F',color='red',lw=1.5)
ax2.fill_between(AK135_innercore['P'], AK135_innercore['vphi max'],
                     AK135_innercore['vphi min'],facecolor='red', alpha=0.3)

# ax2.plot(Attanayaka['P'],Attanayaka['BIN 3 vphi'],':',dash_capstyle='round',
#                    label='AK135',color='black',lw=1.5)   
# ax2.plot(Attanayaka['P'],Attanayaka['BIN 6 vphi'],'-',
#                    label='AK135 + Attanayake 2014 East (Bin 3)',color='black',lw=1.5)

for study in studylist:
    for T in temperatures:
	    EOS_df      = EOS_dfdict[(study,T)]
	    
	    ax2.plot(EOS_df['P'],EOS_df['vphi'],'-',
	             label=labelchoice[study],color=colorchoice[study],lw=1,
	             linestyle=linestylechoice[study])
	    # ax2.fill_between(EOS_df['P'], EOS_df['vphi']+0.071,
	    # 	     EOS_df['vphi']-0.071, facecolor=colorchoice[study], alpha=0.3)

ax2.set_xlabel(r'Pressure (GPa)', fontsize=14)
ax2.set_ylabel(r'$v_\phi$ (km/s)', fontsize=14)
ax2.set_xlim(329.1,364)
ax2.set_ylim(9.7,10.6)
ax2.tick_params(direction='in',right='on',top='on')

# ax0.legend(fontsize=14,bbox_to_anchor=(1.6, 1.45),handles = hlabels)
ax1.legend(fontsize=12,bbox_to_anchor=(0.85, 2.6),handles = h2)

fig.subplots_adjust(hspace=0.05)

fig.savefig('Plots/HighT_EOS_Plot_compare.pdf', format='pdf', bbox_inches='tight')