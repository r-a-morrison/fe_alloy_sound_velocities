# Front matter
import os
import glob
import re
import pandas as pd
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


# Rename files for consistency
##############################
# Run one time only

# Find the filepath of all .xy XRD pattern files
filepath_list = [filepath for filepath in glob.glob('*/*.xy')]

# for filepath in filepath_list:
# 	# Get the first part of the label
# 	labelpart1 = re.findall('/([A-Za-z0-9]+)_',filepath)[0]
# 	# Get the second part of the label
# 	labelpart2 = re.findall('/[A-Za-z0-9]+_([A-Za-z0-9]+)_',filepath)[0]
# 	# If they're listed in reverse order, swap for consistency
# 	if re.findall('(Fe)',labelpart2):
# 		prefix = re.findall('([A-Za-z0-9/_.]+)'+labelpart1,filepath)[0]
# 		suffix = re.findall('_'+labelpart2+'([A-Za-z0-9/_.-]+)',filepath)[0]
# 		newfilepath = prefix+labelpart2+'_'+labelpart1+suffix
# 		os.rename(filepath, newfilepath)

# Create folders and move .xy pattern files into folders
########################################################
# Run one time only

# for filepath in filepath_list:
# 	filename = re.findall('/([A-Za-z0-9/_.]+)',filepath)[0]
# 	directory = re.sub('.xy', '', filepath)
# 	if not os.path.exists(directory):
# 		os.makedirs(directory)
# 	if os.path.exists(filepath):
# 		os.rename(filepath, directory+'/'+filename)

# Create overview plots of XRD patterns
#######################################

# Find the filepath of all .xy XRD pattern files
patternfilepath_list = [filepath for filepath in glob.glob('*/*/*.xy')]

for patternfilepath in patternfilepath_list:
	filepath = re.findall('([A-Za-z0-9/_.]+/)[A-Za-z0-9_]+.xy',patternfilepath)[0]
	patternfilename = re.findall('/([A-Za-z0-9_]+.xy)',patternfilepath)[0]
	filename = re.sub('.xy', '', patternfilename)
	pattern_df = pd.read_csv(patternfilepath, delimiter='\s+', comment='#',
                      header = None, engine='python', names=['2Theta','Intensity'])

	print('Plotting pattern for '+filename)
	fig, ax0 = plt.subplots(nrows = 1, ncols=1, figsize=(10, 6))
	ax0.plot(pattern_df['2Theta'],pattern_df['Intensity'],lw=1)

	ax0.xaxis.set_minor_locator(AutoMinorLocator(5))
	ax0.yaxis.set_minor_locator(AutoMinorLocator(5))

	ax0.xaxis.set_ticks_position('both')
	ax0.yaxis.set_ticks_position('both')

	ax0.set_xlabel(r'2$\theta$ ($^\circ$)',fontsize=16)
	ax0.set_ylabel('Intensity',fontsize=16)
	ax0.set_title(filename,fontsize=16)
	fig.savefig(filepath+filename+'_fullpattern.pdf', format='pdf', bbox_inches='tight')
	plt.close()