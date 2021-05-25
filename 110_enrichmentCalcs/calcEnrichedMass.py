# Calculates enriched and natural isotopic molecular masses in g/mol

# Front matter
##############
import re
import time
import pandas as pd
import numpy as np
from scipy import constants

start_time = time.time()


# Define list of compositions to calculate molecular mass of
############################################################

# Enter values in atomic percentage (not wt%, etc.)

comp_list = ['Fe1.00',
			 'Fe0.91Ni0.09',
			 'Fe0.80Ni0.10Si0.10',
			 'Fe0.92Ni0.08',
			 'Fe0.85Si0.15',
			 'Fe0.868Ni0.086Si0.046',
			 'Fe0.89Ni0.04Si0.07',
			 'Fe0.91Si0.09',
			 'Fe0.89Si0.11']


# Define relative atomic masses of natural and enriched isotopes
################################################################

am_Fe_nat = 55.845 # natural molecular mass of Fe
am_Fe_enr = 56.942 # 95% enriched in 57Fe

am_nat = dict()
am_nat['Fe'] = am_Fe_nat
am_nat['Ni'] = 58.693
am_nat['Si'] = 28.085
am_nat['O']  = 15.999
am_nat['S']  = 32.06
am_nat['C']  = 12.011
am_nat['H']  =  1.008

am_enr = am_nat.copy()
am_enr['Fe'] = am_Fe_enr


# Functions
###########

def calcMolMass(comp_string, atomic_mass):
	elements = re.findall('([A-Za-z]+)[0-9.]+',comp_string)
	percentages = re.findall('[A-Za-z]+([0-9.]+)',comp_string)
	percentages = [float(percentage) for percentage in percentages]
	M = sum(atomic_mass[element]*percentage for element, percentage 
		in zip(elements,percentages))
	return M


# Calculate natural and enriched molecular masses
#################################################

M_nat_list = []
M_enr_list = []

for comp_string in comp_list:
	M_nat = calcMolMass(comp_string, am_nat)
	M_enr = calcMolMass(comp_string, am_enr)

	print(comp_string+':')
	print('\tNatural enrichment  '+str(round(M_nat, 3))+' g/mol')
	print('\tEnriched in 57Fe    '+str(round(M_enr, 3))+' g/mol')

	# Store in lists
	M_nat_list.append(M_nat)
	M_enr_list.append(M_enr)

# Save results to a csv file
M_results = pd.DataFrame({'Composition':comp_list,
						  'Natural':M_nat_list,
						  '57Fe Enriched':M_enr_list})
M_results = M_results.reindex_axis(['Composition','Natural',
						  '57Fe Enriched'],axis=1)
M_results = M_results.round(decimals=3)
M_results.to_csv('MolecularMasses.csv',index=False)