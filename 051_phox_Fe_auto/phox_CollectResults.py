import warnings
warnings.simplefilter(action = 'ignore', category = UserWarning)

# Front matter
import os
import glob
import re
import pandas as pd
import numpy as np
import scipy.constants as constants

# Find the filepath of all .res NRIXS files
resfilepath_list = [filepath for filepath in glob.glob('*/*.res')]

# Consistency dictionary: For saving space in df
consist_dict = {'ok': 'O', 'acceptable': 'A', 'concerning': 'S'}

# Initialize df structures to store all values from phox .ptl files in
all_fitParam_df        = pd.DataFrame()
all_fitQuality_df      = pd.DataFrame()
all_valsFromData_df    = pd.DataFrame()
all_valsFromRefData_df = pd.DataFrame()
all_valsFromPDOS_df    = pd.DataFrame()


# resfilepath = resfilepath_list[27]
for resfilepath in resfilepath_list:
	filepath = re.findall('([A-Za-z0-9/_]+/)[A-Za-z0-9_]+.res',resfilepath)[0]
	filename = re.findall('/([A-Za-z0-9_]+).res',resfilepath)[0]
	ptlfilepath = filepath+'Output/'+filename+'_phox_ptl.txt'
	psthfilepath = filepath+'Output/'+filename+'_psth_ptl.txt'

	folder = re.findall('([A-Za-z0-9/_]+)/',filepath)[0]
	print(folder)

	# Get date information from directory names
	datetag = re.findall('([A-Za-z0-9]+)_',filepath)[0]
	month = re.findall('[A-Za-z]+',datetag)[0]
	year = re.findall('[0-9]+',datetag)[0]

	# Initialize df structure to store values from phox .ptl file in
	fitParam_df        = pd.DataFrame({'Date': [month+' '+year], 'Folder': [folder], 'Index': [filename]})
	fitQuality_df      = pd.DataFrame({'Date': [month+' '+year], 'Folder': [folder], 'Index': [filename]})
	valsFromData_df    = pd.DataFrame({'Date': [month+' '+year], 'Folder': [folder], 'Index': [filename]})
	valsFromRefData_df = pd.DataFrame({'Date': [month+' '+year], 'Folder': [folder], 'Index': [filename]})
	valsFromPDOS_df    = pd.DataFrame({'Date': [month+' '+year], 'Folder': [folder], 'Index': [filename]})

	# Get NRIXS results from phox results file (_phox_ptl.dat)
	# For version 3.0beta8
	with open(ptlfilepath, 'r') as ptlfile:
		inResults = False
		inElasticPeak = False
		inFromData = False
		inFromRefData = False
		inConsistencyTests = False
		inFromPDOS = False
		normFixed = True
		for line in ptlfile:
			# Switches to determine where to take information from
			if 'Results' in line:
				inResults = True
			if 'Fit of the elastic peak' in line:
				inElasticPeak = True
			if 'Quantities derived directly from the data' in line:
				inFromData = True
				inElasticPeak = False
			if 'Quantities derived after refinement' in line:
				inFromRefData = True
				inFromData = False
			if 'Consistency tests using the refined data' in line:
				inConsistencyTests = True
				inFromRefData = False
			if 'Quantities calculated from the partial DOS' in line:
				inFromPDOS = True
				inConsistencyTests = False
			if ('norm' in line) and inResults and not inConsistencyTests:
				normFixed = False

			# Load data from ptl file into df
			# Fitting parameters
			if 'Temperature of the material' in line:
				temperature = re.findall(':\s+([0-9.]+)',line)[0]
				fitParam_df['T'] = float(temperature)
			if 'Constant background' in line:
				bg = re.findall(':\s+([0-9.]+)',line)[0]
				fitParam_df['Background'] = float(re.sub('D','E',bg))
			if 'Normalization correction' in line:
				norm_initial = re.findall(':\s+([\-0-9.]+)',line)[0]
				norm_initial = float(norm_initial)
			if ('norm' in line) and inResults and not inConsistencyTests:
				norm = re.findall('([\-0-9.D+-]+)\s+[+-]',line)[0]
				fitParam_df['Normalization'] = float(re.sub('D','E',norm))
			if ('background' in line) and inResults and not inElasticPeak:
				bg = re.findall('([0-9.D+-]+)\s+[+-]',line)[0]
				fitParam_df['Background'] = float(re.sub('D','E',bg))
			if 'final' in line:
				chi2 = re.findall(':\s+([0-9.]+)',line)[0]
				fitQuality_df['chi2'] = float(chi2)
				asymmetry = re.findall('\s+([\-0-9.]+)',line)[5]
				fitParam_df['Asymmetry'] = float(asymmetry)

			# Quantities derived directly from the data
			if ('Lamb-Moessbauer factor' in line) and inFromData:
				fLM = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromData_df['fLM'] = float(fLM)
				dfLM = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromData_df['dfLM'] = float(dfLM)
			if ('kinetic energy / atom' in line) and inFromData:
				KE = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromData_df['KE'] = float(KE)
				dKE = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromData_df['dKE'] = float(dKE)
			if ('mean force constant' in line) and inFromData:
				MFC = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromData_df['MFC'] = float(MFC)
				dMFC = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromData_df['dMFC'] = float(dMFC)
			
			# Quantities derived after refinement
			if ('Lamb-Moessbauer factor' in line) and inFromRefData:
				fLM = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromRefData_df['fLM'] = float(fLM)
				dfLM = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromRefData_df['dfLM'] = float(dfLM)
			if ('kinetic energy / atom' in line) and inFromRefData:
				KE = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromRefData_df['KE'] = float(KE)
				dKE = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromRefData_df['dKE'] = float(dKE)
			if ('mean force constant' in line) and inFromRefData:
				MFC = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromRefData_df['MFC'] = float(MFC)
				dMFC = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromRefData_df['dMFC'] = float(dMFC)
			if ('isotope fractionation' in line) and inFromRefData:
				IF = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromRefData_df['IF'] = float(IF)
				dIF = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromRefData_df['dIF'] = float(dIF)
			if ('high T isotope frac.' in line) and inFromRefData:
				highTIF = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromRefData_df['highTIF'] = float(highTIF)
				dhighTIF = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromRefData_df['dhighTIF'] = float(dhighTIF)
			
			# Consistency tests
			if ('detailed balance' in line) and inConsistencyTests:
				detbal = re.findall('[+-]\s*[0-9.]+\s+([0-9.]+)',line)[0]
				fitQuality_df['Det balance'] = float(detbal)
				detbal_stat = re.findall('[+-]\s*[0-9.]+\s+[0-9.]+\s+([a-z]+)',line)[0]
				fitQuality_df['Det bal stat'] = consist_dict[detbal_stat]
			if ('energy/temp. calib.' in line) and inConsistencyTests:
				ETcalib = re.findall('[+-]\s*[0-9.*]+\s+([0-9.]+)',line)[0]
				fitQuality_df['E/T calib'] = float(ETcalib)
				ETcalib_stat = re.findall('[+-]\s*[0-9.*]+\s+[0-9.]+\s+([a-z]+)',line)[0]
				fitQuality_df['E/T calib stat'] = consist_dict[ETcalib_stat]
			if ('negativity of DOS' in line) and inConsistencyTests:
				negDOS = re.findall('[+-]\s*[0-9.]+\s+([0-9.]+)',line)[0]
				fitQuality_df['Negativity of DOS'] = float(negDOS)
				negDOS_stat = re.findall('[+-]\s*[0-9.]+\s+[0-9.]+\s+([a-z]+)',line)[0]
				fitQuality_df['Neg DOS stat'] = consist_dict[negDOS_stat]
			if ('norm of DOS' in line) and inConsistencyTests:
				normDOS = re.findall('[+-]\s*[0-9.]+\s+([0-9.]+)',line)[0]
				fitQuality_df['Norm of DOS'] = float(normDOS)
				normDOS_stat = re.findall('[+-]\s*[0-9.]+\s+[0-9.]+\s+([a-z]+)',line)[0]
				fitQuality_df['Norm DOS stat'] = consist_dict[normDOS_stat]
			if ('Lamb-Moessbauer factor' in line) and inConsistencyTests:
				fLM = re.findall('[+-]\s*[0-9.]+\s+([0-9.]+)',line)[0]
				fitQuality_df['LM factor'] = float(fLM)
				fLM_stat = re.findall('[+-]\s*[0-9.]+\s+[0-9.]+\s+([a-z]+)',line)[0]
				fitQuality_df['fLM stat'] = consist_dict[fLM_stat]
			if ('kinetic energy / atom' in line) and inConsistencyTests:
				KE = re.findall('[+-]\s*[0-9.]+\s+([0-9.]+)',line)[0]
				fitQuality_df['KE'] = float(KE)
				KE_stat = re.findall('[+-]\s*[0-9.]+\s+[0-9.]+\s+([a-z]+)',line)[0]
				fitQuality_df['KE stat'] = consist_dict[KE_stat]
			if ('mean force constant' in line) and inConsistencyTests:
				MFC = re.findall('[+-]\s*[0-9.]+\s+([0-9.]+)',line)[0]
				fitQuality_df['MFC'] = float(MFC)
				MFC_stat = re.findall('[+-]\s*[0-9.]+\s+[0-9.]+\s+([a-z]+)',line)[0]
				fitQuality_df['MFC stat'] = consist_dict[MFC_stat]
			if ('rms average' in line) and inConsistencyTests:
				avg = re.findall('\s+([0-9.]+)',line)[0]
				fitQuality_df['Average'] = float(avg)
				avg_stat = re.findall('[0-9.]+\s+([a-z]+)',line)[0]
				fitQuality_df['Avg stat'] = consist_dict[avg_stat]

			# Quantities calculated from the partial DOS
			if ('Lamb-Moessbauer factor        :' in line) and inFromPDOS:
				fLM = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['fLM'] = float(fLM)
				dfLM = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dfLM'] = float(dfLM)
			if ('kinetic energy                :' in line) and inFromPDOS:
				KE = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['KE'] = float(KE)
				dKE = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dKE'] = float(dKE)
			if ('mean force constant' in line) and inFromPDOS:
				MFC = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['MFC'] = float(MFC)
				dMFC = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dMFC'] = float(dMFC)
			if ('Lamb-Moessbauer factor at T=0' in line) and inFromPDOS:
				fLM0 = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['fLM T=0'] = float(fLM0)
				dfLM0 = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dfLM T=0'] = float(dfLM0)
			if ('kinetic energy         at T=0' in line) and inFromPDOS:
				KE0 = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['KE T=0'] = float(KE0)
				dKE0 = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dKE T=0'] = float(dKE0)
			if ('vibrational specific heat' in line) and inFromPDOS:
				cvib = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['Cvib'] = float(cvib)
				dcvib = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dCvib'] = float(dcvib)
			if ('vibrational entropy' in line) and inFromPDOS:
				Svib = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['Svib'] = float(Svib)
				dSvib = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dSvib'] = float(dSvib)
			if ('resilience' in line) and inFromPDOS:
				resil = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['Resilience'] = float(resil)
				dresil = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dResilience'] = float(dresil)
			if ('Lamb-Moessbauer temperature' in line) and inFromPDOS:
				TLM = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['TLM'] = float(TLM)
				dTLM = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dTLM'] = float(dTLM)
			if ('isotope fractionation' in line) and inFromPDOS:
				IF = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['IF'] = float(IF)
				dIF = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dIF'] = float(dIF)
			if ('high T isotope frac.' in line) and inFromPDOS:
				highTIF = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['highTIF'] = float(highTIF)
				dhighTIF = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dhighTIF'] = float(dhighTIF)

		# If normalization value is missing, replace with initial normalization value
		if normFixed:
			fitParam_df['Normalization'] = norm_initial
		fitParam_df['Norm Fixed'] = normFixed

	# Get NRIXS results from psth results file (_psth_ptl.dat)
	# For version 3.0beta8
	with open(psthfilepath, 'r') as psthfile:
		for line in psthfile:
			if ('free energy' in line):
				Fvib = re.findall('([0-9.]+)\s+[+-]',line)[0]
				valsFromPDOS_df['Fvib'] = float(Fvib)
				dFvib = re.findall('[+-]\s+([0-9.]+)',line)[0]
				valsFromPDOS_df['dFvib'] = float(dFvib)

	all_fitParam_df        = all_fitParam_df.append(fitParam_df)
	all_fitQuality_df      = all_fitQuality_df.append(fitQuality_df)
	all_valsFromData_df    = all_valsFromData_df.append(valsFromData_df)
	all_valsFromRefData_df = all_valsFromRefData_df.append(valsFromRefData_df)
	all_valsFromPDOS_df    = all_valsFromPDOS_df.append(valsFromPDOS_df)
			
all_fitParam_df        = all_fitParam_df[['Date','Index','T','Asymmetry',
						'Background','Normalization','Norm Fixed']].reset_index(drop=True)
all_fitQuality_df      = all_fitQuality_df.reset_index(drop=True)
all_valsFromData_df    = all_valsFromData_df.reset_index(drop=True)
all_valsFromRefData_df = all_valsFromRefData_df.reset_index(drop=True)
all_valsFromPDOS_df    = all_valsFromPDOS_df.reset_index(drop=True)

all_fitParam_df.to_csv('Results/phox_fitparam.csv',index=False)
all_fitQuality_df.to_csv('Results/phox_fitQuality.csv',index=False)
all_valsFromData_df.to_csv('Results/phox_valsFromData.csv',index=False)
all_valsFromRefData_df.to_csv('Results/phox_valsFromRefData.csv',index=False)
all_valsFromPDOS_df.to_csv('Results/phox_valsFromPDOS.csv',index=False)