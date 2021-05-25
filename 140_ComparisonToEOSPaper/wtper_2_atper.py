# Front matter
import pandas as pd
import numpy as np

# Take selection of valid compositions from EOS paper and convert to atomic 
# percentages, because all of the compositions in this study are in terms of 
# atomic percentages

# The compositions to look at:
# Panel a: 5 wt% Ni, ~5 wt% Si, and up to 2% O
# Panel c: 5 wt% Ni, 1 wt% C, 6 wt% Si
# Panel e: 5 wt% Ni, 4 wt% Si, 0.5 wt% C, 3.0 wt% S
# Panel g: 5 wt% Ni and up to 7.5 wt% S

# Element molecular weights
Fe_mw = 55.845
Ni_mw = 58.693
Si_mw = 28.085
O_mw = 15.999
C_mw = 12.01
S_mw = 32.07

# # Convert to a list in atomic fraction (at%/100) composition space
# # a_i = (w_i/M_i) / sum_j(w_j/M_j); where a_i = at frac for element i,
# #                                         w_i = wt frac for element i,
# #                                         M = mol. wt. for element i
# element_mw = np.fromiter(element_mw_dict.values(),
# 	 dtype=float) # Convert from ordereddict to np.array
# # Note: Very important to specificy axis on np.sum(),
# # as this can operate a variety of ways
# atfrac_all = np.divide(np.divide(wtfrac_all,element_mw), 
# 	np.sum(np.divide(wtfrac_all,element_mw),axis=1))

def wtper_2_atper(wtfrac_array,element_mw_array):
	atfrac_array = np.divide(np.divide(wtfrac_array,element_mw_array), 
		np.sum(np.divide(wtfrac_array,element_mw_array),axis=0))
	return atfrac_array

# Panel a: 5 wt% Ni, ~5 wt% Si, and up to 2% O
wtfrac_array = [0.88, 0.05, 0.05, 0.02]
element_mw_array = [Fe_mw, Ni_mw, Si_mw, O_mw]
atfrac_array = wtper_2_atper(wtfrac_array,element_mw_array)
print(atfrac_array)
# [ 0.8023298   0.04337488  0.09064631  0.06364902]

# Panel c: 5 wt% Ni,  6 wt% Si, 1 wt% C
wtfrac_array = [0.88, 0.05, 0.06, 0.01]
element_mw_array = [Fe_mw, Ni_mw, Si_mw, C_mw]
atfrac_array = wtper_2_atper(wtfrac_array,element_mw_array)
# print(atfrac_array)
# [ 0.80484499  0.04351085  0.10911657  0.0425276 ]

# Panel e: 5 wt% Ni, 4 wt% Si, 3.0 wt% S, 0.5 wt% C
wtfrac_array = [0.875, 0.05, 0.04, 0.03, 0.005]
element_mw_array = [Fe_mw, Ni_mw, Si_mw, S_mw, C_mw]
atfrac_array = wtper_2_atper(wtfrac_array,element_mw_array)
# print(atfrac_array)
# [ 0.81198906  0.04414791  0.07380945  0.04847845  0.02157513]

# Panel g: 5 wt% Ni and up to 7.5 wt% S
wtfrac_array = [0.875, 0.05, 0.075]
element_mw_array = [Fe_mw, Ni_mw, S_mw]
atfrac_array = wtper_2_atper(wtfrac_array,element_mw_array)
# print(atfrac_array)
# [ 0.83082121  0.04517181  0.12400698]
