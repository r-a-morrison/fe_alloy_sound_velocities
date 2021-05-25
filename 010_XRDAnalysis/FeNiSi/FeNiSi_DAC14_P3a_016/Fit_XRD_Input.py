# What file are we fitting?
file_name = 'FeNiSi_DAC14_P3a_016'

# What NRIXS experiment does this correspond to?
NRIXSexp = 'FeNiSi_DAC14_P3'

# How many background parameters do we need to fit the background?
bg_par_num = 3

# What range should we fit data in?
range_2theta = [24,26]
primaryphase = 'bcc FeNiSi'

# What phases are present, and what is their structure?
phasedict = {'bcc FeNiSi': 'bcc'}

# What starting peak intensities should we use for each phase? (in counts)
peakintdict = {'bcc FeNiSi': 100} # Rough guess

# How many peaks for each phase?
peaknumdict = {'bcc FeNiSi': 1}

# What starting lattice parameters should we use for each phase?
# hcp should have a and c while bcc should have just a
latpardict = {'bcc FeNiSi': {'a':2.863}}
#			{'hcp FeNi': {'a':2.460,'c':3.970}}

# What starting peak width should we use for each phase? (in 2theta degrees)
peakwidthdict = {'bcc FeNiSi': 0.05}

# What starting peak shape should we use for each phase? (from 0 to 1, where 0 is
# pure Gaussian and 1 is pure Lorentzian)
peakshapedict = {'bcc FeNiSi': 0.5}

# What was the x-ray wavelenth?
wavelength = 0.860255