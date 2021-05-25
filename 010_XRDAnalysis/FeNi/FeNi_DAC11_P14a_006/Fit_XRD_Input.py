# What file are we fitting?
file_name = 'FeNi_DAC11_P14a_006'

# What NRIXS experiment does this correspond to?
NRIXSexp = 'FeNi_DAC11_P14'

# How many background parameters do we need to fit the background?
bg_par_num = 4

# What range should we fit data in?
range_2theta = [23,29.5]
primaryphase = 'hcp FeNi'

# What phases are present, and what is their structure?
phasedict = {'hcp FeNi': 'hcp'}

# What starting peak intensities should we use for each phase? (in counts)
peakintdict = {'hcp FeNi': 100} # Rough guess

# How many peaks for each phase?
peaknumdict = {'hcp FeNi': 3}

# What starting lattice parameters should we use for each phase?
# hcp should have a and c while bcc should have just a
latpardict = {'hcp FeNi': {'a':2.30,'c':3.64}}
#			{'hcp FeNi': {'a':2.460,'c':3.970}}

# What starting peak width should we use for each phase? (in 2theta degrees)
peakwidthdict = {'hcp FeNi': 0.05}

# What starting peak shape should we use for each phase? (from 0 to 1, where 0 is
# pure Gaussian and 1 is pure Lorentzian)
peakshapedict = {'hcp FeNi': 0.5}

# What was the x-ray wavelenth?
wavelength = 0.860255