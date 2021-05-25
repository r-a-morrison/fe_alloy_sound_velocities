# What file are we fitting?
file_name = 'FeNi_DAC11_P7b_026'

# What NRIXS experiment does this correspond to?
NRIXSexp = 'FeNi_DAC11_P7'

# How many background parameters do we need to fit the background?
bg_par_num = 3

# What range should we fit data in?
range_2theta = [22,29]
primaryphase = 'hcp FeNi'

# What phases are present, and what is their structure?
phasedict = {'hcp FeNi': 'hcp',
             'extrapeak': 'bcc'}

# What starting peak intensities should we use for each phase? (in counts)
peakintdict = {'hcp FeNi': 100,
               'extrapeak': 100} # Rough guess

# How many peaks for each phase?
peaknumdict = {'hcp FeNi': 3,
               'extrapeak': 1}

# What starting lattice parameters should we use for each phase?
# hcp should have a and c while bcc should have just a
latpardict = {'hcp FeNi': {'a':2.36,'c':3.75},
               'extrapeak': {'a':2.965}}
#			{'hcp FeNi': {'a':2.460,'c':3.970}}

# What starting peak width should we use for each phase? (in 2theta degrees)
peakwidthdict = {'hcp FeNi': 0.05,
                 'extrapeak': 0.05}

# What starting peak shape should we use for each phase? (from 0 to 1, where 0 is
# pure Gaussian and 1 is pure Lorentzian)
peakshapedict = {'hcp FeNi': 0.5,
                 'extrapeak': 0.5}

# What was the x-ray wavelenth?
wavelength = 0.860255