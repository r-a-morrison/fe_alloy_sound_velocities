# What file are we fitting?
file_name = 'FeNiSi_DAC13_P2b_031'

# What NRIXS experiment does this correspond to?
NRIXSexp = 'FeNiSi_DAC13_P2'

# How many background parameters do we need to fit the background?
bg_par_num = 4

# What range should we fit data in?
range_2theta = [22,28.5]
primaryphase = 'hcp FeNiSi'

# What phases are present, and what is their structure?
phasedict = {'hcp FeNiSi': 'hcp',
             'extrapeak': 'bcc'}

# What starting peak intensities should we use for each phase? (in counts)
peakintdict = {'hcp FeNiSi': 100,
               'extrapeak': 100} # Rough guess

# How many peaks for each phase?
peaknumdict = {'hcp FeNiSi': 3,
               'extrapeak': 1}

# What starting lattice parameters should we use for each phase?
# hcp should have a and c while bcc should have just a
latpardict = {'hcp FeNiSi': {'a':2.409,'c':3.858},
               'extrapeak': {'a':2.98}}
#			{'hcp FeNi': {'a':2.460,'c':3.970}}

# What starting peak width should we use for each phase? (in 2theta degrees)
peakwidthdict = {'hcp FeNiSi': 0.05,
                 'extrapeak': 0.05}

# What starting peak shape should we use for each phase? (from 0 to 1, where 0 is
# pure Gaussian and 1 is pure Lorentzian)
peakshapedict = {'hcp FeNiSi': 0.5,
                 'extrapeak': 0.5}

# What was the x-ray wavelenth?
wavelength = 0.860255