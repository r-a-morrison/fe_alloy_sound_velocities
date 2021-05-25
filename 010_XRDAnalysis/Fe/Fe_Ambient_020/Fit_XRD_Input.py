# What file are we fitting?
file_name = 'Fe_Ambient_020'

# What NRIXS experiment does this correspond to?
NRIXSexp = 'Fe_Ambient'

# How many background parameters do we need to fit the background?
bg_par_num = 3

# What range should we fit data in?
range_2theta = [23,26]
primaryphase = 'bcc Fe'

# What phases are present, and what is their structure?
phasedict = {'bcc Fe': 'bcc'}

# What starting peak intensities should we use for each phase? (in counts)
peakintdict = {'bcc Fe': 100} # Rough guess

# How many peaks for each phase?
peaknumdict = {'bcc Fe': 1}

# What starting lattice parameters should we use for each phase?
# hcp should have a and c while bcc should have just a
latpardict = {'bcc Fe': {'a':2.863}}
#			{'hcp FeNi': {'a':2.460,'c':3.970}}

# What starting peak width should we use for each phase? (in 2theta degrees)
peakwidthdict = {'bcc Fe': 0.05}

# What starting peak shape should we use for each phase? (from 0 to 1, where 0 is
# pure Gaussian and 1 is pure Lorentzian)
peakshapedict = {'bcc Fe': 0.5}

# What was the x-ray wavelenth?
wavelength = 0.860255

# Are there extra non-sample peaks to fit?
extrapeaknum = 0
extrapeakloc = []
extrapeakint = []
extrapeakwidth = []
extrapeakshape = []