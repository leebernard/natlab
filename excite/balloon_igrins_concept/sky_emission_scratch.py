import numpy as np
import matplotlib
matplotlib.use("qt5agg")
import matplotlib.pyplot as plt
# import astropy.units as u
# from argparse import ArgumentParser
import lowtran
# from lowtran.plot import irradiance

c1 = {
    "model": 1,  # model selection
    "h1": 4.2,  # height
    "angle": 45,  # zenith angle
    "wlshort": 800,  # nm
    "wllong": 3500,  # nm
    "wlstep": 10,
}

#models: 1 tropical 2 midlatitude summer 3 midlatitude winter
# 4 subarctic summer 5 subarctic winter 6 1976 US standard
# mid latitude: 23.5-66.5

test_trans = lowtran.transmittance(c1)

irradiance(irr, c1, True)

plt.show()

seeing = 0.35




