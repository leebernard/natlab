"""
This file is for calculating the optical parameters of the hi-res instrument
First, calculate the parameters for the igrins instrument.
Then translate those parameters to two different versions of hi-rest: one with a superbit (0.5 m telescope) platform, and the
other with a gigabit (1.35 m telescope) platform.
Reuse appropiate parameters from igrins, such as the grating and the resolution. Tweak as little as possible, to start.

Igrins has R~45,000, but only ~ 20,000 resolution elements. This is due to the gap between H and K, caused by the
atmospheric absorption peak. They require a 2K x 2K detector to image both spectrums

Number of illuminated grooves determines the number of resolution elements per order. IGRINs uses a lithographic
Si immersed grating. Si immersion increases the effective resolution by n (3.4). Uses a single detector.
Immersion grating has a 60 mK stability requirement, to prevent instablility in the index of refraction of Si.
"""

import numpy as np

d_gem = 8.1
d_excite = 0.5
d_giga = 1.35

phi_gem = np.radians(0.34/3600)  # arcseconds to radians. Typicak k-band Gemini seeing is 0.6-0.4 arcseconds
F_gem = 10
n_grating = 3.43  # at 130 Kelvin, longward of 1 um
# Need:
# grating spacing sigma
# blaze angle
# alternate, what modes they are using
alpha_gem = np.radians(71.56)  # they call this R3
grating_density = 36.5  # number of lines/mm
kband_orders = np.arange(98, 122+1)  # number of orders to cover k band
hband_orders = np.arange(72, 92+1)  # number of orders to cover h band

f_camera = 126.6  # mm
px_sampling = 3.66  # pixels per wavelength element
d1_gem = 25  # mm

