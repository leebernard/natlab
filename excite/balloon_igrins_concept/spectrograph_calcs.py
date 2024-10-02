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

wl_igrins = np.linspace(1.47, 2.5, num=40000)  # bandpass wavelengths, um

phi_gem = np.radians(0.34/3600)  # arcseconds to radians. Typicak k-band Gemini seeing is 0.6-0.4 arcseconds
F_igrins = 10
n_grating = 3.44  # at 130 Kelvin, longward of 1 um
delta_igrins = np.radians(71.56)  # R3 grating
# delta_r1 = np.arctan(1)  # R1 grating example
grating_density = 36.5  # number of lines/mm
hband_orders = np.arange(98, 122+1, step=2)  # number of orders to cover k band
kband_orders = np.arange(72, 92+1)  # number of orders to cover h band
lam_b_hband = np.array((1823.59, 1788.03, 1753.89, 1721.08, 1689.52, 1659.15, 1629.15, 1601.71, 1574.53, 1548.30, 1522.97, 1498.50, 1474.85))  # nm
lam_b_kband = np.array((2467.88, 2402.10, 2339.80, 2224.62, 2171.27, 2120.49, 2072.09, 2025.92, 1981.81, 1939.65))  # nm

f_camera = 126.6  # mm
px_sampling = 3.66  # pixels per wavelength element
d1_igrins = 25  # mm

F_2 = f_camera/d1_igrins

'''Calculate the space the spectrum takes up'''
# assume theta is 1 degree
theta = np.radians(1)
alpha_igrins = delta_igrins + theta
beta_blaze = delta_igrins - theta
r = np.cos(alpha_igrins)/np.cos(beta_blaze)
# r should be ~1
print(f'r={r:.2f}')

W_igrins = d1_igrins/np.sin(alpha_igrins)
sigma_igrins = 1/grating_density

lam_blaze = 2*sigma_igrins/hband_orders * n_grating*np.sin(delta_igrins) * 1e6  # convert from mm to nm
theta_igrins = np.arccos(lam_b_hband*1e-9 * hband_orders / (2*sigma_igrins*1e-3 * n_grating * np.sin(delta_igrins)))
# R_igrins = wl_igrins/d_gem *

