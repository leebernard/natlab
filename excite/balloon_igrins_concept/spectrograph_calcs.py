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

D_gem = 8.1
D_super = 0.5
D_giga = 1.35

wl_igrins = np.linspace(1.47, 2.5, num=40000)  # bandpass wavelengths, um

phi_gem = np.radians(0.34/3600)  # arcseconds to radians. Typicak k-band Gemini seeing is 0.6-0.4 arcseconds
F_igrins = 8.8
n_grating = 3.435  # at 130 Kelvin, longward of 1 um
delta_igrins = np.radians(71.56)  # R3 grating
# delta_r1 = np.arctan(1)  # R1 grating example
grating_density = 36.5  # number of lines/mm
hband_orders = np.arange(98, 122+1)  #, step=2)  # number of orders to cover k band
kband_orders = np.arange(72, 92+1)  # number of orders to cover h band
lam_b_hband = np.array((1823.59, 1788.03, 1753.89, 1721.08, 1689.52, 1659.15, 1629.15, 1601.71, 1574.53, 1548.30, 1522.97, 1498.50, 1474.85))  # nm
lam_b_kband = np.array((2467.88, 2402.10, 2339.80, 2224.62, 2171.27, 2120.49, 2072.09, 2025.92, 1981.81, 1939.65))  # nm

f_camera = 126.6  # mm
px_sampling = 3.66  # pixels per wavelength element
d1_igrins = 25  # mm

F_2 = f_camera/d1_igrins

'''Calculate the space the spectrum takes up'''
# assume theta is 1 degree
theta = np.radians(1.4)

alpha_igrins = delta_igrins + theta
beta_blaze = delta_igrins - theta
r_igrins = np.cos(alpha_igrins)/np.cos(beta_blaze)
# r should be ~1
print(f'r={r_igrins:.2f}')
fold1_clearance = 2 * F_igrins * d1_igrins * theta  # small angle approximation
print(f'f1 fold mirror clearance: {fold1_clearance:.2f} mm')

W_igrins = d1_igrins/np.sin(alpha_igrins)
sigma_igrins = 1/grating_density

R_gem = 2*d1_igrins*1e-3 * n_grating * np.tan(delta_igrins) / (phi_gem * D_gem * r_igrins)
print(f'Resolution R={R_gem:g}')

lam_blaze_igrins = 2*sigma_igrins/hband_orders * n_grating*np.sin(delta_igrins) * 1e6  # convert from mm to nm
# note: the above calculation is slightly off from IGRINS published numbers, with the magnitude of error depending on
# the wavelength. This is most likely due to wavelength dependency of n

free_spectral_range = lam_blaze_igrins/hband_orders

'''
Calculate the specs for the proposed hi-res instrument. 
Start with the direct specs of IGRINS, but with full band coverage (no gap between H and K band)
Calculate how much detector space I need.
Then start tweaking parameters to see if I can fit it onto one detector.
A good place to start is noting the spectrum in IGRINS is over-sampled, at 3.66 pixels per wavelength element
'''

hires_orders = np.arange(72, 122+1)
R_hires = 40000
delta_hires = delta_igrins
sigma_hires = sigma_igrins
alpha_hires = delta_hires + theta
beta_b_hires = delta_hires - theta
px_pitch = 18.  # um

wl_crit = 1.5e-6  # wavelength driving slit criteria
m_crit = np.round(n_grating * sigma_hires*1e-3 / wl_crit * (np.sin(beta_b_hires) + np.sin(alpha_hires)))
print(f'criteria order: {m_crit}')

phi_super = wl_crit/D_super  # set slit wide to diffract limit at 1.5 um
print(f'hires slit wide: {np.degrees(phi_super)*3600: .2f} arcseconds')

d1_super = R_hires * sigma_hires/m_crit * np.cos(alpha_hires)
print(f'Hires collimated beam diameter, half meter telescope: {d1_super:.2f} mm')

f2_super = d1_super / (D_super) * 2*px_pitch*1e-6 / phi_super
print(f'f camera, {D_super} m telescope: {f2_super} mm')

phi_giga = wl_crit/D_giga
d1_giga = R_hires * sigma_hires/m_crit * np.cos(alpha_hires)
print(f'Hires collimated beam diameter, 1.35 meter telescope: {d1_giga:.2f} mm')

f2_giga = d1_giga / (D_giga) * 2*px_pitch*1e-6 / phi_giga
print(f'f camera, {D_giga} m telescope: {f2_giga} mm')






