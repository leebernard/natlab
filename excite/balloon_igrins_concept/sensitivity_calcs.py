'''
Sources:
https://www.gemini.edu/observing/resources/near-ir-resources/spectroscopy
https://www.gemini.edu/observing/telescopes-and-sites/sites#Transmission
https://www.iac.es/sites/default/files/documents/2018-06/pwv_sat.pdf
'''

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

from astropy.constants import h, c, k_B


def B_lambda(wavelength, T):

    return 2*h*c**2/wavelength**5 * 1/(np.exp(h*c/(wavelength*k_B*T)) - 1)


def B_nu(freq, T):

    return 2*h*freq**3/c**2 * 1/(np.exp(h*freq/(k_B*T)) - 1)


# calculate mirror contribution
wl = 2159e-9*u.m
wl_delta = 390e-9*u.m
T = 280*u.K
emiss = 0.02
area_gem = np.pi*(8.5*u.m)**2
omega_px = (.5*u.arcsecond)**2
print(f'{B_lambda(wl, T)*1e-9:1.2}')

print(f'{B_lambda(wl, T) * wl * emiss:1.2}')
px_flux = B_lambda(wl, T)/u.steradian * emiss  #omega_px.to(u.steradian)
px_phots = px_flux * wl/(h*c)
keck_mirror_sig = px_phots.to(1/u.s*1/u.arcsecond**2*1/u.um*1/u.m**2)
print(f'{keck_mirror_sig:.2e} photons')

px_phots = B_lambda(wl, T=250*u.K)/u.steradian * emiss * wl/(h*c)
balloon_mirror_sig = px_phots.to(1/u.s*1/u.arcsecond**2*1/u.um*1/u.m**2)

# calculate sky contribution
# pulled from plot
sky_sig = 1e5/78 * 1/(u.s*u.um*u.arcsecond**2*u.m**2)

# calculate host star signal
# wasp 76
kmag = 8.243  # k band magnitude


# https://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/pet/magtojy/
flux_star = 2.16e-13 * u.J/(u.s * u.m**2 * u.um)
# flux_star = B_lambda(wl, T=6115*u.K) * (1.97 * 6.95700e8 /(190 * 3.0857e16))**2
star_sig = flux_star.to(u.J/(u.s*u.um*u.m**2)) * wl/(h*c)



# calculate the signal for a pixel
R_lamba = 40000
wl_bin = wl.to(u.um)/R_lamba
px_bin = 1/2*wl_bin  # assume nyquist sampling

n_spacial_px = 2
px_star_sig = star_sig * px_bin * n_spacial_px
# px_planet_sig = planet_sig * px_bin * n_spacial_px


# calculate typical doppler shift from the target list
# calculate the time for that doppler shift to move by 0.5 delta-lambda

# c = 299792458  # m/s
M_j = 1.89813  # kg/jupiter mass
GM_sol = 1.327124e20  # m^3 s^-2 standard graviational parameter in solar masses

# R = lambda/delta-lambda = 40,000
# delta-V_max = 1/2*c/R
# dV = K_p * 2pi * dphase
# dphase = dt/t_orbit
# delta_t_max = 1/2 * 1/(2pi K_p) * c/R * t_orbit


def delta_t_max(orbital_velocity, t_orbit, R=40000):
    # this is technically modulated by a sin(i), where i is the inclination angle of the exoplanet orbit
    return 0.5 * 1/(2*np.pi * orbital_velocity) * c.value/R * t_orbit


# def B_lambda_nounit(wavelength, T):
#     return 2*h*c**2/wavelength**5 * 1/(np.exp(h*c/(wavelength*k_B*T)) - 1)

# test parameters
# WASP 76 b
period = 1.810  # days
i = 88  # inclination of orbit, +1.3 -1.6 degrees. Close enough to 90
T_eq = 2160  # kelvin
m_p = 0.92  # jupiter masses
m_star = 0.78  # solar masses
rp_rstar = 0.109

t_orb = period*86400  # convert to seconds from days
K_p = (GM_sol * m_star * 2*np.pi/t_orb)**(1/3)  # in m/s
print(f'Orbital velocity: {K_p/1000:.1f} km/s')
max_exp_t = delta_t_max(K_p, t_orb)
print(f'length of time bin: {max_exp_t:.1f} s')


# calculate the planet emission
# assume a black-body with no spectral features or phase-dependency
# i.e., a blank, well-mixed atmosphere
planet_flux_approx = B_lambda(wl, T_eq*u.K) * (1.97 * 6.95700e8 /(190 * 3.0857e16) * rp_rstar)**2 * 3 # fudge factor to make it match star signal
planet_sig = planet_flux_approx.to(u.J/(u.s*u.um*u.m**2)) * wl/(h*c)
px_planet_sig = planet_sig * px_bin * n_spacial_px

extinction = 0.98**15
# telescope area
sn_target = 250

# telescope_area = sn_target**2/(px_star_sig * max_exp_t*u.s)
#
# telescope_dim = np.sqrt(telescope_area/np.pi)*2
telescope_dim = 1.2*u.m
telescope_area = np.pi * (telescope_dim/2)**2
diff_lim = (wl/telescope_dim).to(u.arcsecond, equivalencies=u.dimensionless_angles())

seeing_area = np.pi**2/32400 * 1/60**4
background_change = 8.5**2 * seeing_area / (4*np.pi * (2.2e-6)**2)

background_mag_delta = -2.5*np.log10(background_change)

star_background_contrast = star_sig/(sky_sig + balloon_mirror_sig)
planet_background_contrast = planet_sig/(sky_sig + balloon_mirror_sig)
planet_background_contrast_igrins = planet_sig/((sky_sig + keck_mirror_sig) * background_change )

print('Results:')
print(f'maximum integration time: {max_exp_t:.2e} s, {max_exp_t/60:.2e} mins')
print(f'host star/background contrast: {star_background_contrast:.2f}')
print(f'planet signal/background contrast: {planet_background_contrast:.2f}')
print(f'planet signal/background contrast with igrins: {planet_background_contrast_igrins:.2f}')
print(f'total signal per pixel: {px_star_sig* telescope_area:.2f}')
print(f'planet signal per pixel (esitmate): {px_planet_sig* telescope_area:.2f}')
print(f'source S/N {np.sqrt(px_star_sig * telescope_area * max_exp_t*u.s)}')
print(f'planet emission S/N {px_planet_sig* telescope_area * max_exp_t*u.s/np.sqrt(px_star_sig * telescope_area * max_exp_t*u.s)}')


print(background_mag_delta)



