import numpy as np
import astropy.units as u

from astropy.constants import h, c, k_B


def B_lambda(wavelength, T):

    return 2*h*c**2/wavelength**5 * 1/(np.exp(h*c/(wavelength*k_B*T)) - 1)


def B_nu(freq, T):

    return 2*h*freq**3/c**2 * 1/(np.exp(h*freq/(k_B*T)) - 1)


# calculate mirror contribution
wl = 2199e-9*u.m
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
sky_sig = 1e5/78 * 1/(u.s*u.um*u.arcsecond**2*u.m**2)

# calculate host star signal
# wasp 76
T_star = 6250 * u.K
r_star = 1.73 * u.R_sun
distance = 195 * u.pc
R = distance.to(u.R_sun)
star_flux = B_lambda(wl, T_star)/u.steradian * wl/(h*c)
star_sig = (r_star/R)**2 * star_flux.to(1/(u.s*u.um*u.arcsecond**2*u.m**2))


# calculate typical doppler shift from the target list
# calculate the time for that doppler shift to move by 0.5 delta-lambda

c = 299792458  # m/s
M_j = 1.89813  # kg/jupiter mass
GM_sol = 1.327124e20  # m^3 s^-2 standard graviational parameter in solar masses

# R = lambda/delta-lambda = 40,000
# delta-V_max = 1/2*c/R
# dV = K_p * 2pi * dphase
# dphase = dt/t_orbit
# delta_t_max = 1/2 * 1/(2pi K_p) * c/R * t_orbit


def delta_t_max(orbital_velocity, t_orbit, R=40000):
    # this is technically modulated by a sin(i), where i is the inclination angle of the exoplanet orbit
    return 0.5 * 1/(2*np.pi * orbital_velocity) * c/R * t_orbit


def B_lambda(wavelength, T):
    return 2*h*c**2/wavelength**5 * 1/(np.exp(h*c/(wavelength*k_B*T)) - 1)

# test parameters
# WASP 76 b
period = 1.810  # days
i = 88  # +1.3 -1.6, degrees. Close enough to 90
T_eq = 2160  # kelvin
m_p = 0.92  # jupiter masses
m_star = 0.78  # solar masses

t_orb = period*86400  # convert to seconds from days
K_p = (GM_sol * m_star * 2*np.pi/t_orb)**(1/3)
print(f'Orbital velocity: {K_p/1000:.1f} km/s')
max_exp_t = delta_t_max(K_p, t_orb)
print(f'length of time bin: {max_exp_t:.1f} s')

print(max_exp_t)





