import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

from astropy.constants import h, c, k_B

rp_rstar = 0.109

def B_lambda(wavelength, T):

    return 2*h*c**2/wavelength**5 * 1/(np.exp(h*c/(wavelength*k_B*T)) - 1)

# open files

gemini_trans_file = 'excite/balloon_igrins_concept/mktrans_zm_16_15.dat'
mk_trans_data = np.loadtxt(gemini_trans_file)

gemini_em_file = 'excite/balloon_igrins_concept/mk_skybg_zm_16_15_ph.dat'
mk_emi_data = np.loadtxt(gemini_em_file)

oh_em_file = 'excite/balloon_igrins_concept/rousselot2000.dat'
oh_em_data = np.loadtxt(oh_em_file)

bandpass = np.array([1.5, 2.6])

id0 = np.where(mk_trans_data[:, 0] == bandpass[0])[0][0]
id1 = np.where(mk_trans_data[:, 0] == bandpass[1])[0][0]
wavelengths = mk_trans_data[id0:id1+1, 0]

mk_em = mk_emi_data[id0:id1+1, 1] * 1000 # convert from 1/nm to 1/um
mk_trans = mk_trans_data[id0:id1+1, 1]

idx0 = np.absolute(oh_em_data[:, 0] - 15000).argmin()
idx1 = np.absolute(oh_em_data[:, 0] - 26000).argmin()
oh_em = oh_em_data[id0:id1+1]

planet_flux = B_lambda(wavelengths*1e-6 * u.m, T=2100*u.K).to(u.J/(u.s*u.um*u.m**2)) * (1.97 * 6.95700e8 /(190 * 3.0857e16) * rp_rstar)**2 * 3 # fudge factor to make it match star signal
star_flux = B_lambda(wavelengths*1e-6 * u.m, T=6100*u.K).to(u.J/(u.s*u.um*u.m**2)) * (1.97 * 6.95700e8 /(190 * 3.0857e16))**2

# convert to photons
wl_width = 2.e-5*u.um
planet_spectrum = planet_flux* wavelengths*1e-6*u.m /(h*c)
star_blackbody = star_flux * wavelengths*1e-6*u.m /(h*c)


fig, ax = plt.subplots(tight_layout=True)
ax.plot(wavelengths, mk_em, label='Sky Emission, Mauna Kea')
# ax.plot(wavelengths, mk_trans, label='Sky absorption ')
# ax.plot(wavelengths, star_blackbody*mk_trans, label='stellar blackbody, Mauna Kea')
ax.plot(wavelengths, planet_spectrum, label='exoplanet blackbody, top of atmosphere')
ax.plot(wavelengths, planet_spectrum*mk_trans, label='exoplanet blackbody, Mauna Kea')
ax.plot(wavelengths, star_blackbody, label='Stellar blackbody, top of atomsphere')
ax.set_ylabel('ph/sec/arcsec^2/um/m^2')
ax.set_xlabel('Wavelength (um), R~100,000')
ax.set_yscale('log')
ax.legend()

fig, ax = plt.subplots(tight_layout=True)
ax.plot(wavelengths, mk_trans)

fig, ax = plt.subplots(tight_layout=True)
ax.plot(oh_em[:, 0]/10000, oh_em[:, 1]*1e2)



