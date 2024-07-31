import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter
from astropy.constants import h, c, k_B

debug = False


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

mcmurdo_trans_file = 'excite/balloon_igrins_concept/transmmsday_40p0km.dat'
mcmurdo_trans_data = np.flip(np.loadtxt(mcmurdo_trans_file), axis=0)  # data file is in long to short wavelength order

mcmurdo_em_file = 'excite/balloon_igrins_concept/radmmsday_40p0km.dat'
mcmurdo_em_data = np.flip(np.loadtxt(mcmurdo_em_file), axis=0)

peter_mk_trans_file = 'excite/balloon_igrins_concept/transmaunanite_4p2km.dat'
peter_mk_trans_data = np.flip(np.loadtxt(peter_mk_trans_file), axis=0)

if debug:
    fig, ax = plt.subplots()
    ax.plot(mcmurdo_trans_data[:, 0]/1000, mcmurdo_trans_data[:, 1], label='Sky transmission, 40 km above McMurdo')
    ax.plot(mk_trans_data[:, 0], mk_trans_data[:, 1], label='Sky transmission, Mauna Kea', linewidth=2.0)
    ax.plot(peter_mk_trans_data[:, 0]/1000, peter_mk_trans_data[:, 1], label='Mauna Kea, Peter\'s version', linewidth=2.0)
    ax.set_xlabel('Wavelength (um)')

    ax.set_xlim(1.45, 2.6)
    # ax.set_ylim(0.5)
    ax.legend()



short_limit = 1.45
long_limit = 2.6

bandpass = np.array([short_limit, long_limit])


idmk0 = np.where(mk_trans_data[:, 0] == bandpass[0])[0][0]
idmk1 = np.where(mk_trans_data[:, 0] == bandpass[1])[0][0]
wavelengths = mk_trans_data[idmk0:idmk1+1, 0]

sample_resolution = np.mean(wavelengths/np.diff(wavelengths, prepend=mk_trans_data[idmk0-1, 0]))  # average R value of sample

bin_width = np.mean(np.diff(wavelengths))

desired_r = 360
filter_sigma = sample_resolution/desired_r

mk_em = gaussian_filter(mk_emi_data[idmk0:idmk1+1, 1], sigma=filter_sigma) * 1000 # convert from 1/nm to 1/um
mk_trans = gaussian_filter(mk_trans_data[idmk0:idmk1+1, 1], sigma=filter_sigma)  # assumes mk_trans uses the same wavelengths as mk_em

idmm0 = np.absolute(mcmurdo_trans_data[:, 0] - short_limit*1000).argmin()
idmm1 = np.absolute(mcmurdo_trans_data[:, 0] - long_limit*1000).argmin()
wavelengths_mcmurdo = mcmurdo_trans_data[idmm0:idmm1+1, 0] / 1000  # convert to um from nm
mcmurdo_trans = mcmurdo_trans_data[idmm0:idmm1+1, 1] # slice the data array to the bandpass
# convert from uW cm^-2 nm^-1 sr^-1
# to photons/s arcsec^-2 um^-1 m^-2          J/uW   m/um  phots/(J m)        cm^2/m^2 um/m   sr/arcsec^2
mcmurdo_em = mcmurdo_em_data[idmm0:idmm1+1, 1] * 1e-6 * 1e-6/(h.value*c.value) * 100**2 * 1000 * 1/4.25e10

peter_sample_r = np.mean(wavelengths_mcmurdo/np.diff(wavelengths_mcmurdo, prepend=mcmurdo_trans_data[idmm0-1, 0] / 1000))  # average R value of sample
# test = wavelengths_em == wavelengths_mcmurdo

# idx0 = np.absolute(oh_em_data[:, 0] - short_limit*10000).argmin()  # 1.5 um, in angstroms
# idx1 = np.absolute(oh_em_data[:, 0] - long_limit*10000).argmin()  # 2.6 um, in angstroms
# oh_em = oh_em_data[id0:id1+1]

planet_flux = B_lambda(wavelengths*1e-6 * u.m, T=2100*u.K).to(u.J/(u.s*u.um*u.m**2)) * (1.97 * 6.95700e8 /(190 * 3.0857e16) * rp_rstar)**2 * 3 # fudge factor to make it match star signal
mcmurdo_planet_flux = B_lambda(wavelengths_mcmurdo*1e-6 * u.m, T=2100*u.K).to(u.J/(u.s*u.um*u.m**2)) * (1.97 * 6.95700e8 /(190 * 3.0857e16) * rp_rstar)**2 * 3 # fudge factor to make it match star signal
star_flux = B_lambda(wavelengths*1e-6 * u.m, T=6100*u.K).to(u.J/(u.s*u.um*u.m**2)) * (1.97 * 6.95700e8 /(190 * 3.0857e16))**2


# convert to photons
# wl_width = 2.e-5*u.um
planet_spectrum = planet_flux * wavelengths*1e-6*u.m /(h*c)
mcmurdo_planet_spectrum = mcmurdo_planet_flux * wavelengths_mcmurdo * 1e-6*u.m / (h*c)
star_blackbody = star_flux * wavelengths*1e-6*u.m /(h*c)


fig, ax = plt.subplots(tight_layout=True, figsize=(10, 6))
ax.plot(wavelengths, mk_em, label='Sky Emission, Mauna Kea')
ax.plot(wavelengths, planet_spectrum*mk_trans, label='exoplanet blackbody, from Mauna Kea', color='C2')

ax.plot(wavelengths_mcmurdo, mcmurdo_em, label='Sky Emission, 40 km above McMurdo', color='tab:purple', linewidth=2.5)
ax.plot(wavelengths_mcmurdo, mcmurdo_planet_spectrum*mcmurdo_trans, label='exoplanet blackbody, 40 km above McMurdo', color='C7', linewidth=2.5)

ax.plot(wavelengths, planet_spectrum, label='exoplanet blackbody, top of atmosphere', color='C1', linewidth=2.5, linestyle='dotted')
ax.plot(wavelengths, star_blackbody, label='Stellar blackbody (for reference)', color='C3', linewidth=2.5, linestyle='dotted')

ax.set_xlim(1.45, 2.5)
ax.set_ylim(10)
ax.set_ylabel('Flux (photons/sec/arcsec^2/um/m^2)')
ax.set_xlabel(f'Wavelength (um), R~{desired_r}')
ax.set_yscale('log')
ax.legend()


fig2, (ax1, ax2) = plt.subplots(2, tight_layout=True)

ax1.plot(wavelengths, mk_trans, label='Sky Transmission, Mauna Kea')
ax1.plot(wavelengths_mcmurdo, mcmurdo_trans, label='Sky Transmission, 40 km above McMurdo')
# ax1.plot(peter_mk_trans_data[:, 0] / 1000, peter_mk_trans_data[:, 1], label='Mauna Kea, Peter\'s version', linewidth=2.0)

ax1.set_xlim(1.45, 2.5)
ax1.set_xlabel('Wavelength (um)')
ax1.set_ylabel('Transmissivity')
ax1.legend()

ax2.plot(wavelengths, mk_trans, label='Sky Transmission, Mauna Kea')
ax2.plot(wavelengths_mcmurdo, mcmurdo_trans, label='Sky Transmission, 40 km above McMurdo')

ax2.set_xlim(1.45, 2.5)
y_ticks = [.98, .99, 1.00]
ax2.set_yticks(y_ticks)
# ax2.yaxis.set_major_formatter('{y_ticks:.2f}')
ax2.set_ylim(0.978, 1.002)

ax2.set_xlabel('Wavelength (um)')
ax2.set_ylabel('Zoomed to 98% region')
ax2.legend()




# fig, ax = plt.subplots(tight_layout=True)
# ax.plot(oh_em[:, 0]/10000, oh_em[:, 1]*1e2)



