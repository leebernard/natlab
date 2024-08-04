import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter
from astropy.constants import h, c, k_B

debug = False

# constants
rp_rstar = 0.109  # ratio of the planet radius to star radius
r_instrument = 2500  # R value of the instrument

def B_lambda(wavelength, T):
    # produces units of J / (s m3)
    return 2*h*c**2/wavelength**5 * 1/(np.exp(h*c/(wavelength*k_B*T)) - 1)


"""open files"""
gemini_trans_file = 'excite/balloon_igrins_concept/mktrans_zm_16_15.dat'
mk_trans_data = np.loadtxt(gemini_trans_file)

gemini_em_file = 'excite/balloon_igrins_concept/mk_skybg_zm_16_15_ph.dat'
mk_emi_data = np.loadtxt(gemini_em_file)  # ph/sec/arcsec^2/nm/m^2. Assumes 273 K

oh_em_file = 'excite/balloon_igrins_concept/irlinespec1.txt'
oh_em_data = np.loadtxt(oh_em_file)  # units are (um, W cm^-2 str^-1 um^-1)

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

"""end open files"""

"""begin spectrum generation"""
# define bandpass
short_limit = 1.45
long_limit = 2.6

bandpass = np.array([short_limit, long_limit])

idmk0 = np.where(mk_trans_data[:, 0] == bandpass[0])[0][0]
idmk1 = np.where(mk_trans_data[:, 0] == bandpass[1])[0][0]
wavelengths = mk_trans_data[idmk0:idmk1+1, 0]


# calculate the resolution of the sample
sample_resolution = np.mean(wavelengths/np.diff(wavelengths, prepend=mk_trans_data[idmk0-1, 0]))  # average R value of sample
# bin_width = np.mean(np.diff(wavelengths))

filter_sigma = sample_resolution/r_instrument

# filter the data to the resolution specified by r_instrument
mk_em = gaussian_filter(mk_emi_data[idmk0:idmk1+1, 1], sigma=filter_sigma) * 1000 # convert from 1/nm to 1/um
mk_trans = gaussian_filter(mk_trans_data[idmk0:idmk1+1, 1], sigma=filter_sigma)  # assumes mk_trans uses the same wavelengths as mk_em

idmm0 = np.absolute(mcmurdo_trans_data[:, 0] - short_limit*1000).argmin()
idmm1 = np.absolute(mcmurdo_trans_data[:, 0] - long_limit*1000).argmin()
wavelengths_mcmurdo = mcmurdo_trans_data[idmm0:idmm1+1, 0] / 1000  # convert to um from nm
mcmurdo_trans = mcmurdo_trans_data[idmm0:idmm1+1, 1] # slice the data array to the bandpass

peter_wl = mcmurdo_em_data[idmm0:idmm1+1, 0]  # nm
# convert from uW cm^-2 nm^-1 sr^-1
# to photons/s arcsec^-2 um^-1 m^-2            J/uW   m/nm    nm       phots/(J m)         cm^2/m^2 nm/um   sr/arcsec^2
peter_em = mcmurdo_em_data[idmm0:idmm1+1, 1] * 1e-6 * 1e-9 * peter_wl/(h.value*c.value) * 100**2 * 1000 * 1/4.25e10  # photons/s arcsec^-2 um^-1 m^-2

peter_sample_r = np.mean(wavelengths_mcmurdo/np.diff(wavelengths_mcmurdo, prepend=mcmurdo_trans_data[idmm0-1, 0] / 1000))  # average R value of sample
# test = wavelengths_em == wavelengths_mcmurdo

idoh0 = np.absolute(oh_em_data[:, 0] - short_limit).argmin()  # 1.5 um, in angstroms
idoh1 = np.absolute(oh_em_data[:, 0] - long_limit).argmin()  # 2.6 um, in angstroms
oh_wl = oh_em_data[idoh0:idoh1+1, 0]  # um
# Watt cm^(-2) str^(-1) micron^(-1)   (cm/m)^2  sr/arcsec^2  m/um  um    phots/(J m)
oh_em = oh_em_data[idoh0:idoh1+1, 1] * 100**2 * 1/4.25e10 * 1e-6 * oh_wl/(h.value*c.value)  # Phtons/s arcsec^-2 um^-1 m^-2

# upsample
peter_upsample = np.interp(oh_wl, peter_wl/1000, peter_em)

mcmurdo_em = peter_upsample + oh_em


planet_flux = B_lambda(wavelengths*1e-6 * u.m, T=2100*u.K).to(u.J/(u.s*u.um*u.m**2)) * (1.97 * 6.95700e8 /(190 * 3.0857e16) * rp_rstar)**2 * 3 # fudge factor to make it match star signal
mcmurdo_planet_flux = B_lambda(wavelengths_mcmurdo*1e-6 * u.m, T=2100*u.K).to(u.J/(u.s*u.um*u.m**2)) * (1.97 * 6.95700e8 /(190 * 3.0857e16) * rp_rstar)**2 * 3 # fudge factor to make it match star signal
star_flux = B_lambda(wavelengths*1e-6 * u.m, T=6100*u.K).to(u.J/(u.s*u.um*u.m**2)) * (1.97 * 6.95700e8 /(190 * 3.0857e16))**2

if debug:
    fig, ax = plt.subplots()
    ax.plot(oh_em[:, 0], oh_em[:, 1])
"""End spectrum generation"""

"""Begin flux calculations"""
mirror_d_mk = 8.1  # m
seeing_mk = np.radians(0.34/3600)  # 0.34 arcsecond slit. Slit width of igrins at Gemini South. Units: rad/bin
mirror_d_mcmurdo = 1.2  # m
diffraction_mcmurdo = wavelengths[0]/mirror_d_mcmurdo * 1e-6  # convert um to m. Units rad/bin
px_sampling = 2

# wl_bin = 1.45/r_instrument  # resolution element width. Units are delta-um/bin
px_radians_mk = 1/px_sampling*seeing_mk  # units: rad/px
px_radians_mcmurdo = 1/px_sampling*diffraction_mcmurdo

sky_area_mk = np.pi/4 * px_radians_mk**2  # sr of a cone with angle equal to the pixel angle at Mauna Kea
sky_area_mk = sky_area_mk*4.25e10  # convert to arcsec^2
sky_area_mcmurdo = np.pi/4 * px_radians_mcmurdo**2  # sr of a cone with angle equal to the balloon pixel angle
sky_area_mcmurdo = sky_area_mcmurdo*4.25e10  # convert to arcsec^2

if debug:
    fig, ax = plt.subplots()
    ax.plot(wavelengths, sky_area_mcmurdo)
    ax.axhline(sky_area_mk, color='red')





# convert to photons
# wl_width = 2.e-5*u.um
planet_spectrum = planet_flux * wavelengths*1e-6*u.m /(h*c) * 1/px_sampling**2
mcmurdo_planet_spectrum = mcmurdo_planet_flux * wavelengths_mcmurdo * 1e-6*u.m / (h*c) * 1/px_sampling**2
star_blackbody = star_flux * wavelengths*1e-6*u.m /(h*c)  * 1/px_sampling**2

R = 40000
fig, (axmk, axmcmurdo) = plt.subplots(2, sharex=True, tight_layout=True, figsize=(12,12))

gem_area = np.pi * (mirror_d_mk/2)**2
integration = 70
axmk.plot(wavelengths, mk_em * sky_area_mk * wavelengths/(2*R), label='Sky Background, Mauna Kea')
axmk.plot(wavelengths, planet_spectrum * mk_trans * wavelengths/(2*R), label='exoplanet blackbody, from Mauna Kea', color='C1')
axmk.plot(wavelengths, star_blackbody * wavelengths/(2*R), label='Stellar blackbody at top of atmosphere', color='C3', linewidth=2.5, linestyle='dotted')

# axmk.set_xlim(1.45, 2.5)
# axmk.set_ylim(1e-4, 15)
axmk.set_ylabel('Flux (photons/sec/pixel/m^2)')
# axmk.set_xlabel(f'Wavelength (um), R~{r_instrument}')
axmk.set_yscale('log')
axmk.legend()



axmcmurdo.plot(oh_wl, mcmurdo_em * sky_area_mcmurdo * oh_wl/(2*R), label='Sky Background, 40 km above McMurdo', color='mediumblue')  #, color='tab:purple', linewidth=2.5)
axmcmurdo.plot(wavelengths_mcmurdo, mcmurdo_planet_spectrum * mcmurdo_trans * wavelengths_mcmurdo/(2*R), label='exoplanet blackbody, 40 km above McMurdo', color='C1')
# ax.plot(wavelengths, planet_spectrum, label='exoplanet blackbody, top of atmosphere', color='C1', linewidth=2.5, linestyle='dotted')
axmcmurdo.plot(wavelengths, star_blackbody * wavelengths/(2*R), label='Stellar blackbody at top of atmosphere', color='C3', linewidth=2.5, linestyle='dotted')

axmcmurdo.set_xlim(1.45, 2.5)
# axmcmurdo.set_ylim(1e-4, 15)
axmcmurdo.set_ylabel('Flux (photons/sec/pixel/m^2)')
axmcmurdo.set_xlabel(f'Wavelength (um), smoothed to R~{r_instrument}')
axmcmurdo.set_yscale('log')
axmcmurdo.legend()


# transmission comparision zoom
fig2, (ax1, ax2) = plt.subplots(2, tight_layout=True, figsize=(12, 12))
linewidth=2.0

ax1.plot(wavelengths, mk_trans, label='Sky Transmission, Mauna Kea', linewidth=linewidth)
ax1.plot(wavelengths_mcmurdo, mcmurdo_trans, label='Sky Transmission, 40 km above McMurdo', linewidth=linewidth)
# ax1.plot(peter_mk_trans_data[:, 0] / 1000, peter_mk_trans_data[:, 1], label='Mauna Kea, Peter\'s version', linewidth=2.0)

ax1.set_xlabel('Wavelength (um)')
ax1.set_ylabel('Transmissivity')
ax1.legend()

ax2.plot(wavelengths, mk_trans, label='Sky Transmission, Mauna Kea', linewidth=linewidth)
ax2.plot(wavelengths_mcmurdo, mcmurdo_trans, label='Sky Transmission, 40 km above McMurdo', linewidth=linewidth)

ax2.set_xlim(1.45, 2.5)
y_ticks = [.98, .99, 1.00]
ax2.set_yticks(y_ticks)
# ax2.yaxis.set_major_formatter('{y_ticks:.2f}')
ax2.set_ylim(0.978, 1.002)

ax2.set_xlabel('Wavelength (um)')
ax2.set_ylabel('Zoomed to 98% region')
# ax2.legend()




# fig, ax = plt.subplots(tight_layout=True)
# ax.plot(oh_em[:, 0]/10000, oh_em[:, 1]*1e2)



