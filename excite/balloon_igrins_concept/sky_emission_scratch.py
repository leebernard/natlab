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

trans = lowtran.transmittance(c1)
wl = np.array(trans.transmission['wavelength_nm']) * 1.0e-3  # convert to microns

# keck transmission
t_k = np.array(trans.transmission[0, :, 0])
sc_k = np.array(lowtran.scatter(c1)['pathscatter'][0, :, 0])
ra_k = np.array(lowtran.radiance(c1)['radiance'][0, :, 0])

# excite transmission
c1['model'] = 4
c1['h1'] = 38.
trans_mc = lowtran.transmittance(c1)
t_e = np.array(trans_mc.transmission[0, :, 0])
sc_e = np.array(lowtran.scatter(c1)['pathscatter'][0, :, 0])
ra_e = np.array(lowtran.radiance(c1)['radiance'][0, :, 0])

# irradiance(irr, c1, True)
fig, (ax_trans, ax_em) = plt.subplots(2, tight_layout=True, figsize=(9,6))

ax_trans.plot(wl, t_k, label='Mauna Kea')
ax_trans.plot(wl, t_e, label=f'McMurdo, {c1["h1"]} km')
ax_trans.set_xlabel('wavelength (um)')
ax_trans.set_ylabel('Atmosphere transmission')
ax_trans.set_xlim(1.5, 2.5)
ax_trans.legend()

# ax_em.plot(wl, sc_k, label='Mauna Kea scattered light')
ax_em.plot(wl, ra_k, label='MK atmospheric emission')
ax_em.plot(wl, sc_e, label='Antartic scattered light')
ax_em.plot(wl, ra_e, label = 'Antartic Atmospheric emission')
ax_em.set_xlabel('Wavelength (um)')
ax_em.set_ylabel('Radiance (W /ster /cm^2 /um)')
ax_em.set_xlim(1.5, 2.5)
ax_em.legend()

plt.show()

seeing = 0.35




