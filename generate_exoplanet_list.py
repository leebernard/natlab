import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np

from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time

# location geo:34.401944,-104.194722, 4032 ft (1220 m)
fort_sumner = EarthLocation(lat=34.401944*u.deg, lon=-104.194722*u.deg, height=1220*u.m)

utcoffset = -6*u.hour  # Eastern Daylight Time
time = Time('2024-9-15 18:00:00') - utcoffset

sun_loc = get_sun(time)

cols_needed = ('pl_name, hostname, hip_name, sy_vmag, sy_zmag, pl_bmassj, pl_bmasse, ra, dec, '
               'pl_orbper, pl_orbpererr1, pl_orbpererr2,'
               'pl_trandep, pl_trandeperr1, pl_trandeperr2, pl_trandeplim,'
               'pl_orbper, pl_orbperstr')

criteria = 'pl_trandep > 0 and sy_vmag < 11 and dec < 77.5 and dec > -33.5 and pl_orbper < 1.5 and pl_trandep > .1'
exolist = NasaExoplanetArchive.query_criteria(table='pscomppars', select=cols_needed, where=criteria)

test_loc = SkyCoord(exolist['ra'][0], exolist['dec'][0], frame='icrs')

seperation = sun_loc.separation(test_loc)

print('Seperation from sun:', sun_loc.separation(test_loc).deg)

min_sep = 60*u.deg  # minimum target separation from the sun in degrees
dist_to_sun = np.array([sun_loc.separation(SkyCoord(target['ra'], target['dec'], frame='icrs')) > min_sep for target in exolist])

print(len(exolist))

print(exolist['pl_trandep'])

fig, ax = plt.subplots(tight_layout=True)

ax.hist(exolist['pl_trandep'].data, bins=100)



