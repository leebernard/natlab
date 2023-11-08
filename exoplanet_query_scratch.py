'''
This is for testing querying the NASA Exoplanet Archive using astroquery.
'''
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time

# location geo:34.401944,-104.194722, 4032 ft (1220 m)
fort_sumner = EarthLocation(lat=34.401944*u.deg, lon=-104.194722*u.deg, height=1220*u.m)

utcoffset = -6*u.hour  # Eastern Daylight Time
time = Time('2024-9-15 18:00:00') - utcoffset

# altitude limits
alt_lims = (57*u.deg, 22*u.deg)
alt_lowerlim = alt_lims[-1]

# lower limit due south at midnight
sky_point = SkyCoord(alt=alt_lowerlim, az=180*u.deg, obstime=time, frame = 'altaz', location=fort_sumner)
print(sky_point.transform_to('icrs'))

# upper limit due north at midnight
sky_point = SkyCoord(alt=alt_lowerlim, az=0.0*u.deg, obstime=time, frame = 'altaz', location=fort_sumner)
print(sky_point.transform_to('icrs'))

NasaExoplanetArchive.TAP_TABLES

test = NasaExoplanetArchive.query_object("WASP-18 b", table='pscomppars', select='*')

test.colnames

cols_needed = ('pl_name, hostname, hip_name, sy_vmag, sy_zmag, pl_bmassj, pl_bmasse, ra, dec, '
               'pl_orbper, pl_orbpererr1, pl_orbpererr2,'
               'pl_trandep, pl_trandeperr1, pl_trandeperr2, pl_trandeplim,'
               'pl_orbper, pl_orbperstr')

test2 = NasaExoplanetArchive.query_object("WASP-18 b", table='pscomppars', select=cols_needed)

test2.pprint_all()

criteria = 'pl_trandep > 0 and sy_vmag < 11 and dec < 77.5 and dec > -33.5 and pl_orbper < 1.5'
test3 = NasaExoplanetArchive.query_criteria(table='pscomppars', select=cols_needed, where=criteria)

print(len(test3))

print(test3['pl_trandep'])

import matplotlib.pyplot as plt

fig, ax = plt.subplots(tight_layout=True)

ax.hist(test3['pl_trandep'].data)


