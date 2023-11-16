import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np

from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time

cols_needed = ('pl_name, hostname, hip_name, sy_vmag, sy_zmag, pl_bmassj, pl_bmasse, ra, dec, '
               'pl_orbper, pl_orbpererr1, pl_orbpererr2, pl_orbperstr,'
               'pl_trandep, pl_trandeperr1, pl_trandeperr2, pl_trandeplim')

criteria = 'pl_trandep > 0.0 and sy_vmag < 11 and pl_orbper < 7'
exolist = NasaExoplanetArchive.query_criteria(table='pscomppars', select=cols_needed, where=criteria)

test_loc = SkyCoord(exolist['ra'][0], exolist['dec'][0], frame='icrs')

# location geo:34.401944,-104.194722, 4032 ft (1220 m)
fort_sumner = EarthLocation(lat=34.401944*u.deg, lon=-104.194722*u.deg, height=1220*u.m)
alt_upper_lim = 57*u.deg
alt_lower_lim = 22*u.deg

utc_offset = -6*u.hour  # Eastern Daylight Time
time = Time('2024-9-1 09:00:00') - utc_offset

sept_hours = time + np.arange(30*24)*u.hour

sun_loc = get_sun(time)
print('sun location', sun_loc)

# seperation should be +- 50 degrees from anti-sun azimuth
sun_altaz = sun_loc.transform_to(AltAz(obstime=time,location=fort_sumner))
print(f"sun's Altitude = {sun_altaz.alt:.5}")
print(f"sun's Azimuth = {sun_altaz.az:.5}")

anti_sun = sun_altaz.az - 180*u.deg
range = 50*u.deg

min_az = anti_sun - range
max_az = anti_sun + range

exolist_azalt = [SkyCoord(target['ra'], target['dec'], frame='icrs').transform_to(AltAz(obstime=time,location=fort_sumner)) for target in exolist]
is_anti_sun = np.array([max_az > target.az > min_az for target in exolist_azalt])
is_in_elevation = np.array([alt_upper_lim > target.alt > alt_lower_lim for target in exolist_azalt])

# valid_exolist = [exoplanet for (exoplanet, valid) in zip(exolist, is_anti_sun*is_in_elevation) if valid]
# the above is cute, but wrong. the correct way is below
valid_exolist = exolist[is_in_elevation * is_anti_sun]

# min_sep = (180 - 50)*u.deg  # minimum target separation from the sun in degrees
# has_min_seperation = np.array([sun_loc.separation(SkyCoord(target['ra'], target['dec'], frame='icrs')) > min_sep for target in exolist])

print("Number of planets in database:", len(exolist))
print('Number of planets observable at', time, ':', len(valid_exolist))

print(valid_exolist['pl_trandep'])

fig, ax = plt.subplots(tight_layout=True)

ax.hist(exolist['pl_trandep'].data, bins=100)



