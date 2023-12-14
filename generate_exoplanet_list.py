import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np

from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun, Longitude, Latitude
from astropy.time import Time

verbose = False

cols_needed = ('pl_name, hostname, hip_name, sy_vmag, sy_kmag, pl_bmassj, pl_bmasse, ra, dec, '
               'pl_orbper, pl_orbpererr1, pl_orbpererr2,'
               'pl_trandep, pl_trandeperr1, pl_trandeperr2, pl_trandeplim')

criteria = 'pl_trandep > 0.0 and sy_vmag < 12 and sy_kmag < 11 and pl_orbper < 7'
exotable = NasaExoplanetArchive.query_criteria_async(table='pscomppars', select=cols_needed, where=criteria).to_table()

test_loc = SkyCoord(exotable['ra'][0]*u.deg, exotable['dec'][0]*u.deg, frame='icrs')

# location geo:34.401944,-104.194722, 4032 ft (1220 m)
fort_sumner = EarthLocation(lat=34.401944*u.deg, lon=-104.194722*u.deg, height=1220*u.m)
alt_upper_lim = 57  # *u.deg
alt_lower_lim = 22  # *u.deg


utc_offset = -6*u.hour  # Eastern Daylight Time
start_time = Time('2024-9-6 14:00:00') - utc_offset
end_time = Time('2024-9-7 08:00:00') - utc_offset
flight_time = end_time - start_time  # this gives flight time in days
flight_time = flight_time.to_value(u.hour)  # convert flight time to hours
samples_per_hour = 12
min_time_obsv = 4  # in hours

# generate time samples on a per-hour basis, with the frequency determined by samples_per_hour
sept_hours = np.linspace(start_time, end_time, num=round(flight_time * samples_per_hour))

sun_loc = get_sun(sept_hours)

# seperation should be +- 50 degrees from anti-sun azimuth
# sun_altaz = sun_loc.transform_to(AltAz(obstime=time, location=fort_sumner))
sun_altaz = sun_loc.transform_to(AltAz(obstime=sept_hours, location=fort_sumner))
if verbose:
    print('sun location', sun_loc)
    print(f"sun's Altitude = {sun_altaz.alt}")
    print(f"sun's Azimuth = {sun_altaz.az}")

anti_sun = sun_altaz.az/u.deg - 180
range = 50

min_az = anti_sun - range
max_az = anti_sun + range

# astronomical twilight is defined as the Sun being 12-18° below the horizon. <-18° is complete darkness
# 12° below the horizon (-12°) is the start of observable conditions.
# if alt < -12°, sun position does not matter

exo_skycoord = np.array([SkyCoord(target['ra']*u.deg, target['dec']*u.deg, frame='icrs').
                         transform_to(AltAz(obstime=sept_hours, location=fort_sumner)) for target in exotable])

# # transpose the array, for efficiency
# exo_skycoord = exo_skycoord.transpose()

# sun_altaz.alt > -12
# get_az = np.vectorize(lambda x: x.az.value)
# exo_az = get_az(exo_skycoord)
# extract the altitude and azimuth cooridnates, because arrays of objects sucks
altaz_cube = np.array([[[altaz.alt/u.deg, altaz.az/u.deg] for altaz in altaz_list] for altaz_list in exo_skycoord])

alt = altaz_cube[:,:, 0]
az = altaz_cube[:,:,1]
is_anti_sun = (az < max_az[0]) & (az > min_az[0])

is_night = ~ (sun_altaz.alt > -12*u.deg)

is_in_elevation = (alt_upper_lim > alt) & (alt > alt_lower_lim)


# valid_exotable = [exoplanet for (exoplanet, valid) in zip(exotable, is_anti_sun*is_in_elevation) if valid]
# the above is cute, but wrong. the correct way is below
# combine all the validation checks
valid_altaz = is_in_elevation * (is_anti_sun + is_night[None, :])
# find out which exoplanet candidates survived
if verbose:
    print('time instances valid', valid_altaz.sum(axis=1))
    print('amounut of time observable', valid_altaz.sum(axis=1)/samples_per_hour)

# check if the sum of valid time samples meets the minimum observing time
is_valid = valid_altaz.sum(axis=1) > min_time_obsv*samples_per_hour

'''
I need a way of calculating how long each valid target is observable.
More accurately, I need to calculate when each target enters and exits observability.
If I sample the observing time frequently enough, that provides enter and exit times to a reasonable approximation.
'''

# set all non-valid altaz values to zero
valid_altaz_cube = altaz_cube*valid_altaz[..., None]
# remove all targets that don't meet criteria
valid_altaz_cube = valid_altaz_cube[is_valid, ...]

test = (valid_altaz_cube != 0 ).argmax(axis=1)
test2 = valid_altaz_cube[np.arange(valid_altaz_cube.shape[0]), test[:, 0], :]

# min_sep = (180 - 50)*u.deg  # minimum target separation from the sun in degrees
# has_min_seperation = np.array([sun_loc.separation(SkyCoord(target['ra'], target['dec'], frame='icrs')) > min_sep for target in exotable])

valid_exotable = exotable[is_valid]

verbose = True
if verbose:
    print("Number of planets in database:", len(exotable))
    print('Number of planets observable between', start_time, 'and', end_time, ':', len(valid_exotable))

if verbose:
    print(valid_exotable['pl_trandep'])

fig, ax = plt.subplots(tight_layout=True)

ax.hist(exotable['pl_trandep'].data, bins=100)


'''Compare my list to Peter's list'''
list_path = '/home/lee/natlab/excite_targets/flight_20240906_140000_UT-6phase3A.csv'

from astropy.io import ascii
from astropy.table import join

# open Peter's table
peter_table = ascii.read(list_path, header_start=11, delimiter=',')

print(peter_table['System'])
print('compared to my list:')
print(valid_exotable['hostname'])





