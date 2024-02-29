import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import re

from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun, Longitude, Latitude
from astropy.time import Time
from numpy import sin, cos, arccos


def great_angle(ra1, dec1, ra2, dec2):
    # assume radians
    return arccos(sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra2 - ra1))


verbose = False

'''query database for potential exoplanet targets'''
cols_needed = ('pl_name, hostname, hip_name, sy_vmag, sy_kmag, pl_bmassj, pl_bmasse, ra, dec, '
               'disc_facility, disc_telescope, disc_instrument, disc_refname,'
               'pl_orbper, pl_orbpererr1, pl_orbpererr2,'
               'pl_trandep, pl_trandeperr1, pl_trandeperr2, pl_trandeplim')

criteria = 'pl_trandep > 0.0 and sy_vmag < 12 and sy_kmag < 11 and pl_orbper < 7'
exotable = NasaExoplanetArchive.query_criteria_async(table='pscomppars', select=cols_needed, where=criteria).to_table()

# test_loc = SkyCoord(exotable['ra'][0]*u.deg, exotable['dec'][0]*u.deg, frame='icrs')


''' set constants'''
# location geo:34.401944,-104.194722, 4032 ft (1220 m)
fort_sumner = EarthLocation(lat=34.401944*u.deg, lon=-104.194722*u.deg, height=1220*u.m)
alt_upper_lim = 57  # *u.deg
alt_lower_lim = 23  # *u.deg


utc_offset = -6*u.hour  # Eastern Daylight Time
start_time = Time('2024-9-6 14:00:00') - utc_offset
end_time = Time('2024-9-7 08:00:00') - utc_offset
flight_time = end_time - start_time  # this gives flight time in days
flight_time = flight_time.to_value(u.hour)  # convert flight time to hours
samples_per_hour = 12
min_day_obsv = 2  # in hours
min_night_obsv = 4  # in hours


'''Filter the table according to observability'''
# generate time samples on a per-hour basis, with the frequency determined by samples_per_hour
sept_hours = np.linspace(start_time, end_time, num=round(flight_time * samples_per_hour))

sun_loc = get_sun(sept_hours)

sun_altaz = sun_loc.transform_to(AltAz(obstime=sept_hours, location=fort_sumner))
if verbose:
    print('sun location', sun_loc)
    print(f"sun's Altitude = {sun_altaz.alt}")
    print(f"sun's Azimuth = {sun_altaz.az}")

# seperation should be +- 50 degrees from anti-sun azimuth
az_range = 50  # plus/minus from the antisun postion

anti_sun = Longitude(sun_altaz.az - 180*u.deg)  # using astropy Longitude to take advantage of wrapping
min_az = Longitude(anti_sun - az_range*u.deg)/u.deg

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
# extract the altitude and azimuth cooridnates, because arrays of objects suck
altaz_cube = np.array([[[altaz.alt/u.deg, altaz.az/u.deg] for altaz in altaz_list] for altaz_list in exo_skycoord])

alt = altaz_cube[:,:, 0]
az = altaz_cube[:,:,1]

# validation checks: all should be bool arrays

is_anti_sun = (az - min_az)%360 < az_range*2  # check if the az is within 100 degrees of min_antisun

is_night = ~ (sun_altaz.alt > -12*u.deg)  # check if it's astronomical night

is_in_elevation = (alt < alt_upper_lim) & (alt > alt_lower_lim)  # check if the telescope can point at the target


# valid_exotable = [exoplanet for (exoplanet, valid) in zip(exotable, is_anti_sun*is_in_elevation) if valid]
# the above is cute, but wrong. the correct way is below
# combine all the validation checks
is_valid_daytime_altaz = is_in_elevation * is_anti_sun * ~is_night[None, :]

is_valid_nighttime_altaz = is_in_elevation * is_night[None, :]

# check if the observation window is long enough
is_day_obstime = is_valid_daytime_altaz.sum(axis=1) > min_day_obsv*samples_per_hour

# check continuity
day_altaz_bool = is_valid_daytime_altaz[is_day_obstime]
dayobs_start_index = day_altaz_bool.argmax(axis=1)
dayobs_end_index = (day_altaz_bool.shape[1] - np.fliplr(day_altaz_bool).argmax(axis=1) - 1)
is_continuous = []
for n in range(day_altaz_bool.shape[0]):
    start_i = dayobs_start_index[n]
    end_i = dayobs_end_index[n]
    is_continuous.append(day_altaz_bool[n, start_i:end_i].any() == False)


is_night_valid = is_valid_nighttime_altaz.sum(axis=1) > min_night_obsv*samples_per_hour
# find out which exoplanet candidates survived
verbose = False
if verbose:
    print('Number observable during the day', is_day_valid.sum())
    print('Number observable during the night', is_night_valid.sum())

is_valid = is_day_valid+is_night_valid


# extract the start and stop times of observability
# test = is_valid_daytime_altaz[is_day_valid]
#
# # extract the windows at which the valid targets are observable
# test = (is_valid_daytime_altaz[is_day_valid] != 0 ).argmax(axis=1)


daytime_start_times = np.ma.masked_array(sept_hours[is_valid_daytime_altaz.argmax(axis=1)], mask=~is_day_valid, fill_value=np.ma.masked)
daytime_end_times = np.ma.masked_array(np.flip(sept_hours)[(np.fliplr(is_valid_daytime_altaz)).argmax(axis=1)], mask=~is_day_valid, fill_value=np.ma.masked)
exotable['daytime_window_start_time'] = daytime_start_times
exotable['daytime_window_end_time'] = daytime_end_times

# should I combine daytime and nighttime hours?

# test = daytime_end_times - daytime_start_times
# test = test.to_value(u.hour)
# print(test)

nighttime_start_times = np.ma.masked_array(sept_hours[(is_valid_nighttime_altaz).argmax(axis=1)], mask=~is_night_valid, fill_value=np.ma.masked)
nighttime_end_times = np.ma.masked_array(np.flip(sept_hours)[(np.fliplr(is_valid_nighttime_altaz)).argmax(axis=1)], mask=~is_night_valid, fill_value=np.ma.masked)
exotable['nighttime_window_start_time'] = nighttime_start_times
exotable['nighttime_window_end_time'] = nighttime_end_times
# test = nighttime_end_times - nighttime_start_times
# print(test.to_value(u.hour))


# generate filtered list of targets
valid_exotable = exotable[is_valid]

is_kepler = valid_exotable['disc_telescope'] == '0.95 m Kepler Telescope'
is_canon = valid_exotable['disc_telescope'] == 'Canon 200mm f/1.8L'
min_depth = 0.1  # in percentage
is_validdepth = valid_exotable['pl_trandep'] >= min_depth

if verbose:
    valid_exotable['disc_facility'].pprint(max_lines=-1)
    valid_exotable['disc_telescope'].pprint(max_lines=-1)

# extract array of ra and dec
valid_ra = valid_exotable['ra']
valid_dec = valid_exotable['dec']

non_kepler_ra = valid_exotable['ra'][~is_kepler]
non_kepler_dec = valid_exotable['dec'][~is_kepler]

verbose = False
if verbose:
    print("Number of planets in database:", len(exotable))
    print('Number of planets observable between', start_time, 'and', end_time, ':', len(valid_exotable))

if verbose:
    print(valid_exotable['pl_trandep'])


fig, ax = plt.subplots(tight_layout=True)

ax.hist(exotable['pl_trandep'].data, bins=100)


'''Compare my list to Peter's list'''
# note:
# Peter is using TEPCat (https://www.astro.keele.ac.uk/jkt/tepcat/)
# I am using the NASA Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu/)
list_path = '/home/lee/natlab/excite_targets/flight_20240906_140000_UT-6phase3A.csv'

from astropy.io import ascii
from astropy.table import join

# open Peter's table
peter_table = ascii.read(list_path, header_start=11, delimiter=',')
valid_3day_table = ascii.read('target_list_3dayorbits.csv')

if verbose:
    print(peter_table['System'])
    print('compared to my list:')
    print(valid_exotable['hostname'])

# careful comparision between the lists using RA and Dec

peter_ra = peter_table['RA [deg]']
peter_dec = peter_table['dec [deg]']

ra_3day = valid_3day_table['ra']
dec_3day = valid_3day_table['dec']

closest_match = []
for ra, dec in zip(peter_ra, peter_dec):
    closest = great_angle(np.radians(ra), np.radians(dec), np.radians(valid_ra), np.radians(valid_dec))
    closest_match.append(np.degrees(closest.min()))
closest_match = np.array(closest_match)

# assume anything less than 1 arcmin is the same object
is_overlap = closest_match < 1/60
overlap_targets = closest_match[is_overlap]
overlap_table = peter_table[is_overlap]

fig, ax = plt.subplots(tight_layout=True)

ax.scatter(valid_ra, valid_dec, label='My targets')
# ax.scatter(non_kepler_ra, non_kepler_dec, label='My targets (excluding Kepler)')
# ax.scatter(valid_ra[is_kepler], valid_dec[is_kepler], label='Kepler discoveries', marker='s', s=60, color='tab:red', facecolors='none')
# ax.scatter(valid_ra[is_canon], valid_dec[is_canon], label='Canon 200mm f/1.8L', marker='s', s=60, color='tab:red', facecolors='none')
ax.scatter(valid_ra[is_validdepth], valid_dec[is_validdepth], label=f'Transit depth > {min_depth}', marker='s', s=60, color='tab:red', facecolors='none')

ax.scatter(peter_ra, peter_dec, label='Peter\'s targets', marker='+', color='tab:orange')

ax.legend()
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.set_title('Sky projection')
ax.set_ylim(-90, 90)
ax.set_xlim(0, 360)

'''write out the results'''

valid_exotable.write('target_list.csv', format='csv')

with np.printoptions(threshold=np.inf, linewidth=200):
    print(day_altaz_bool)

