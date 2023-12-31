import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time

debug = True

path = '/home/lee/natlab/excite_targets/'
# file = 'simbad_engineering_2023-09-04.tsv'

id_file = 'engineering_2023-09-04.txt'

with open(path + id_file) as file:
    target_list = [line.rstrip() for line in file]

# hip_id = target_list[0]
# hip_id = 'HIP 86201'
hip_id = 'HIP 94376'
target_coord = SkyCoord.from_name(hip_id)

# geo:34.401944,-104.194722, 4032 ft (1220 m)
fort_sumner = EarthLocation(lat=34.401944*u.deg, lon=-104.194722*u.deg, height=1220*u.m)
utcoffset = -6*u.hour  # Eastern Daylight Time
time = Time('2023-9-20 14:00:00') - utcoffset

target_altaz = target_coord.transform_to(AltAz(obstime=time,location=fort_sumner))
print(f"Target's Altitude = {target_altaz.alt:.5}")
print(f"Target's Azimuth = {target_altaz.az:.5}")




