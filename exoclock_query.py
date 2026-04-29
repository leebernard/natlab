"""
This script is for accessing the ExoClock database, and pulling exoplanet data
for EXCITE's target list.



ExoClock:
https://www.exoclock.space/database/planets

The top level dictionary is a list of exoplanets, with each key being the exoplanet name.
"""

import urllib
import json
import re

# from astropy.io import ascii
from astropy.table import Table

verbose = True

exoclock_planets = json.loads(urllib.request.urlopen('https://www.exoclock.space/database/planets_json').read())

file_path = '/home/lee/nat_lab/excite_mission/excite_targets/'

targets_filename = 'flight_2026122118.00.00_UT-6phase2A.csv'

tepcat_table = Table.read(file_path+targets_filename, format='ascii.csv', header_start=10, data_start=11)

print('exoclock planet names')
print(exoclock_planets.keys())

print('EXCITE target table keys')
print(tepcat_table.colnames)

exoclock_system_names = ' '.join(exoclock_planets.keys())
excite_target_systems = tepcat_table['System']

# check how many EXCITE targets are in the ExoClock list
# regular expression pattern
is_leading_0 = re.compile(r'(?<![0-9])0')
for target_system in excite_target_systems:
    if verbose:
        print(target_system)
    # print(type(target_system))
    # strip underscore characters
    target_system = target_system.replace('_', '')
    target_system = re.sub(is_leading_0, '', target_system)
    if verbose:
        print('After removing characters:', target_system)
        print('Is in exoclock?', target_system in exoclock_system_names)

# Now pull the locations where ExoClock has EXCITE targets
excite_keys = []
for target_system in excite_target_systems:
    target_system = target_system.replace('_', '')
    target_system = re.sub(is_leading_0, '', target_system)
    if verbose:
        print('EXCITE target name after removing characters:', target_system)
    exoclock_sys = [system for system in exoclock_planets.keys() if target_system in system]
    print(exoclock_sys)
    if exoclock_sys:
        excite_keys.append(exoclock_sys[0])

excite_target_ephemeris = [exoclock_planets[excite_key] for excite_key in excite_keys]



