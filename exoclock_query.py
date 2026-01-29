"""
This script is for accessing the ExoClock database, and comparing it to EXCITE's target list.

ExoClock:
https://www.exoclock.space/database/planets

The top level dictionary is a list of exoplanets, with each key being the exoplanet name.
"""

import urllib
import json

# from astropy.io import ascii
from astropy.table import Table

exoclock_planets = json.loads(urllib.request.urlopen('https://www.exoclock.space/database/planets_json').read())

file_path = '/home/lee/nat_lab/excite_mission/excite_targets/'

targets_filename = 'flight_2026122118.00.00_UT-6phase2A.csv'

tepcat_table = Table.read(file_path+targets_filename, format='ascii.csv', header_start=10, data_start=11)


