"""
This script is for accessing the ExoClock database, and comparing it to EXCITE's target list.

ExoClock:
https://www.exoclock.space/database/planets

The top level dictionary is a list of exoplanets, with each key being the exoplanet name.
"""

import urllib
import json

exoclock_planets = json.loads(urllib.request.urlopen('https://www.exoclock.space/database/planets_json').read())



