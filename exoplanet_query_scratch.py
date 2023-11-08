'''
This is for testing querying the NASA Exoplanet Archive using astroquery.
'''
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

NasaExoplanetArchive.TAP_TABLES

test = NasaExoplanetArchive.query_object("WASP-18 b", table='pscomppars', select='*')

test.colnames

cols_needed = ('pl_name, hostname, hip_name, sy_vmag, sy_zmag,'
               'pl_trandep, pl_trandeperr1, pl_trandeperr2, pl_trandeplim,'
               'pl_orbper, pl_orbperstr')

test2 = NasaExoplanetArchive.query_object("WASP-18 b", table='pscomppars', select=cols_needed)

test2.pprint_all()

criteria = 'pl_trandep > 0.5 and sy_vmag < 11'
test3 = NasaExoplanetArchive.query_criteria(table='pscomppars', select=cols_needed, where=criteria)

print(len(test3))

test3['pl_trandep']


