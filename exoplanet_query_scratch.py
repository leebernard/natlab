'''
This is for testing querying the NASA Exoplanet Archive using astroquery.
'''
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

test = NasaExoplanetArchive.query_object("WASP-18 b", table='pscomppars', select='*')

test.colnames

cols_needed = ('pl_name, hostname, hip_name, sy_vmag, sy_zmag,'
               'pl_trandep, pl_trandeperr1, pl_trandeperr2, pl_trandeplim,'
               'pl_orbper, pl_orbperstr')

test2 = NasaExoplanetArchive.query_object("WASP-18 b", table='pscomppars', select=cols_needed)

test2.pprint_all()

test3 = NasaExoplanetArchive.query_criteria_async(table='pscomppars', select=cols_needed)
