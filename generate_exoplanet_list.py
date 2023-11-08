import matplotlib.pyplot as plt

from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

cols_needed = ('pl_name, hostname, hip_name, sy_vmag, sy_zmag, pl_bmassj, pl_bmasse, ra, dec, '
               'pl_orbper, pl_orbpererr1, pl_orbpererr2,'
               'pl_trandep, pl_trandeperr1, pl_trandeperr2, pl_trandeplim,'
               'pl_orbper, pl_orbperstr')

criteria = 'pl_trandep > 0 and sy_vmag < 11 and dec < 77.5 and dec > -33.5 and pl_orbper < 1.5'
exolist = NasaExoplanetArchive.query_criteria(table='pscomppars', select=cols_needed, where=criteria)

print(len(exolist))

print(exolist['pl_trandep'])

fig, ax = plt.subplots(tight_layout=True)

ax.hist(exolist['pl_trandep'].data)




