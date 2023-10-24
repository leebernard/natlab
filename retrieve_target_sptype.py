import numpy as np
# import re
# from astropy.table import setdiff

from astroquery.simbad import Simbad

file_path = '/home/lee/Downloads/FTSnight_20230919/FTSnight_20230919/'

file_name = 'visstars_FTS_202309020_000000_UT-6.csv'

midnight_targets = np.loadtxt(file_path+file_name, skiprows=12, delimiter=',')

hip_ids = midnight_targets[:, 0]

hip_ids = np.char.mod('%d', hip_ids)

hip_ids = np.char.add('HIP ', hip_ids)

from filter_star_targets import retrieve_targetdata

target_table = retrieve_targetdata(hip_ids,
                                   add_fields=['typed_id', 'ra(d;A;ICRS;J2017.5;2000)', 'dec(d;D;ICRS;J2017.5;2000)',
                                               'sp','flux(V)', 'flux(I)', 'otype', 'otypes'],
                                   remove_fields=['coordinates'])


target_table.pprint_all()
