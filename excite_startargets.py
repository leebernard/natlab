import numpy as np
import re
from astropy.table import setdiff

from astroquery.simbad import Simbad
# import csv

debug = False


def retrieve_targetdata(filename,
                        add_fields=['typed_id', 'ra(d;A;ICRS;J2017.5;2000)', 'dec(d;D;ICRS;J2017.5;2000)', 'otype', 'otypes'],
                        remove_fields=['coordinates'],
                        debug=False):

    with open(filename) as file:
        id_list = [line.rstrip() for line in file]

    query = Simbad()
    query.remove_votable_fields(*remove_fields)
    query.add_votable_fields(*add_fields)
    result_table = query.query_objects(id_list)

    if debug:
        if debug:
            query.get_votable_fields()
            result_table.pprint_all()

    return result_table


def validate_targets(target_table, otype_flags, debug=False, inclusive_filter=True):

    pattern = re.compile('|'.join(otype_flags))
    is_variable = []
    for row in target_table:
        if debug:
            print(row)

        check = bool(re.search(pattern, row['OTYPES']))
        if debug:
            print('Checking if star matches flag', check)
        is_variable.append(check)

    is_valid = np.asarray(is_variable)
    if inclusive_filter is False:
        # return the targets that Do Not have the flags
        return result_table[~is_valid]
    else:
        # return targets that have the flags
        return result_table[is_valid]


path = '/home/lee/natlab/excite_targets/'

file = 'simbad_engineering_2023-09-04.tsv'

id_file = 'engineering_2023-09-04.txt'

test = retrieve_targetdata(path+id_file)


with open(path+id_file) as file:
    excite_list = [line.rstrip() for line in file]

excite_query = Simbad()
excite_query.add_votable_fields('typed_id')
excite_query.remove_votable_fields('coordinates')
excite_query.add_votable_fields('ra(d;A;ICRS;J2017.5;2000)', 'dec(d;D;ICRS;J2017.5;2000)')
excite_query.add_votable_fields('otype', 'otypes')
result_table = excite_query.query_objects(excite_list)

if debug:
    excite_query.get_votable_fields()
    result_table.pprint_all()

# data = np.loadtxt(path+file, dtype='str', delimiter='\t', skiprows=7, max_rows=21)
#
# data_list = data.tolist()

# re_pattern = r'(V\*)|(Ir\*)|(Er\*)|(Ro\*)|(Pu\*)'
var_flags = r'V\*|Ir\*|Er\*|Ro\*|Pu\*'
bin_flags = r'El\*|EB\*'
var_pattern = re.compile(var_flags)
bin_pattern = re.compile(bin_flags)

debug = True
is_variable = []
is_binary = []
for row in result_table:
    if debug:
        print(row)

    check = bool(re.search(var_pattern, row['OTYPES']))
    if debug:
        print('Checking if star is variable', check)
    is_variable.append(check)

    check = bool(re.search(bin_pattern, row['OTYPES']))
    if debug:
        print('Checking if star is binary variable', check)
    is_binary.append(check)

# clean up some memeory
# del data_list
is_variable = np.asarray(is_variable)
is_binary = np.asarray(is_binary)

validated_targets = result_table[~is_variable]

binary_variable = result_table[is_binary]




# # HIP ID, Vmag, RA (deg), DEC (deg), AZ (deg), EL (deg)
# hip_id, Vmag, Ra, DEC, AZ, EL = np.loadtxt(path+file, delimiter=',').transpose()
#
# hip_id = hip_id.astype('int')
#
# for id, mag, alt in zip(hip_id, Vmag, EL):
#     print(id, mag, alt)


