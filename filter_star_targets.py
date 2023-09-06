import numpy as np
import re
# from astropy.table import setdiff

from astroquery.simbad import Simbad
# import csv


def retrieve_targetdata(id_list, add_fields, remove_fields, debug=False):

    query = Simbad()
    query.remove_votable_fields(*remove_fields)
    query.add_votable_fields(*add_fields)
    result_table = query.query_objects(id_list)


    if debug:
        query.get_votable_fields()
        result_table.pprint_all()

    return result_table


def is_valid(target_table, otype_regex, debug=False):

    pattern = re.compile(otype_regex)
    if debug:
        print('Pattern', pattern.pattern)
    is_variable = []
    for row in target_table:
        if debug:
            print(row)

        check = bool(re.search(pattern, row['OTYPES']))
        if debug:
            print('Checking if star matches flag', check)
        is_variable.append(check)

    return np.asarray(is_variable)


def validate_targets(target_list, otype_flags,
                     add_fields=['typed_id', 'ra(d;A;ICRS;J2017.5;2000)', 'dec(d;D;ICRS;J2017.5;2000)', 'otype',
                                 'otypes'],
                     remove_fields=['coordinates'],
                     debug=False
                     ):

    data_table = retrieve_targetdata(target_list, add_fields=add_fields, remove_fields=remove_fields, debug=debug)

    return is_valid(data_table, otype_flags, debug=debug)


def is_target_variable(target_list, otype_flags=r'V\*|Ir\*|Er\*|Ro\*|Pu\*',
                     add_fields=['typed_id', 'ra(d;A;ICRS;J2017.5;2000)', 'dec(d;D;ICRS;J2017.5;2000)', 'otype',
                                 'otypes'],
                     remove_fields=['coordinates'],
                     debug=False):
    return validate_targets(target_list, otype_flags, add_fields=add_fields, remove_fields=remove_fields, debug=debug)

if __name__ == "__main__":
    debug = False

    path = '/home/lee/natlab/excite_targets/'
    # file = 'simbad_engineering_2023-09-04.tsv'

    id_file = 'engineering_2023-09-04.txt'

    with open(path+id_file) as file:
        target_list = [line.rstrip() for line in file]

    # retrieve a table of target data from Simbad
    target_table = retrieve_targetdata(target_list,
                                       add_fields=['typed_id', 'ra(d;A;ICRS;J2017.5;2000)', 'dec(d;D;ICRS;J2017.5;2000)', 'otype','otypes'],
                                       remove_fields=['coordinates'])
    if debug:
        target_table.pprint_all()

    # generate a binary filter according to variable flags
    var_flags = r'V\*|Ir\*|Er\*|Ro\*|Pu\*'
    bin_flags = r'El\*|EB\*'

    is_variable_star = is_valid(target_table, var_flags)
    is_binary_variable = is_valid(target_table, bin_flags)

    # is_target_variable() combines the above two functions:
    test3 = is_target_variable(target_list)



# # HIP ID, Vmag, RA (deg), DEC (deg), AZ (deg), EL (deg)
# hip_id, Vmag, Ra, DEC, AZ, EL = np.loadtxt(path+file, delimiter=',').transpose()
#
# hip_id = hip_id.astype('int')
#
# for id, mag, alt in zip(hip_id, Vmag, EL):
#     print(id, mag, alt)


