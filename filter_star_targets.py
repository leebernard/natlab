import numpy as np
import re
# from astropy.table import setdiff

from astroquery.simbad import Simbad
# import csv


def retrieve_targetdata(id_list, add_fields, remove_fields, debug=False, **keywords):
    '''
    Queries Simbad for object data, using a list of identifiers. The data is returned as an astropy table, with table
    columns corresponding to the fields queried from Simbad. The default fields are 'main_id' (the id Simbad uses) and
    'coordinates' (a shorthand for all fields involving Ra and Dec, including Ra/Dec in sexagesimal format, error
    information, and bibcodes).

    Parameters
    ----------
    id_list : list, str
        a list of star identifiers
    add_fields : list, str
        Simbad field identifiers to be added to the query
    remove_fields : list, str
        Simbad field identifiers to be removed from the query
    debug : bool
    keywords : keywords to pass through to the query_objects() method

    Returns
    -------
    An astropy table containing data retrieved from Simbad.
    '''
    query = Simbad()
    query.remove_votable_fields(*remove_fields)
    query.add_votable_fields(*add_fields)
    result_table = query.query_objects(id_list, **keywords)


    if debug:
        print('current votablle fields')
        query.get_votable_fields()
        print('retrieved data table')
        result_table.pprint_all()

    return result_table


def is_valid(target_table, flag_regex, field='OTYPES', debug=False):
    '''
    Checks if the regex pattern matches any entries in the given field.
    Uses re.search(pattern) to search the given field for matches.

    Parameters
    ----------
    target_table : `~astropy.table.Table`
        Table of object data retrieved from Simbad
    flag_regex : str
        String to be compiled into a regex pattern
    field : str
        identifier of the Simbad field to be searched
    debug : bool

    Returns
    -------
    Bool array. True means a match was found.
    '''
    pattern = re.compile(flag_regex)
    if debug:
        print('Pattern', pattern.pattern)
    is_variable = []
    for row in target_table:
        if debug:
            print(row)

        check = bool(re.search(pattern, row[field]))
        if debug:
            print('Checking if star matches flag', check)
        is_variable.append(check)

    return np.asarray(is_variable)


def is_target_variable(target_list, otype_flags=r'V\*|Ir\*|Er\*|Ro\*|Pu\*',
                     add_fields=['typed_id', 'ra(d;A;ICRS;J2017.5;2000)', 'dec(d;D;ICRS;J2017.5;2000)', 'otype',
                                 'otypes'],
                     remove_fields=['coordinates'],
                     debug=False):
    '''
    Wrapper function for retrieve_targetdata() and is_valid(). Pass this function a list of identifiers, and it will
    validate whether the targets are variable stars.

    It uses the list of target identifiers to query the Simbad database, and check if any of the targets are variable
    stars. It then returns a bool array.

    Parameters
    ----------
    target_list : list, str
        List of target identifiers to check for if they are variable.
    otype_flags : str, optional
        A string to be compiled as a regex pattern. This regex pattern searches the 'otypes' field for matches.
    add_fields : list, str, optional
        Fields to be added to the Simbad data dump
    remove_fields : list, str, optional
        Fields to be removed from the Simbad data dump
    debug : bool, optional
        turns on debugging mode

    Returns
    -------
    A bool numpy array, where True means the corresponding target is a variable star according to Simbad.
    Returned array corresponds positionally to the input array.

    Notes
    -----
    See http://simbad.cds.unistra.fr/guide/sim-fscript.htx for flag and field identifiers.
    '''

    data_table = retrieve_targetdata(target_list, add_fields=add_fields, remove_fields=remove_fields, debug=debug)

    return is_valid(data_table, otype_flags, debug=debug)


def retrieve_spectral_types(target_list, debug=False, **keywords):
    data_table = retrieve_targetdata(target_list, add_fields=['sp'], remove_fields=['coordiantes'], debug=debug, **keywords)

    return np.asarray(data_table['SP_TYPE'])


def retrieve_sptype_and_variable(target_list, otype_flags=r'V\*|Ir\*|Er\*|Ro\*|Pu\*',
                                add_fields=['typed_id',
                                            'ra(d;A;ICRS;J2017.5;2000)',
                                            'dec(d;D;ICRS;J2017.5;2000)',
                                            'sp',
                                            'otype',
                                            'otypes'],
                                remove_fields=['coordinates'],
                                debug=False):
    data_table = retrieve_targetdata(target_list, add_fields=add_fields, remove_fields=remove_fields, debug=debug)

    return is_valid(data_table, otype_flags, debug=debug), np.asarray(data_table['SP_TYPE'])



if __name__ == "__main__":
    debug = True

    path = '/home/lee/natlab/excite_targets/'
    # file = 'simbad_engineering_2023-09-04.tsv'

    id_file = 'engineering_2023-09-04.txt'

    with open(path+id_file) as file:
        target_list = [line.rstrip() for line in file]

    # retrieve a table of target data from Simbad
    target_table = retrieve_targetdata(target_list,
                                       add_fields=['typed_id', 'ra(d;A;ICRS;J2017.5;2000)', 'dec(d;D;ICRS;J2017.5;2000)', 'flux(V)', 'sp','otype','otypes'],
                                       remove_fields=['coordinates'])
    if debug:
        target_table.pprint_all()

    # generate a binary filter according to variable flags
    spect_class_filter = r'V'  # any 'V' luminosity star is main sequence
    var_flags = r'V\*|Ir\*|Er\*|Ro\*|Pu\*'
    bin_flags = r'El\*|EB\*'

    is_variable_star = is_valid(target_table, var_flags)
    is_binary_variable = is_valid(target_table, bin_flags)

    # is_target_variable() is a wrapper function combines the above two functions:
    shortcut_to_is_variable_star = is_target_variable(target_list)
    print('Ran sucessfully!')
    print(shortcut_to_is_variable_star)
    print(is_variable_star == shortcut_to_is_variable_star)

    spectral_type = retrieve_spectral_types(target_list)

    is_variable_star2, spectral_type2 = retrieve_sptype_and_variable(target_list)
    print('Same variable star result:', is_variable_star == is_variable_star2)
    print('Same spectral types:', spectral_type == spectral_type2)

# # HIP ID, Vmag, RA (deg), DEC (deg), AZ (deg), EL (deg)
# hip_id, Vmag, Ra, DEC, AZ, EL = np.loadtxt(path+file, delimiter=',').transpose()
#
# hip_id = hip_id.astype('int')
#
# for id, mag, alt in zip(hip_id, Vmag, EL):
#     print(id, mag, alt)


