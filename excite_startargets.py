import numpy as np
import re
# import csv

debug = True

path = '/home/lee/natlab/excite_targets/'

file = 'simbad_engineering_2023-09-04.tsv'

data = np.loadtxt(path+file, dtype='str', delimiter='\t', skiprows=7, max_rows=21)

data_list = data.tolist()

# re_pattern = r'(V\*)|(Ir\*)|(Er\*)|(Ro\*)|(Pu\*)'
var_flags = r'V\*|Ir\*|Er\*|Ro\*|Pu\*'
bin_flags = r'El\*|EB\*'
var_pattern = re.compile(var_flags)
bin_pattern = re.compile(bin_flags)

debug = False
is_variable = []
is_binary = []
for row in data_list:
    if debug:
        print(row)
    check = bool(re.search(var_pattern, ''.join(row)))
    if debug:
        print('Checking if variable', check)
    is_variable.append(check)

    check = bool(re.search(bin_pattern, ''.join(row)))
    if debug:
        print('Checking if binary variable', check)
    is_binary.append(check)

# clean up some memeory
del data_list
is_variable = np.asarray(is_variable)
is_binary = np.asarray(is_binary)

debug = True

validated_targets = data[~is_variable]

binary_variable = data[is_binary]




# # HIP ID, Vmag, RA (deg), DEC (deg), AZ (deg), EL (deg)
# hip_id, Vmag, Ra, DEC, AZ, EL = np.loadtxt(path+file, delimiter=',').transpose()
#
# hip_id = hip_id.astype('int')
#
# for id, mag, alt in zip(hip_id, Vmag, EL):
#     print(id, mag, alt)


