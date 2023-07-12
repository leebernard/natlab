import numpy as np
from toolkit import spectrum_slicer

# filename = '/home/lee/natlab/excite_optic_data/CaF2_Transmission.xlsx'
filename = '/home/lee/natlab/excite_optic_data/caf2_trans.dat'

caf2_rawdata = np.loadtxt(filename, skiprows=2)

resolution = 0.010

start_wl = 0.400
end_wl = 5.6
wavelengths = np.arange(start_wl, end_wl, resolution)

# convert to um from nm
caf2_rawdata[:, 0] *= 1/1000

# convert from percentage to fraction
caf2_rawdata[:, 1] *= 1/100


# slice data
temp = spectrum_slicer(wavelengths[0], wavelengths[-1], caf2_rawdata)

transmission = np.interp(wavelengths, temp[:, 0], temp[:, 1])

output_data = np.stack((wavelengths, transmission), axis=1)

outfilename = '/home/lee/natlab/excite_optic_data/caf2_trans_downsampled_v2.dat'

np.savetxt(outfilename, output_data)


