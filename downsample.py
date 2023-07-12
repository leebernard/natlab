import numpy as np

from toolkit import spectrum_slicer

wavelength_range = [400, 5600]  # nm

path = '/home/lee/natlab/excite_optic_data/'
file = 'caf2_trans.dat'

data = np.loadtxt(path + file, skiprows=2)

sliced_data = spectrum_slicer(*wavelength_range, data)
