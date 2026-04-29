"This script is for examining the science detector data from the first light"

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

from zscale import zscale

file_name = '/home/lee/nat_lab/excite_mission/goddard_test_data/Dev_Obs21757_Exp54_firstlight_complete.fits'

test_hdul = fits.open(file_name)

test_data = test_hdul[2].data.field('ScienceImage') * 1.0  # cast to float

test_frame = test_data[-5] - test_data[0]

z_min, z_max = zscale(test_frame)

plt.imshow(test_frame, norm='log')

slit_h = 3.0  # mm
f1 = 101.6  # mm
f2 = 272.25  # mm
m = f2/f1  # mm
px_pitch = 18.0e-3  # mm
n_slit_px = slit_h * m / px_pitch

slit_top = 1157
slit_bot = 1489
slit_showing = slit_bot - slit_top

fract_slit_showing = slit_showing/n_slit_px

apt_meas_r = np.array([245, 225, 258, 237])  # pixels
aperture_r = np.mean(apt_meas_r)

aperture_dia = 2*aperture_r  # *px_pitch

delta_x = (aperture_dia - n_slit_px)/2

shift_needed = 82 * px_pitch  # mm

shift_angle = 82*px_pitch/(f2)


