"""
This script is meant to test the optical quality of a point source on the FGS camera of the EXCITE telescope.
It does this by finding the profile of the PSF, and characterizing any asymmetry.
Eventually, it will also calculate Strehl ratio.

Current plan is to identify what appears to be astigmatism in the artificial star image.
The image is an elliptical profile, with clear reflection symmetry along perpendicular axis. The axis are rotated by
some angle with respect to the pixel grid.
1. Find the orientation of the FGS relative to the telescope
2. What wavelength are we working at?
3. What is the angle of rotation relative to the dichroic?
"""
import matplotlib
from IPython.core.pylabtools import figsize
from imageio.plugins.spe import roi_array_to_dict

matplotlib.use('TkAgg')


import numpy as np
import matplotlib.pyplot as plt

# consider using photutils in the future
from astropy.io import fits

debug = True

# open the data
file_path = '~/nat_lab/excite_mission/telescope_testing/'
file_name = 'excite_bench_test_20250430.fits'
spot_data = fits.getdata(file_path+file_name)

# calcualte background
source_mask = np.zeros(spot_data.shape)

background = np.ma.masked_array(spot_data, mask=source_mask)
background_median = np.ma.median(background)

# make a postage stamp containing the spot image
spot_stamp = spot_data[798:998,1088:1288]
# spot_stamp = spot_data[828:968,1118:1258]
if debug:
    # display the postage stamp
    fig, ax = plt.subplots()
    ax.imshow(spot_stamp, cmap='viridis', norm='log', vmin=background_median, vmax=200, origin='lower')
    plt.show()

# # plot the histogram
# hist_fig, hist_ax = plt.subplots()
# hist_ax.hist(spot_stamp, bins=100, range=(spot_stamp.min(), 200))
# plt.show()

# manually rotate the image
from scipy.ndimage import rotate
rotated_spot = rotate(spot_stamp, -14)

# pull the center 150 pixels
x_lims = int(rotated_spot.shape[0]/2 - 55), int(rotated_spot.shape[0]/2 + 55)
y_lims = int(rotated_spot.shape[1]/2 - 55), int(rotated_spot.shape[1]/2 + 55)

rotated_stamp = rotated_spot[x_lims[0]:x_lims[1], y_lims[0]:y_lims[1]]
if debug:
    # display the postage stamp
    fig, ax = plt.subplots()
    ax.imshow(rotated_stamp, cmap='viridis', norm='log', vmin=background_median, vmax=200, origin='lower')
    plt.show()


# plot the profiles
stamp_center = int(rotated_stamp.shape[0]/2), int(rotated_stamp.shape[1]/2)
y_slice = rotated_stamp[:,stamp_center[1]] - background_median
x_slice = rotated_stamp[stamp_center[1], :] - background_median

profile_fig, (im_ax, profile_ax) = plt.subplots(2, figsize=(6, 12))

im_ax.imshow(rotated_stamp, cmap='viridis', norm='log', vmin=background_median, vmax=200, origin='lower')

profile_ax.plot(x_slice/x_slice.max(), label='Long direction')
profile_ax.plot(y_slice/y_slice.max(), label='Short direction')
profile_ax.axhline(y=0.5, label='half-max', color='r')
profile_ax.legend()

plt.tight_layout()
plt.show()


