"""
The purpose of this script is to analyze the test flight data from the Sept 2025 engineering flight of EXCITE

We have data from the science (H2RG), Slit Viewing Camera (SVC), and the Fine Guidance Camera (FGS). We want to evaluate
the performance of the science instrument, or at least check if we got any science signal. The issue is, there were
issues with the readout electronics receiving signal interference, which forced Steve to reset the H2RG periodically.
Each reset saturates the H2RG, ruining the detector data for a certain amount of time. I need to identify which of the
data frames contain useful, unsaturated data.

Steve Maher uploaded some sample data, which can be found at
https://drive.google.com/drive/folders/1KiMbniW6K3ctkcy83xl21KohxAAxbmPe

First step is to examine this sample data, and gain understanding of the headers and data organization.

Once I understand how to access the data, I need to determine which frames are saturated and eliminate them
My idea is to iterate through every single H2RG frame, taking the median of each frame, and look for patters. Hopefully,
the saturated (bad) frames will stand out if I plot these medians as a function of flight time.
"""

# temporary hack to fix backend issues
import matplotlib
matplotlib.use('TkAgg')

import numpy as np
import matplotlib.pyplot as plt
import os.path as path

from os import listdir
from astropy.io import fits

# file_path = '/home/lee/nat_lab/excite_mission/test_flight_data/example_fits_files/'
top_dir = '/home/lee/nat_lab/excite_mission/test_flight_data_ft_sumner24/science_detector_data'

dir_list = listdir(top_dir)

# sort the directory
dir_list.sort()

# open the first file in the path
test_file = path.join(top_dir, dir_list[500])

test_hdul = fits.open(test_file)
print('testing hdul:', test_file)
print('header')
print(test_hdul[2].header)
print('data python type', type(test_hdul[2].data))

test_data = test_hdul[2].data.field('ScienceImage') * 1.0  # cast to float
test_header = test_hdul[2].header
print(test_data.shape)
# test_sci_data = test_data.field('ScienceImage') * 1.0  # cast to float
print(test_data.shape)
print(test_data.dtype)
print('number of ramp samples:', test_hdul[2].header['NAXIS2'])

test_signal = test_data[-1] - test_data[0]
print('mean of last - first:', np.mean(test_signal))

# plot the image frame
test_fig, (test_sig_ax, first_ax, last_ax) = plt.subplots(3)
test_sig_ax.imshow(test_signal)
first_ax.imshow(test_data[2])
last_ax.imshow(test_data[-1])
test_fig.tight_layout()

'''
Data structure appears to be a HDUL with primary header (metadata), SVC HDU, and Science Image HDU
The Science Image HDU is a table in the astropy fitsrec format.
The table has 6 fields. The first field (TTYPE1 = ScienceImage) is a (n, 2048, 768) array. I believe this is a window of
the detector, with 'n' being the number of up-the-ramp samples.
The rest of the frames seem to be housekeeping data.
Of note is TTYPE4 = 'resetFrame', which appears to be a flag for if this frame was taken after the Acadia was reset.
'''
num_ramps = []
sig_mean = []
for dir_entry in dir_list:
    if path.splitext(dir_entry)[1] == '.fits':
        file_path = path.join(top_dir, dir_entry)
        print('file', dir_entry)
        with fits.open(file_path) as hdul:
            science_data = hdul[2].data.field('ScienceImage') * 1.0  # cast to float
            science_header = hdul[2].header
            ramps = science_header['NAXIS2']
            print('number of up-the-ramp reads:', ramps)
            signal = np.mean(science_data[-1] - science_data[0])
            print('Mean of last minus first:', signal)
            num_ramps.append(ramps)
            sig_mean.append(signal)

    else:
        print(dir_entry, 'is not a .fits file')


# interesting files
file_1 = '24-08-31_16_29_52_targetName2_excite_complete.fits'  # has 22 ramp samples and mean -38586.77




