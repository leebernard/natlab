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
# remove the zip files, quick and dirty style
dir_list = dir_list[:-23]

# # open the first file in the path
# test_file = path.join(top_dir, dir_list[500])
#
# test_hdul = fits.open(test_file)
# print('testing hdul:', test_file)
# print('header')
# print(test_hdul[2].header)
# print('data python type', type(test_hdul[2].data))
#
# test_data = test_hdul[2].data.field('ScienceImage') * 1.0  # cast to float
# test_header = test_hdul[2].header
# print(test_data.shape)
# # test_sci_data = test_data.field('ScienceImage') * 1.0  # cast to float
# print(test_data.shape)
# print(test_data.dtype)
# print('number of ramp samples:', test_hdul[2].header['NAXIS2'])
#
# test_signal = test_data[-1] - test_data[0]
# print('mean of last - first:', np.mean(test_signal))
#
# # plot the image frame
# test_fig, (test_sig_ax, first_ax, last_ax) = plt.subplots(3)
# test_sig_ax.imshow(test_signal)
# first_ax.imshow(test_data[2])
# last_ax.imshow(test_data[-1])
# test_fig.tight_layout()


flight_list = dir_list[1096:]  # quick and dirty removeal of the ground data.
i_ascent = [i for i, str in enumerate(flight_list) if 'ascent' in str]
float_list = flight_list[i_ascent[-1]+1:]


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
for dir_entry in float_list:
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

# below has positive means
file_1 = '24-08-31_12_35_18_fullarray0_excite_complete.fits'  # mean 103.0 counts.  Weird artifact along the top
file_2 = '24-08-31_13_08_37_fullarray0_excite_complete.fits'  # mean 283.9 counts. Only two ramp samples
file_3 = '24-08-31_13_45_29_fullarray0_excite_complete.fits'  # mean 526.8.  this one is full of noise. Also shows the readout amps are vertical oriented
file_4 = '24-08-31_14_20_53_targetName2_excite_complete.fits'  # mean 14.31 counts. This one also has the artifact along the top.

'''
After looking at the below 4 files, it's clear they are all blank images, with no signal.
'''
file_5 = '24-08-31_14_34_53_targetName2_excite_complete.fits'  # mean is 1.62 counts
file_6 = '24-08-31_14_37_00_targetName2_excite_complete.fits'  # mean 2.57 counts
file_7 = '24-08-31_14_39_02_targetName2_excite_complete.fits'  # mean 1.15

# the above files were one right after another. Below was separate
file_8 = '24-08-31_15_44_41_targetName2_excite_complete.fits'  # mean 0.478

# open the file
file_path = path.join(top_dir, file_4)

file_hdul = fits.open(file_path)
print('interesting file 2 hdul:',file_path)
print('header')
print(file_hdul[2].header)
print('data python type', type(file_hdul[2].data))

file_data = file_hdul[2].data.field('ScienceImage') * 1.0  # cast to float
file_header = file_hdul[2].header
print(file_data.shape)
# test_sci_data = test_data.field('ScienceImage') * 1.0  # cast to float
print(file_data.shape)
print(file_data.dtype)
print('number of ramp samples:', file_hdul[2].header['NAXIS2'])

file_signal = file_data[-1] - file_data[0]
print('mean of last - first:', np.mean(file_signal))


# plot the image frame
fig, (sig_ax, first_ax, last_ax) = plt.subplots(ncols=3)
sig_ax.imshow(file_signal)
first_ax.imshow(file_data[2])
last_ax.imshow(file_data[-1])
fig.tight_layout()

# zoom in on the frame
zoom_fig, zoom_ax = plt.subplots()
zoom_ax.imshow(file_signal[:,:], norm='log')
zoom_fig.tight_layout()

# make a histogram
hist_fig, hist_ax = plt.subplots()
hist_ax.hist(file_signal.flatten(), bins=281, range=(-300, -20))
hist_fig.tight_layout()


