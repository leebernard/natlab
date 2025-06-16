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

for dir_entry in dir_list:
    if path.splitext(dir_entry)[1] == '.fits':
        print(dir_entry)

# open the first file in the path
test_file = path.join(top_dir, dir_list[0])

with fits.open(test_file) as test_hdul:
    print('Info from the complete sample')
    test_hdul.info()
    # print(test_hdul[0].header)
    # print(test_hdul[1].header)
    # print(test_hdul[2].header)
    data = test_hdul[2]





