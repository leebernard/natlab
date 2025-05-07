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

from astropy.io import fits

file_path = '/home/lee/nat_lab/excite_mission/test_flight_data/example_fits_files/'

file_complete = '24-08-31_06_26_34_targetName0_excite_complete.fits'

with fits.open(file_path+file_complete) as test_hdul:
    print('Info from the complete sample')
    test_hdul.info()
    # print(test_hdul[0].header)
    # print(test_hdul[1].header)
    # print(test_hdul[2].header)
    data = test_hdul[2]

file_cds = '24-08-31_17_10_27_targetName2_excite_CDS.fits'

with fits.open(file_path+file_cds) as test2_hdul:
    print('And now from the CDS sample')
    test2_hdul.info()




