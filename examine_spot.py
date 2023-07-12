import numpy as np
# import matplotlib
import matplotlib.pyplot as plt
import pickle
import os

# import re
from astropy.io import fits
# import pyds9
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from scipy.optimize import curve_fit
# from scipy.special import gammaincc
# import needed functions from the toolbox
# from ccd_tools import bias_subtract, background_subtract, get_regions


filename = '/home/lee/natlab/Zoe/Z 3.14 8.fit'
hdul = fits.open(filename)
imdata = hdul[0].data

darkfilename = '/home/lee/natlab/Zoe/30s_darkframe.fit'
darkframe_hdul = fits.open(darkfilename)
darkframe = darkframe_hdul[0].data

subtracted_data = imdata.astype(np.int32) - darkframe.astype(np.int32)

plt.imshow(subtracted_data, origin='lower', cmap='viridis')
xcenter = 130
ycenter = 524
ywidth = 24
xwidth = 26
subframe = subtracted_data[int(xcenter - xwidth/2):int(xcenter + xwidth/2),
                           int(ycenter - ywidth/2):int(ycenter + ywidth/2)]

plt.imshow(subframe, origin='lower', cmap='viridis')
plt.colorbar()

# save the result
newheader = hdul[0].header
newheader['bitpix'] = 32
# newheader['xorgsubf'] = int(xcenter - xwidth/2)
# newheader['yorgsubf'] = int(ycenter - ywidth/2)

newheader.set('comment', 'background subtracted with 30s dark frame')
newheader.set('comment', 'datatype set to int32, from uint16')


newhdu = fits.PrimaryHDU(data=subtracted_data, header=newheader)

newfilename = 'processed_spotdata.fit'
newhdu.writeto(newfilename)
