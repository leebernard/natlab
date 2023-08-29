"""
This is for a mapping function between the Fine Guidance System (FGS) and the Slit Viewing Camera (SVC).
This should include the mapping function, and a procedure for calibrating the mapping funciton
"""

import numpy as np

from skimage.transform import rescale

ps_fgs = 3.15  # plate scale of the fgs
ps_svc = 8.5  # plate scale of the svc

# open data here
from astropy.io import fits
from astropy.utils.data import download_file
image_file = download_file('http://data.astropy.org/tutorials/FITS-images/HorseHead.fits', cache=True )
fgs_im = fits.getdata(image_file)

# fgs is right angle to the svc
# rotate clockwise 90 degrees
# rotate and downscale the fgs image to the scv scale
scale = ps_fgs/ps_svc
sub_im = rescale(np.rot90(fgs_im, axes=(1, 0)), scale, preserve_range=True, anti_aliasing=True)

# pad the downscaled image to the size of the svc image, with zeros
# upper left corner
offset = (307, 423)  # y,x
# use hstack and vstack to pad out the image
sub_im = np.vstack

'''
svc_im = np.insert(np.zeros(im_size),  # generate an array of zeros
                       int(im_size[0]/2),                               # location to insert data, the middle
                       ,                                       # spectrum to be inserted
                       axis=0)                                             # axis the spectrum is inserted along
'''


