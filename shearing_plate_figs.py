'''
This script is for generating figures for the shearing plate paper
'''


import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

# play around with oap rotation
from shear_plate_ideal import shear_plate

cmap='twilight'
crop_idn = 25

# ideal case
oap_f = 101.6
oap_F = 4
oap_alpha = np.radians(45)
d_theta = np.radians(0.1)
test_wavelength = 632.8e-6  # HeNe laser wavelength in mm
wf_error = 0  # d_theta*oap_f*np.tan(oap_alpha/2) / (8 * oap_F**2 * np.sqrt(6)) * np.sqrt(1) / test_wavelength

rot_angle = 0.0  # in degrees
ideal_pattern = shear_plate(Rc=1.e10, atype='oap', wfe=wf_error, wfe_phi=rot_angle)

plt.imshow(ideal_pattern.T[crop_idn:-crop_idn, crop_idn:-crop_idn], cmap=cmap)
plt.tight_layout()


# defocus case
wf_error = 0.25
defocus_pattern = shear_plate(Rc=1.e10, atype='defocus', wfe=wf_error, wfe_phi=rot_angle)
plt.imshow(defocus_pattern.T[crop_idn:-crop_idn, crop_idn:-crop_idn], cmap=cmap)
plt.tight_layout()

# oap case
oap_fig, (oap0, oap45, oap90) = plt.subplots(1, 3)
rot_angle = 0
oap_pattern_0 = shear_plate(Rc=1.e10, atype='oap', wfe=wf_error, wfe_phi=rot_angle)
oap0.imshow(oap_pattern_0.T[crop_idn:-crop_idn, crop_idn:-crop_idn], cmap=cmap)
# plt.tight_layout()

rot_angle = 45
oap_pattern_0 = shear_plate(Rc=1.e10, atype='oap', wfe=wf_error, wfe_phi=rot_angle)
oap45.imshow(oap_pattern_0.T[crop_idn:-crop_idn, crop_idn:-crop_idn], cmap=cmap)
# plt.tight_layout()

# oap case
rot_angle = 90
oap_pattern_0 = shear_plate(Rc=1.e10, atype='oap', wfe=wf_error, wfe_phi=rot_angle)
oap90.imshow(oap_pattern_0.T[crop_idn:-crop_idn, crop_idn:-crop_idn], cmap=cmap)
plt.tight_layout()


# Spherical Abberation
rot_angle = 0.0  # in degrees
sphere_pattern = shear_plate(Rc=1.e10, atype='sa', wfe=wf_error, wfe_phi=rot_angle)

plt.imshow(sphere_pattern.T[crop_idn:-crop_idn, crop_idn:-crop_idn], cmap=cmap)
plt.tight_layout()

# coma
rot_angle = 0.0  # in degrees
coma_pattern = shear_plate(Rc=1.e10, atype='coma', wfe=wf_error, wfe_phi=rot_angle)

plt.imshow(coma_pattern.T[crop_idn:-crop_idn, crop_idn:-crop_idn], cmap=cmap)
plt.tight_layout()

