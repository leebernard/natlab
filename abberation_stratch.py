"""
This script is to support learning diffraction theory of optical aberrations.
It is largely based on Basic Wavefront Aberration Theory for Optical Metrology, by Kingslake and Shannon 1992

Sign conventions:
+z is the direction of propagation.
+x is the vertical direction (meridional)
+y is the right-hand direction relative to the direction of propagation (sagittal)

Note: Born and Wolf uses y as the meridional and x as the sagittal

Optical Path Difference (OPD):
Positive if the aberrated wavefront leads the ideal wavefront
Negative if the aberrated wavefront trails the ideal wavefront
Positive if the wavefront curvature is tighter than ideal
    E.g., a negative focal shift (tighter curvature, image displaced in the -z direction) produces a positive aberration

Angles increase in the CCW direction, relative to the z-axis
The zero angle in the pupil plane is defined as the x-axis

"""

import numpy as np
import matplotlib
matplotlib.use("qt5agg")
import matplotlib.pyplot as plt

from scipy import integrate
# import bessel function
from scipy.special import jv

# perfect lens with beam diameter 8 mm, focal length 101.4 mm

d = 8  # mm
f = 101.4 # mm

# define wavefront in the pupil
def w_perfect(x, y, R):
    return (x**2 + y**2) / (2*R)


def w_defocus(x, y, R, epsilon_z):
    return (x**2 + y**2)/(2*R) - epsilon_z * (x**2 + y**2)/(2*R)


def phi_defocus(roh, m):
    '''
    Phase error due to defocus, defined in number of wavelengths m.
    Note that defocus is rotationally symmetric, so this phase difference is a function of radius roh.

    Parameters
    ----------
    m: the focus error, in units of wavelengths
    roh: radial position

    Returns
    -------
    The difference in phase from an ideal wavefront
    '''

    return m*2*np.pi*roh**2


def pupil_defocus(roh):
    return jv(0, roh)


def aperture_wavefront(roh, pupil_func, r, phi, m):
    return pupil_func(np.pi*r*roh) * np.exp(complex(0, phi(roh, m))) * 2*roh

defocus = 0
r = 0
intensity, _ = integrate.quad(aperture_wavefront, 0, 1, args=(pupil_defocus, r, phi_defocus, defocus), complex_func=True)


r = np.linspace(0, 16)
intensity_dist = np.zeros(r.shape, dtype=np.cdouble)
for idr, dr in enumerate(r):
    intensity_dist[idr], _ = integrate.quad(aperture_wavefront, 0, 1, args=(pupil_defocus, dr, phi_defocus, defocus), complex_func=True)

irradiance = np.abs(intensity_dist)**2

fig, ax = plt.subplots()
ax.plot(r, irradiance)

plt.show()


