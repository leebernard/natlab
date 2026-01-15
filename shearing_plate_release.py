"""
Module for simulating the interference pattern of a linear shearing interferometer.
This code is based upon the Thorlabs shearing plate interferometer. This plate comes in several sizes, with the
model number and specs given below.

SI254P: 10-25.4 mm beam diameter (R=21.5), ang=18, T=6.35 (plate=0)
SI100P: 5-10 mm beam diameter (R=11.), ang=40, T=2.6 (plate=1)
SI500P: 25.4-50 mm beam diameter (R=30.), ang=10, T=13 (plate=2)
                Note, R is more like 70-80.

This version is intended for release to the public domain as open source software.


"""

from numpy import arange, zeros, ones, sqrt, sin, cos, round, pi, arctan2, sign, exp, ceil, mean, diff

def shear_plate(Rc=1.e5, D=25.4, Dmax=36.0, eps=0, Rs=50, alpha=45, plate=-1, T=6.35, wedge_ang=18., n_idx=1.46, lam=632.8e-6, atype='coma', wfe=0, wfe_phi=0, acenter=[0., 0.], N=203):
    """

    Parameters
    ----------
    Rc: radius of curvature of the wavefront
    D: Diameter of the beam (in the crossplane direction)
    Dmax: Diameter of the shearing plate. This also sets the image size
    eps: fraction of central obscuration, assuming a centered, circular obscuration.
    Rs: Distance from the shearing plate to the observation plane
    alpha: angle of incidence of the beam upon the shearing plate
    plate: option for a preset shearing plate geometry
    T: thickness of the shearing plate
    wedge_ang: wedge angle of the shearing plate
    n_idx: index of refraction of the shearing plate glass
    lam: wavelength of the beam of light
    atype: option for choosing which wavefront aberration to display
    wfe: magnitude of the wavefront error due to aberration
    wfe_phi: angle of the wavefront error in the cross-plane
    N: Parameter controlling resolution of the sample. Number of rays generated is NxN

    Returns
    -------

    """

    # default shearing plate geometry's. This will over-ride passed parameters
    if (plate==0): Dmax, wedge_ang, T = 15., 10., 2.6
    if (plate==1): Dmax, wedge_ang, T = 36., 18., 6.35
    if (plate==2): Dmax, wedge_ang, T = 70., 40., 13.
    print("beam radius:", D)
    print("Radius of curvature:", Rc)

    # calculate beam divergence at the edge of the beam
    Rw = D/2
    bdiv = 2*sign(Rc) * arctan2(Rw, sign(Rc)*Rc)
    print(f"beam divergence: {bdiv: .2e}")

    alpha *= pi/180
    wedge_ang *= pi/180/3600
    shear = T * sin(2*alpha) / sqrt(n_idx**2 - sin(alpha)**2)
    theta = 2*wedge_ang * sqrt(n_idx**2 - sin(alpha)**2)
    print (f"Ray tilt angle: {theta*180/pi*3600:.3f} arcsec")
    print (f"Expected line spacing: {lam/theta:.3f}")
    print (f"Shear distance: {shear:.3f}")
    phi = arctan2(shear,theta*Rc)*180/pi
    print (f"Angle of fringes: {phi:.2f} (degrees)")

    # create some starting ray positions
    D0 = Dmax + int(ceil(shear))
    x0 = (arange(N) - (N-1)/2) * D0 * 1/(N-1)
    # create the primary beam, shifting it so the center of the shearing pattern is zero
    x = zeros((N, N), dtype='float32') + x0 + shear/2
    y = zeros((N, N), dtype='float32')
    yt = y.T; yt += x0
    print (f"Image size: {x.max()-x.min():.2f}x{y.max()-y.min()} (pixel scale: {mean(diff(x)):.3f})")

    # Generate the sheared reflected beam
    xp = x - shear
    yp = y + theta*Rs
    # calculate the combined phase at the observation plane for defocus
    phase = 2*pi/lam*( (theta*(1.-Rs/Rc))*y + (shear/Rc)*x ) - pi*(shear**2 + (theta*Rs)**2) / (lam*Rc)

    # give the rays some abberations
    # if ()

    # calculate the pattern
    # set the outer radius of the reflected beams
    # R1 = Dmax/2
    # R2 = Dmax/2*cos(alpha)
    if (D<Dmax):
        R1 = D/2
    else:
        R1 = Dmax/2

    if (D<Dmax/2*cos(alpha)):
        R2 = D/2
    else:
        R2 = Dmax/2*cos(alpha)
    print("R1, R2:", R1, R2)

    # mask the aperture
    aperture_mask = (x/R1)**2 + (y/R2)**2 > 1
    # mask the central obscuration, if present
    if (eps>0):
        aperture_mask[(x/(R1*eps))**2 + (y/(R2*eps))**2 <= 1] = True
    norm1 = ones(x.shape)
    norm1[aperture_mask] = 0
    print("shape of norm1", norm1.shape)

    # and again for the sheared beam
    aperture_mask = ((x - shear)/R1)**2 + (y/R2)**2 > 1
    if (eps>0):
        aperture_mask[((x - shear)/(R1*eps))**2 + (y/(R2*eps)**2) <= 1] = True
    norm2 = ones(x.shape)
    norm2[aperture_mask] = 0

    # convert phase to intensity
    s = 0.25*(norm1**2 + norm2**2) - 0.5*norm1*norm2*cos(phase)

    return s

