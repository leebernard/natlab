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

from numpy import arange, zeros, sqrt, sin, cos, round, pi, arctan2, sign, exp

def shear_plate(Rc, R, eps, Rs, alpha=45, plate=-1, T=18.0, wedge_ang=6.35, n_idx=1.46, lam=632.8e-6, atype='coma', wfe=0, wfe_phi=0, N=203):
    """

    Parameters
    ----------
    Rc: radius of curvature of the wavefront
    R: radius of the beam (in the crossplane direction)
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
    if (plate==0): wedge_ang, T = 18.,6.35
    if (plate==1): wedge_ang, T = 40.,2.6
    if (plate==2): wedge_ang, T = 10.,13.
    print("beam radius:", R)
    print("Radius of curvature:", Rc)

    # calculate beam divergence at the edge of the beam
    bdiv = 2*sign(Rc) * arctan2(R, sign(Rc)*Rc)
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
    R0 = R*1.5
    x0 = (arange(N) - (N-1)/2) * R0 * 2/(N-1)
    x = zeros((N, N), dtype='float32') + x0
    y = zeros((N, N), dtype='float32')
    ty = y.T + x0
    print(f"Image pixel scale: {R0*2/(N-1)}")




