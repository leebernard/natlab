import numpy as np

# constants
lam = 632.8e-9  # wavelength of light
k = 2*np.pi/lam

delta_25mm = np.radians(18/3600)  # 18 arcsecs converted to radians

def R(z):
    # plane wave approximation for now
    return z


def phi(x, y, z):
    # phase of wavefront
    # x is the separation direction
    # y is the shear direction
    # R (z) is local radius of curvature
    return k/(2*R(z)) * (x**2 + y**2)

def shear(t, alpha, n):
    '''
    shear distance between the beams (the spacial displacement)
    t: plate thickness
    alpha: angle between principle ray and first surface, about the x-axis (typically, 45 degrees)
    n: index of refraction

    returns
     spacial displacement between the beams
    '''
    return t*np.sin(2*alpha)/np.sqrt(n**2 - np.sin(alpha)**2)

def theta(delta, alpha=np.radians(45), n=1.46):
    '''
    angle of secondary reflected beam ray, about the y-axis. Assumes primary ray angle is 0
    this angle is caused by the shear plate wedge
    delta: the angle of the wedge. If delta=0, then theta=0
    alpha: angle between principle ray and first surface, about the x-axis (typically, 45 degrees)
    n: index of refraction, default n=1.46 fused silica

    returns
    theta: angle of the secondary beam
    '''
    return 2*delta*np.sqrt(n**2 - np.sin(alpha)**2)


def transform_to_second(x, y, z, theta, shear):
    '''
    transform to secondary beam coordinate system
    Parameters
    ----------
    x
    y
    z
    theta: angle of the secondary beam, relative to the primary
    s: path length difference

    Returns
    -------

    '''
    x_prime = x - shear
    y_z_prime = np.matmul(np.array([[1, theta], [-theta, 1]]), np.array([y, z]))
    return np.stack([x_prime, *y_z_prime])


def path_length_diff(t, alpha, n):
    '''
    Parameters
    ----------
    t: plate thickness (meters)
    alpha: angle of shear plate
    n: index of refraction

    Returns
    -------
    Path length difference for the two beams
    '''
    return 2*t*np.sqrt(n**2 - np.sin(alpha)**2)


def pattern(x, y, z, delta, alpha, n, t):
    '''
    Parameters
    ----------
    x
    y
    z
    delta: angle of the wedge
    alpha: angle of shear plate
    n: index of refraction
    t: thickness of plate

    Returns
    -------

    '''
    primes = transform_to_second(x, y, z, theta(delta, alpha, n),  shear(t, alpha, n))
    x_prime = primes[0]
    y_prime = primes[1]
    z_prime = primes[2]
    D = path_length_diff(t, alpha, n)
    return np.sin(0.5*(k*z - k*z_prime - k*D + phi(x, y, z) - phi(x_prime, y_prime, z_prime + D)))**2


shear_angle = theta(delta_25mm)


