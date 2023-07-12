import numpy as np

# constants
lam = 632.8e-9  # wavelength of light
k = 2*np.pi/lam


def R(z):
    return z


def phi(x, y, z):
    # phase of wavefront
    # x is the separation direction
    # y is the shear direction
    # R (z) is local radius of curvature
    return k/(2*R(z)) * (x**2 + y**2)

def s(t, alpha, n):
    # shear distance between the beams (the spacial displacement)
    # t: plate thickness
    # alpha: angle between principle ray and first surface, about the x-axis (typically, 45 degrees)
    # n: index of refraction
    return t*np.sin(2*alpha)/np.sqrt(n**2 - np.sin(alpha)**2)

def theta(delta, alpha, n):
    # angle of secondary reflected beam ray, about the y-axis. Assumes primary ray angle is 0
    # this angle is caused by the shear plate wedge
    # delta is the angle of the wedge. If delta=0, then theta=0
    # alpha is angle between principle ray and first surface, about the x-axis (typically, 45 degrees)
    # n: index of refraction
    return 2*delta*np.sqrt(n**2 - np.sin(alpha)**2)


def transform_to_second(x, y, z, theta, s):
    x_prime = x - s
    y_z_prime = np.matmul(np.array([[1, theta], [-theta, 1]]), np.array([y, z]))
    return np.stack([x_prime, *y_z_prime])


def path_length_diff(t, alpha, n):
    return 2*t*np.sqrt(n**2 - np.sin(alpha)**2)


def pattern(x, y, z, delta, alpha, n, t):
    primes = transform_to_second(x, y, z, theta(delta, alpha, n),  s(t, alpha, n))
    x_prime = primes[0]
    y_prime = primes[1]
    z_prime = primes[2]
    D = path_length_diff(t, alpha, n)
    return np.sin(0.5*(k*z - k*z_prime - k*D + phi(x, y, z) - phi(x_prime, y_prime, z_prime + D)))**2





