import numpy as np

# wavelength grid, in microns
l = np.linspace(0.8,3.5,1024)
l2 = l*l

# for CaF2 prism: B1=5.676e-1, C1=2.526e-3, B2=4.711e-1, C2=1.008e-2, B3=3.848, C3=1.201e3
a0, a1, a2 = 5.676e-1, 4.711e-1, 3.848
b0, b1, b2 = 2.526e-3, 1.008e-2, 1.201e3
n = np.sqrt( 1 + a0*l2/(l2-b0) + a1*l2/(l2-b1) + a2*l2/(l2-b2) )
dn_dl = 1./n * abs( a0*b0*l/(l2-b0)**2 + a1*b1*l/(l2-b1)**2 + a2*b2*l/(l2-b2)**2 )  # min of ~4.67e-3 at l~1.55 microns

# prism geometry, snell's law: n*sin(0.5*alpha) = sin(0.5*(alpha+theta)) is symmetric passage, can be chosen for one wavelength
#  alpha=60 is wedge angle and theta is deflection angle
reference_l = 1.0  # align to this wavelength (um)
i0 = abs(l-reference_l).argmin()
n0 = n[i0]



