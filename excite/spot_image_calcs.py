import numpy as np

F = 12  # f number of telescope
D = 0.5  # diameter of primary mirror, meters

wl = .800  # wavelength, um

dof = 2 * wl * F**2  # depth of focus, in units of whatever wl is
print(f'Excite depth of focus: +/-{dof: .2f} um')





