import numpy as np

f = 101.6  # oap 1
theta = np.radians(45)

# tolerance is 3.8 arcmins

theta_tolerance = np.radians(3.8/60)

# misalignment of a beam being off-center in the horizontal direction
# is given by angle_error = x_error/(f * np.cos(theta))

x_budget = theta_tolerance * f * np.cos(theta)
print(f'x-axis error budget: {x_budget*1000:.2f} microns')

# ruler alignment accuracy is about 250 um.

theta_error1 = 250e-3/(f * np.cos(theta))
print(f'esimated error upper limit from ruler alignment: {np.degrees(theta_error1)*60:.2f} arcmins')

theta_error_upper = 500e-3/(f * np.cos(theta))  # half a millimeter, way worse than what I measured
print(f'esimated error upper limit from ruler alignment: {np.degrees(theta_error_upper)*60:.2f} arcmins')
# what about output angle?
# I don't have a good measure of output angle
# I could measure the precise distance and the precise size of the projected spot
# from the photos I estimate the projected beam is F/10, exactly what it should be



# what about astimatisim?
# I can use the ration of vertical to horizontal width to estimate the asgimatic ratio, and calculate the angle error

height = 220  # in pixels
width = 290  # pixels
# Note width is complicated by the intensity pattern introduced by the slit aperture
# Note height is shorter than it actually is, due to angle of camera
# this can be used to make a upper limit on error though, as both those effects should make it worse, not better

# projection angle is about 2.5/4
proj_angle = np.arctan(2.5/3)
adj_height = height/np.cos(proj_angle)

R = adj_height/width
theta_error2 = (1 - R)/(1 + R) * 1/np.tan(theta/2)
print(f'estimated error upper limit from astigmatism: {np.degrees(theta_error2)*60:.2f} arcmins')
# that is almost a full degree, that is obviously wrong

