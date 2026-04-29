import numpy as np

f = 101.6  # oap 1
alpha = np.radians(45)

# tolerance is 3.8 arcmins

alpha_tolerance = np.radians(3.8/60)

# misalignment of a beam being off-center in the horizontal direction
# is given by angle_error = x_error/(f * np.cos(alpha))

x_budget = alpha_tolerance * f * np.cos(alpha)
print(f'x-axis error budget: {x_budget*1000:.2f} microns')

# ruler alignment accuracy is about 250 um.
# alignments made with the caliper had an accuracy of about 50 um.

alpha_error1 = 250e-3/(f * np.cos(alpha))
print(f'esimated error upper limit from ruler alignment: {np.degrees(alpha_error1)*60:.2f} arcmins')

alpha_error_upper = 500e-3/(f * np.cos(alpha))  # half a millimeter, way worse than what I measured
print(f'esimated error upper limit from ruler alignment: {np.degrees(alpha_error_upper)*60:.2f} arcmins')
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
alpha_error2 = (1 - R)/(1 + R) * 1/np.tan(alpha/2)
print(f'estimated error upper limit from astigmatism: {np.degrees(alpha_error2)*60:.2f} arcmins')
# that is almost a full degree, that is obviously wrong

'''Strehl ratio calculations'''

eps = .38  # fraction of aperture obscured by secondary
crit_wavelength = 1.0e-6  # criteria wavelength in meters
F = 12  # approximate F number of the telescope beam
f = 101.6e-3  # focal length in meters
working_angle = np.radians(45)

defocus_tol = 2.41 * crit_wavelength * F**2
print(f'Defocus tolerance: {defocus_tol*1e6:.2f} um')

angular_tol = (defocus_tol - 100e-6)/(f * np.tan(working_angle/2))
print(f'Angular tolerance: {np.degrees(angular_tol)*60:.2f} arcminutes')

lambda_req = 1.0  # um
lambda_goal = .633 # um
F_req = 12
F_goal = 4
f1 = 101.6  # mm
f2 = 272.25 # mm
mag = f2/f1

err1_focus_req = 8*np.sqrt(12)/(1-.38**2) * .075*lambda_req * F_req**2  # um
err1_focus_goal = 8*np.sqrt(12)/(1-0.0**2) * .075*lambda_goal * F_goal**2  # um

err1_alpha_req = 8*np.sqrt(6)/np.tan(alpha/2) /np.sqrt(1 + .38**2 + .38**4) * 0.075*lambda_req * F_req**2 * 1.0e-3/ f1  # radians
err1_alpha_goal = 8*np.sqrt(6)/np.tan(alpha/2) /np.sqrt(1 + 0.0**2 + 0.0**4) * 0.075*lambda_goal * F_goal**2 * 1.0e-3/ f1  # radians

print('OAP 1 angular error requirement/goal:', np.degrees(err1_alpha_req)*60, '/', np.degrees(err1_alpha_goal)*60, 'arcmins')
print('OAP 1 focus error requirement/goal:', err1_focus_req, '/', err1_focus_goal, 'um')
# now add in the 2.7 factor for OAP 2
# focal length is 272.25 mm


alpha2 = np.radians(30)
err2_focus_req = 8*np.sqrt(12)/(1-.38**2) * .075*lambda_req * (F_req*mag)**2  # um
err2_focus_goal = 8*np.sqrt(12)/(1-0.0**2) * .075*lambda_goal * (F_goal*mag)**2  # um

err2_alpha_req = 8*np.sqrt(6)/np.tan(alpha2/2) /np.sqrt(1 + .38**2 + .38**4) * 0.075*lambda_req * (F_req*mag)**2 * 1.0e-3/ f2  # radians
err2_alpha_goal = 8*np.sqrt(6)/np.tan(alpha2/2) /np.sqrt(1 + 0.0**2 + 0.0**4) * 0.075*lambda_goal * (F_goal*mag)**2 * 1.0e-3/ f2  # radians

print('OAP 2 angular error requirement/goal:', np.degrees(err2_alpha_req)*60, '/', np.degrees(err2_alpha_goal)*60, 'arcmins')
print('OAP 2 focus error requirement/goal:', err2_focus_req, '/', err2_focus_goal, 'um')


# calculate the measured error
meas_alpha_err = np.radians(16/60) # error is 16 arcmins

alpha_wrms = meas_alpha_err * np.tan(alpha/2)/(8*np.sqrt(6)) * np.sqrt(1 + .38**2 + .38**4) * 101.6/12**2 * 2 * 1000

focus_wrms_budget = np.sqrt(0.075**2 - alpha_wrms**2)

err1_focus_budget = 8*np.sqrt(12)/(1-.38**2) * focus_wrms_budget * 12**2

np.sqrt(np.log(1/0.8))/(2*np.pi)



