"""
The angular tolerance of OAP 1 is ~1 degree. However
"""
import numpy as np

path_length = 101.6*4 + 272.2  # approx internal path length of spectrograph

# Final fold mirror is ~10 mm in clear dia.
# Dispersion is ~5 mm. centering is ~1 mm. The shadow of the mount is ~1 mm on each side.
# 10 - 5 - 1 - 2 = 2 mm
# Room for +/- 1mm of misalignment, 2 mm total

# angular field of view.
# assuming the field of view is centered on the final fold mirror, there is 2 mm of space
tolerance = 2/path_length  # tolerance in radians
fov = np.degrees(tolerance)*60  # fov in arminutes
print(fov, 'arcminutes')  # about 10 arcminutes fov. This is not accounting for oap effects.

# distance from slit to D1 dichroic is about 135 mm
d1_distance = 135

print('Centering tolerance on D1:', tolerance*d1_distance)


tiptilt_distance = 300  # estimated distance from slit to tip-tilt mirror

print('Centering tolerance on tip-tilt mirror:', tolerance*tiptilt_distance)



