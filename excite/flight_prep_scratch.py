"""
This is scratch space for Ft Sumner 2024 flight campaign calculations.
"""

import numpy as np

F = 12
D = 0.5  # m
telescope_fl = F * D  # m

square_accuracy = np.degrees(0.178/250) * 60  # arcminutes

tube_accuracy = np.degrees(0.8/300) * 60 - 2  # arcminutes

basler_dim = np.array((2448, 2048))  # width x height of detector
basler_pxpitch = 3.45  # um

dot_loc = np.array((1400, 1100))

camera_center = basler_dim/2

centering_err = np.absolute(camera_center - dot_loc) * basler_pxpitch

basler_xdim = basler_dim[0] * basler_pxpitch
basler_ydim = basler_dim[1] * basler_pxpitch

basler_telescope_fov = np.degrees(basler_dim*basler_pxpitch*1e-6/telescope_fl) * 60  # arcminutes


as_alignment = np.sqrt((np.degrees(centering_err*1e-6/telescope_fl) * 60)**2 + tube_accuracy**2)  # arcminutes
print(f'artifical star alignment: {as_alignment[0]:.2f} arcminutes')

tip_tilt_focal_distance = 330 # mm, approx to about +/- 10 mm

fgs_tiptilt_fov = np.degrees(basler_xdim*1e-3/tip_tilt_focal_distance)

print(f'tip-tilt fov: {fgs_tiptilt_fov:.3f} degrees')


D_artstar = 6.5*25.4  # mm
wl = .8  # wavelength in um
spot_artstar = wl/(D_artstar*1e3) * telescope_fl*1e6  # spot size in um
print(f'artifical star spot size: {spot_artstar:.2f} um')
print(f'In pixels: {spot_artstar/basler_pxpitch: .2f}')

delta_y = 500*basler_pxpitch  # um
delta_angle = np.degrees(delta_y*1e-3/tip_tilt_focal_distance)
print(f'change in angle after adjusting dichroic: {delta_angle*60:.2f} arcminutes')



