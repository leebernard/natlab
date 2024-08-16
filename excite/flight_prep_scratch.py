"""
This is scratch space for Ft Sumner 2024 flight campaign calculations.

slit temperature 24-08-14 14:51
177 K


"""
# import matplotlib
# matplotlib.use("qt5agg")
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


D_artstar = np.array((6.5, 8))*25.4  # mm
wl = .8  # wavelength in um
spot_artstar = wl/(D_artstar*1e3) * telescope_fl*1e6  # spot size in um
print(f'artifical star diffraction fwhm: ({spot_artstar[0]:.2f}, {spot_artstar[1]:.2f}) um')
print(f'In pixels: ({spot_artstar[0]/basler_pxpitch: .2f}, {spot_artstar[1]/basler_pxpitch: .2f})')

# size of the spot as measured
fgc_spot_fwhm = np.array((11, 15))  # in pixels
print(f'measured spot on FGC: ({fgc_spot_fwhm[0]*basler_pxpitch: .2f}, {fgc_spot_fwhm[1]*basler_pxpitch: .2f}) um')


delta_y = 500*basler_pxpitch  # um
delta_angle = np.degrees(delta_y*1e-3/tip_tilt_focal_distance)
print(f'change in angle after adjusting dichroic: {delta_angle*60:.2f} arcminutes')

'''offset due to dichroic'''
fgc_y = 1515
offset = basler_dim[1]/2 - fgc_y
print(f'vertical offset is {offset*basler_pxpitch / 25.4 /1000: .4f} inches')

dichroic_focal_distance = 135  # mm
n = 1.46
wedge_angle = np.radians(25.9/60)  # converting arcminutes to radians
dichroic_angle = np.radians(44.5)
wedge_thickness = 10  # mm

# fgc offset calculations
sgc_angle_offset = np.arcsin(n*np.sin(np.arcsin(np.sin(dichroic_angle/n - wedge_angle))) - dichroic_angle)
print(f'angle of fgc beam: {np.degrees(sgc_angle_offset): .3f} degrees')
print(f'physical displacement at focal plane: {sgc_angle_offset * dichroic_focal_distance:.3f} mm')


thickness_offset = np.sin(dichroic_angle)/n * wedge_thickness

# slit ghost offset calculations
angle_offset = wedge_angle * dichroic_focal_distance
ghost_offset = angle_offset + thickness_offset
print(f'offset of ghost: {ghost_offset: .2f} mm')

# slit width um/slit width px
svc_pxscale = 100/15.7
svc_fwhm_upper = 50  # px
print(f'spot size on slit substrate (estimated FWHM) Upper limit: {svc_fwhm_upper * svc_pxscale: .2f} um')
svc_fwhm_lower = 25  # px
print(f'spot size on slit substrate (estimated FWHM) Lower limit: {svc_fwhm_lower * svc_pxscale: .2f} um')


# esimated worst case of the spot size
# about twice the size of the slit
# slit_spot_size = 200 # um
# that's about 6x the diffraction limit
defocus = svc_fwhm_upper * svc_pxscale * 12  # um
print(f'Worst case, slit is defocused by {defocus/1000:.2f} mm')




