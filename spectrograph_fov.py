"""
The angular tolerance of OAP 1 is ~1 degree. However
"""
import numpy as np
import matplotlib.pyplot as plt


path_length = 101.6*4 + 272.2  # approx internal path length of spectrograph, in mm

mirror_aperture = 11.77  # diameter in mm
mount_depth = 1  # mm
edge_obscuration = mount_depth*np.sin(np.radians(30)) # *np.cos(phi), where phi is the polar angle from the horizontal
effective_mirror_aperture = mirror_aperture - edge_obscuration*2  # effective mirror diameter in mm
print(f'Effect mirror diameter: {effective_mirror_aperture: .2f} mm')
# Final fold mirror is ~11.77 mm in clear dia.
# Dispersion is ~5 mm. centering is ~1 mm. The shadow of the mount is ~0.5 mm on each side.
# 11.77 - 5 - 1 - 1 = 4.8 mm
# Room for +/- 2.4mm of misalignment, 4.8 mm total
error_room_obscured = effective_mirror_aperture - 5 - 1
error_room = mirror_aperture - 5 - 1

# angular field of view.
# assuming the field of view is centered on the final fold mirror, there is 2 mm of space
tolerance = error_room/path_length  # tolerance in radians
fov_vert = np.degrees(tolerance)*60  # fov in arminutes

tolerance_horz = error_room_obscured/path_length
fov_horz = np.degrees(tolerance_horz)*60
print(f'FOV horizontal: {fov_horz: .2f} arcminutes')  # about 24 arcminutes fov. This is not accounting for oap effects.
print(f'FOV vertical: {fov_vert: .2f} arcminutes')

# distance from slit to  dichroic is about 135 mm
d1_distance = 135

print('Centering tolerance on D1:', tolerance*d1_distance)


tiptilt_distance = 300  # estimated distance from slit to tip-tilt mirror

print('Centering tolerance on tip-tilt mirror:', tolerance*tiptilt_distance)


from spot_class import Spot
sp=Spot(N=128)

# distance to mirror
t = 3*25.4e-3  # 3 inches

h=sp.j0
sp.changelambda(0.8e-6)
plt.plot ((sp.z6-t*sp.dz6/sp.dx6)[h],(sp.y6-t*sp.dy6/sp.dx6)[h],'ro',label=r'center slit $\lambda = 0.8 \mu$m')
sp.changelambda(2.5e-6)
plt.plot ((sp.z6-t*sp.dz6/sp.dx6)[h],(sp.y6-t*sp.dy6/sp.dx6)[h],'mo',label='center slit $\lambda = 2.5 \mu$m')

# place source at top of slit
sp.y00 = 1.5e-3
sp.raytrace()

sp.changelambda(0.8e-6)
plt.plot ((sp.z6-t*sp.dz6/sp.dx6)[h],(sp.y6-t*sp.dy6/sp.dx6)[h],'bo',label='top slit $\lambda = 0.8 \mu$m')
sp.changelambda(2.5e-6)
plt.plot ((sp.z6-t*sp.dz6/sp.dx6)[h],(sp.y6-t*sp.dy6/sp.dx6)[h],'go',label='top slit $\lambda = 2.5 \mu$m')

# plot the outside of the mirror
phi=np.linspace(0,2*np.pi,100)
tilt = 20.*np.pi/180
ms = 0.5*25.4e-3 # mirror diameter
plt.plot (ms/2*np.sin(phi)*np.cos(tilt),ms/2*np.cos(phi),'k:',label='mirror edge')
plt.xlim((-ms/2,ms/2)); plt.ylim((-ms/2,ms/2))
plt.xlabel("Position [mm]"); plt.ylabel("Position [mm]")
plt.legend()


