'''
This code is a transcription of Nat's sort-of walkthrough to using spot0
'''
import matplotlib.pyplot as plt
import numpy as np

from spot0 import spot
ps = 18e-6  # pixel pitch
f1, f2 = 0.1016, 0.27224  # OAP focal lengths, meters

# to see documentation, enter help(spot)

# #########################
# #########################
# Things you can do:
# #########################
# #########################

# Get the position of the spot at He-Ne, relative to center on OAP2, looking directly at OAP from the collimation direction:

# position at 1.5-micron (default), step=5 is surface of OAP2
x0,y0,z0,dx0,dy0,dz0,s00 = spot(step=5)
# now specify 632.8 nm
x,y,z,dx,dy,dz,s0 = spot(lam=0.6328e-6,step=5)

# plot x inverted since we're looking opposite to optical axis (z-direction)
plt.plot(-(x[1:]-x0[0]),y[1:]-y0[0],'o')
print("Offset [mm]:",-(x[0]-x0[0])*1.e3) # about 0.92 mm

# Note that the code assumes that the center of the prism is a distance F1 from the surface of OAP2.
# You could also get the angle of deflection after the prism and use that downstream:

x,y,z,dx,dy,dz,s0 = spot(lam=0.6328e-6,step=5)
theta = np.degrees(dx[0])   # about -0.54 degrees
print("Offset [mm]:",-np.radians(theta)*f1*1.e3)  # agrees to about 44 microns, because the distance in the raytrace after the prism to OAP2 is different from the distance to the center of the prism

# What's the offset at the detector from the 1.5 micron ray?
y0,z0,w = spot()
y,z,w = spot(lam=0.6328e-6)
print("Detector Offset [mm]:",(z[0]-z0[0])*1.e3)  # about 2.58 mm (close to but not exactly 0.92*F2/F1; difference is small ~120 microns)

'''You can also find this from the rays leaving OAP2:'''

x,y,z,dx,dy,dz,s0 = spot(lam=0.6328e-6,step=6)
# central ray x[0] will reach focus after path s, such that x[0]+s*dx[0] = 0, so
s = -x[0]/dx[0]   # where
z_detector = z[0]+s*dz[0]   # 2.58 mm as above

'''
This is useful, because it also gives the offset at different distances s along the beam, e.g. corresponding to fold mirrors along the way.  However, what you really
want for these mirrors is the offset relative to the channel center wavelength.  Say we're looking at Channel 1 and define the central wavelength to be 1.65 microns:
'''

x0,y0,z0,dx0,dy0,dz0,s00 = spot(lam=1.65e-6,step=6)
s1 = -x0[0]/dx0[0]
z_detector_center = z0[0]+s1*dz0[0]
print((z_detector-z_detector_center)*1.e3)  # 2.85 mm

# For Channel 2,

x0,y0,z0,dx0,dy0,dz0,s00 = spot(lam=3.e-6,step=6)
s1 = -x0[0]/dx0[0]
z_detector_center = z0[0]+s1*dz0[0]
print((z_detector-z_detector_center)*1.e3)  # 5.82 mm

'''
In this way, you could work out offsets for each fold mirror that allow you to center them using the He-Ne laser with 
offsets that would center the channel central wavelength. Of course, at the detector,
it's a pretty big offset, but it's still smaller than the size of the SP filters (10 mm diameter).
'''

#########################
#########################

'''You can make all the usual spot diagrams.  Note that the meaningful metric is the wavefront error.  Here it is with 
focus set inward of the slit by 250 microns:'''

defoc=250.e-6
y,z,wfe = spot(x00=defoc)
plt.plot ((z[1:]-z[0])/ps,(y[1:]-y[0])/ps,'o',alpha=0.1)  # positive relative to center
plt.xlabel("Dispersion Direction"); plt.ylabel("Cross Dispersion Direction");
# try to correct with an appropriate defocus (moving detector forward)
y,z,wfe = spot(x00=defoc,x11=defoc*(f2/f1)**2)
plt.plot((z[1:]-z[0])/ps,(y[1:]-y[0])/ps,'o',alpha=0.1)
# note that the spot doesn't appear perfect, but the wfe goes to ~zero (so it's perfect)

'''This is essentially the only kind of correction OAP2 can make for OAP1.  Others we saw previously based on the spot 
diagrams (e.g., a shift compensating for a defocus, etc.) are not real.'''

#########################
#########################

'''You can make a plot of the detector position (and image quality) versus wavelength:'''

lam = np.linspace(0.6328,3.5,100)*1.e-6
y,z,wfe = 0.*lam,0.*lam,0.*lam
for i in range(100):
    yy,zz,wfe[i] = spot(lam=lam[i],verbose=False)
    y[i],z[i] = yy[0],zz[0]

plt.plot(z,lam*1e6)  # note that channel2 flips the sign on z (not included here)
# as a result, red end of spectrum is always toward the detector edge (blue toward center)
# like:

#  provided wfe<~0.075 (strehl<~0.8), the spot FWHM ~ 0.95*lam*(Fn*F2/F1)**2
y1 = y + 0.5*0.95*lam*(12*f2/f1)**2
y2 = y - 0.5*0.95*lam*(12*f2/f1)**2
j= lam<2.5e-6
plt.scatter (z[j]-z[j].mean()-512*ps,y1[j],c=lam[j]*1.e6,vmin=lam[0]*1.e6,vmax=lam[-1]*1.e6)
plt.scatter (z[j]-z[j].mean()-512*ps,y2[j],c=lam[j]*1.e6,vmin=lam[0]*1.e6,vmax=lam[-1]*1.e6)
plt.scatter (-z[~j]+z[~j].mean()+512*ps,y1[~j],c=lam[~j]*1.e6,vmin=lam[0]*1.e6,vmax=lam[-1]*1.e6)
plt.scatter (-z[~j]+z[~j].mean()+512*ps,y2[~j],c=lam[~j]*1.e6,vmin=lam[0]*1.e6,vmax=lam[-1]*1.e6)
plt.colorbar(label="Wavelength")
plt.xlabel("Detector Position [Dispersion]"); plt.ylabel("Detector Position [Cross Dispersion]")
plt.xlim((-1024*ps,1024*ps)); plt.ylim((-1024*ps,1024*ps))

plt.plot(lam*1e6,wfe)  #  strehl is ~exp(-(2*pi*wfe)**2)
#  provided wfe<~0.075 (strehl<~0.8), the spot FWHM ~ 0.95*lam*(Fn*F2/F1)**2



