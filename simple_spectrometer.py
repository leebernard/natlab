#########################################
'''
Telescope F#: 11

What is the effective focal length (or F#) of a channel?  Trace a ray at sky
angle theta through the system.

collimator lens F1: F*theta = F1*theta1
camera lens F2: F2*theta1 == Feff*theta

  so,
      Feff = F*F2/F1
  also true,
      F#_eff = F# *F2/F1
      D1 = D*F1/F (size of collimated beam)

( F2/F1 is set by 1.22*lambda0*F#_eff = 2 pixels. )


Now, given prism dispersion dtheta/dlambda, what is R (also sets F1, or D1)?

  F2*(dtheta/dlambda) * dlambda = 1.22*lambda*F#_eff , so

  R = D1*(dtheta/dlambda) / 1.22
    = (size of spectrograph) * (prism dispersion)
'''
#
#
import numpy as np

# How do you get dtheta/dlambda?  For a prism, look up n(lambda),

# wavelength grid, in microns
l = np.linspace(0.8, 3.5, 1000)
l2 = l*l

# see https://refractiveindex.info
# nc=sqrt( 1+0.33973+0.69913*l2/(l2-0.093742**2)+0.11994*l2/(l2-21.182**2)+4.35181*l2/(l2-38.46**2) ) # CaF2
#  some others:
# nb=sqrt(1+2.3330067*l2/(l2-0.0168838419)+0.452961396*l2/(l2-0.0716086325)+1.25172339*l2/(l2-118.707479)) #P-SF68
# ns=sqrt(1+1.4313493*l2/(l2-0.0726631**2)+0.65054713*l2/(l2-0.1193242**2)+5.3414021*l2/(l2-18.028251**2)) # Sapphire
# nz=sqrt( 1+4.45813734*l2/(l2-0.2008598532**2)+0.467216334*l2/(l2-0.3913711662**2)+2.89566290*l2/(l2-47.13621082**2) ) # ZnSe
n = np.sqrt( 1 + 0.6961663*l2/(l2-0.0684043**2)+0.4079426*l2/(l2-0.1162414**2)+0.8974794*l2/(l2-9.896161**2) ) # fused silica

# prism geometry, snell's law: n*sin(0.5*alpha) = sin(0.5*(alpha+theta)) is symmetric passage, can be chosen for one wavelength
#  alpha is wedge angle and theta is deflection angle
alpha = 60*np.pi/180.
i0 = len(l)/2
theta0 = np.arcsin(n[i0]*np.sin(alpha/2))

# trace rays trough prism, theta0 -> theta0p (first surface) and then theta1 -> theta1p (second surface)
#  see, e.g., https://en.wikipedia.org/wiki/Prism
theta0p = np.arcsin(np.sin(theta0)/n) # 1st surface: theta_i is theta0 (index 1), theta_t is theta0p (index n)
theta1 = alpha - theta0p
theta1p = np.arcsin(np.sin(theta1)*n) # 2nd surface: theta_i is theta1 (index n), theta_t is theta1p (index 1)
theta = theta0 + theta1p - alpha

# reflection losses (Fresnel equations), T is total transmission through prism
Rs1 = ( -np.sin(theta0-theta0p)/np.sin(theta0+theta0p) )**2
Rp1 = ( np.tan(theta0-theta0p)/np.tan(theta0+theta0p) )**2
Rs2 = ( -np.sin(theta1p-theta1)/np.sin(theta1+theta1p) )**2
Rp2 = ( np.tan(theta1p-theta1)/np.tan(theta1p+theta1) )**2
T = (1-0.5*(Rs1+Rp1))*(1-0.5*(Rs2+Rp2))

la = (l[1:]+l[:-1])/2.
dl = l[1:]-l[:-1]
dtheta_dl = abs(theta[1:]-theta[:-1])/dl
dtheta1_dl = abs(theta1[1:]-theta1[:-1])/dl

# current design has F1=4 inch collimator, and F2=10.7 inch camera mirror
#   telescope has F-number, Fn=11
F1, F2, Fn = 4*2.54*1.e-2, 10.7*2.54*1.e-2, 11.0

# diffraction limited, using formula 1.028*lambda/D
R = dtheta_dl*F1/(1.028*Fn)*1.e6
D = 0.5
fwhm = 1.028e-6*la/D*180/np.pi*3600   # arcsec

# spread over pixels
ps = 18.e-6  # H2RG pixel size
x = F2*(theta-theta.mean()) / ps
Feff = Fn*D*F2/F1
theta_pixel = ps/Feff*180/np.pi*3600  # arcsec
R0 = R*fwhm/theta_pixel

# note, we can decrease alpha to get better T for the same R, provided F1 increases
# say, alpha=45 degrees, then need F1=0.19 (which is getting big!), and get about 1.5% better transmission
# it's actually D1 = F1/F# which must get bigger (ie, the actual beam size, not just the focal length)



