from numpy import arange,zeros,sqrt,dot,newaxis,sin,cos,array,round,pi,arctan2,abs,sign,exp,where

def shear_plate(Rc=1.e5,R=21.5,Rs=50.,T=6.35,wedge_ang=18.,N=203,n_idx=1.46,lam=632.8e-6,wfe=0.,wfe_phi=0.,w0=21.1,
                alpha=45.,atype='coma',nripple=10,acenter=[0.,0.],plate=-1):
    """
      SI254P: 10-25.4 mm beam diameter (R=21.5), ang=18, T=6.35 (plate=0)
      SI100P: 5-10 mm beam diameter (R=11.), ang=40, T=2.6 (plate=1)
      alpha: working angle for shear plate
      Rc is radius of curvature (distance from shear plate to focus)
      R is beam radius
      Rs is distance from shear plate to observation screen
      T is shear plate thickness
      n_idx is index of refraction
      wedge_ang is wedge angle in arcsec
      lam: wavelength in same units as R,T
      N: create NxN samples
      wfe: rms wavefront error (units of lam)

      Then x = shear_plate(); imshow(x.T,cmap='gray')
      lab setup uses f1=2.75, f2=40; 0.63 mm beam diameter (1/e^2); FWHM~0.74 mm
           (lens f1 focus is 1.91 mm from back of lens holder)
           (lens f2 focus is 31.0 mm from back of lens holder)
           -> 2.76 micron spot; collimated beam size: 9.16 mm
              focus is with f/12 using a 110 mm focal length lens
         pinhole diameter should be ~ 1.3*632.8e-9/0.63e-3*2.75e-3 = 3.6 micron
              25 mm portion?
              something at focus ~1mm in front of slit

        want a f1=1.49 mm lens
        https://www.thorlabs.com/thorproduct.cfm?partnumber=C710TMD-B 
        adapter: https://www.thorlabs.com/thorproduct.cfm?partnumber=E06RMS
         then, pinhole diameter should be ~ 1.3*632.8e-9/0.63e-3*1.49e-3 = 1.95 micron (use 2-micron pinhole)
            can then get beamsize 9.7-16.2 mm
    """
    if (plate==0): R,wedge_ang,T = 21.5,18.,6.35
    if (plate==1): R,wedge_ang,T = 11.,40.,2.6

    if (w0>0): bdiv = sign(Rc)*arctan2(w0/2,sign(Rc)*Rc)
    else: bdiv = sign(Rc)*arctan2(R,sign(Rc)*Rc)
    print (f"Beam divergence at edge: {bdiv:.2e}")

    alpha *= pi/180.  # convert to radians
    shear = T*sin(2*alpha)/sqrt(n_idx**2-sin(alpha)**2)
    th = 2*wedge_ang*pi/180/3600*sqrt(n_idx**2-sin(alpha)**2)
    print (f"Expected line spacing: {lam/th:.3f}")
    print (f"Shear distance: {shear:.3f}")

    # create some starting ray positions
    if (w0>0): R0=w0*1.5
    else: R0=R
    x0 = (arange(N)-(N-1)/2)*R0*2/(N-1)
    x = zeros((N,N),dtype='float32') + x0
    y = zeros((N,N),dtype='float32'); yt=y.T; yt += x0

    x += shear/2
    if (w0>0): norm = exp(-4*(x/w0)**2-4*(y/w0)**2)
    xp = x - shear
    yp = y + th*Rs
    #phase = 2*pi/lam*( th*y + 0.5/Rc*(x**2+y**2) - 0.5/Rc*(xp**2+yp**2) ), expand and drop constant term:
    phase = 2*pi/lam*( (th*(1.-Rs/Rc))*y + (shear/Rc)*x )

    # give the rays some abberations
    if (wfe>0):
        if (w0>0): Rw = w0/2.
        else: Rw = R
        c,s = cos(wfe_phi*pi/180),sin(wfe_phi*pi/180)
        x1,y1 = (x-acenter[0])*c + (y-acenter[1])*s, -(x-acenter[0])*s + (y-acenter[1])*c
        xp1,yp1 = (xp-acenter[0])*c + (yp-acenter[1])*s, -(xp-acenter[0])*s + (yp-acenter[1])*c
        if (atype=='coma'):
            W = 2**1.5*wfe*lam*(x1**2+y1**2)*x1/Rw**3 # coma
            W -= 2**1.5*wfe*lam*(xp1**2+yp1**2)*xp1/Rw**3
        elif (atype=='ripple'):
            W = sqrt(2)*wfe*lam*sin(pi/2*nripple/Rw*x1) # ripple
            W -= sqrt(2)*wfe*lam*sin(pi/2*nripple/Rw*xp1)
        elif (atype=='ripple0'):
            W = sqrt(2)*wfe*lam*sin(pi/2*nripple/Rw*(sqrt(x1**2+y1**2))) # ripple0
            W -= sqrt(2)*wfe*lam*sin(pi/2*nripple/Rw*(sqrt(xp1**2+yp1**2)))
        elif (atype=='coma5'):
            W = 2**1.5*wfe*lam*(x1**2+y1**2)*x1**3/Rw**5 # coma5
            W -= 2**1.5*wfe*lam*(xp1**2+yp1**2)*xp1**3/Rw**5
        elif (atype=='comatilt'):
            W = 2**1.5*wfe*lam*( 3*(x1**2+y1**2)/Rw**2 - 2 )*x1/Rw # comatilt
            W -= 2**1.5*wfe*lam*( 3*(xp1**2+yp1**2)/Rw**2 - 2 )*xp1/Rw
        elif (atype=='asti'):
            W = 4*wfe*lam*x1**2/Rw**2 # astigmatism
            W -= 4*wfe*lam*xp1**2/Rw**2
        elif (atype=='tilt'):
            W = 2*wfe*lam*x1/Rw # tilt
            W -= 2*wfe*lam*xp1/Rw
        elif (atype=='sa'):
            W = 1.5*sqrt(5)*wfe*lam*(x1**2+y1**2)**2/Rw**4 # SA
            W -= 1.5*sqrt(5)*wfe*lam*(xp1**2+yp1**2)**2/Rw**4
        else:
            W = 0.

        phase += 2*pi/lam*W

    R1,R2 = R/sqrt(2),R
    good = (x/R1)**2+(y/R2)**2<=1
    nx = int(round(shear/(x0[1]-x0[0])))
 
    if (w0>0):
        norm1 = 0.*phase
        norm1[:,nx:N+nx] = norm[:,:N-nx]
        norm2 = 1.*norm
    else:
        norm1,norm2=1.,1.
 
    s = 0.5*(norm1**2+norm2**2) + norm1*norm2*cos(phase)

    mask = 1.*good
    mask[:,nx:N+nx] += good[:,:N-nx]
    s[mask<2]=0

    return s.T
