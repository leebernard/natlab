from numpy import arange,zeros,sqrt,dot,newaxis,sin,cos,array,round,pi,arctan2,abs,sign,exp,where,logical_and,ceil
from numpy.random import randn

def shear_plate(Rc=1.e5,D=25.4,Dmax=36.,eps=0.,Rs=50.,T=6.35,wedge_ang=18.,N=203,n_idx=1.46,lam=632.8e-6,wfe=0.,wfe_phi=0.,w0=0.,alpha=45.,atype='coma',nripple=10,acenter=[0.,0.],plate=-1,visual=True):
    """
      alpha: working angle for shear plate
      Rc is radius of curvature (distance from shear plate to focus)
      D is beam diameter
      w0 (if>0) is gaussian 1/e^2 beam width (overrides D)
      eps*D is diameter of central obscuration (default 0)
      Rs is distance from shear plate to observation screen
      Dmax is shear plate aperture (set's image size also)
      T is shear plate thickness
      n_idx is index of refraction
      wedge_ang is wedge angle in arcsec
      lam: wavelength in same units as R,T
      N: image size (NxN samples)
      visual: if True, make like the pattern on ground glass

      plate: pick a shear plate (default -1, otherwise overrides Dmax,wedge_ang,T)
      # implied units are mm
      SI100P: 5-10 mm beam diameter (Dmax=15.), ang=40, T=2.6 (plate=0)
      SI254P: 10-25.4 mm beam diameter (Dmax=36.), ang=18, T=6.35 (plate=1, default)
      SI500P: 25.4-50 mm beam diameter (Dmax=70.), ang=10, T=13 (plate=2)
      
      # adding wavefront errors:
      atype: type of wavefront error (string)
      wfe: rms wavefront error (units of lam)
      wfe_phi: rotation angle of wavefront error (default 0)
      acenter: offset from center of wavefront error (default [0,0])

      Then x = shear_plate(visual=False); imshow(x,cmap='gray'), or
           x = shear_plate(); imshow(x)
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
    if (plate==0): Dmax,wedge_ang,T = 15.,40.,2.6
    if (plate==1): Dmax,wedge_ang,T = 36.,18.,6.35
    if (plate==2): Dmax,wedge_ang,T = 70.,10.,13.
    if (w0>0): print ("Beam diameter:",w0,"(1/e^2)")
    else: print ("Beam diameter:",D)

    print ("Radius of curvature:",Rc)

    if (w0>0): Rw = w0/2.
    else: Rw = D/2

    Wrms = 1/(2*Rc*sqrt(12)*lam)*Rw**2*(1-eps**2)
    print (f"Rms wavefront error associated with Rc: {Wrms:.3f}/lambda")

    if (w0>0): bdiv = sign(Rc)*arctan2(w0/2,sign(Rc)*Rc)
    else: bdiv = 2*sign(Rc)*arctan2(D/2,sign(Rc)*Rc)
    print (f"Beam divergence: {bdiv:.2e}")

    alpha *= pi/180.
    shear = T*sin(2*alpha)/sqrt(n_idx**2-sin(alpha)**2)
    th = 2*wedge_ang*pi/180/3600*sqrt(n_idx**2-sin(alpha)**2)
    print (f"Ray tilt angle: {th*180/pi*3600:.2f} arcsec")
    print (f"Fringe spacing in y: {lam/th/(1-Rs/Rc):.2f}")
    print (f"Shear distance: {shear:.2f}")
    phi = arctan2(shear,th*Rc)*180/pi
    print (f"Rotation angle: {phi:.2f} (degrees)")

    # create some starting ray positions
    #D0 = Dmax/cos(alpha)
    D0 = Dmax + int(ceil(shear))
    x0 = (arange(N)-(N-1)/2)*D0/(N-1)
    x = zeros((N,N),dtype='float32') + x0 + shear/2
    y = zeros((N,N),dtype='float32'); yt=y.T; yt += x0
    print (f"Image size: {D0*N/(N-1):.2f} (pixel scale: {D0/(N-1):.3f})")

    xp, yp = x - shear, y + th*Rs
    #phase = 2*pi/lam*( th*y + 0.5/Rc*(x**2+y**2) - 0.5/Rc*(xp**2+yp**2) ), expand:
    phase = 2*pi/lam*( (th*(1.-Rs/Rc))*y + (shear/Rc)*x ) - pi*(shear**2+(th*Rs)**2)/(lam*Rc)

    # give the rays some abberations
    if (wfe>0):
        c,s = cos(wfe_phi*pi/180),sin(wfe_phi*pi/180)
        x1,y1 = (x-acenter[0])*c + (y-acenter[1])*s, -(x-acenter[0])*s + (y-acenter[1])*c
        xp1,yp1 = (xp-acenter[0])*c + (yp-acenter[1])*s, -(xp-acenter[0])*s + (yp-acenter[1])*c
        if (atype=='coma'):
            eps_fac = sqrt( 1 + eps**2 + eps**4 + eps**6 )
            W = 2**1.5*wfe*lam*(x1**2+y1**2)*x1/Rw**3/eps_fac # coma
            W -= 2**1.5*wfe*lam*(xp1**2+yp1**2)*xp1/Rw**3/eps_fac
        elif (atype=='defocus'):
            eps_fac = 1-eps**2
            W = 2*sqrt(3)*wfe*lam*(x1**2+y1**2)/Rw**2/eps_fac # defocus
            W -= 2*sqrt(3)*wfe*lam*(xp1**2+yp1**2)/Rw**2/eps_fac
        elif (atype=='ripple'):
            W = sqrt(2)*wfe*lam*sin(pi/2*nripple/Rw*x1) # ripple
            W -= sqrt(2)*wfe*lam*sin(pi/2*nripple/Rw*xp1)
        elif (atype=='ripple0'):
            W = sqrt(2)*wfe*lam*sin(pi/2*nripple/Rw*(sqrt(x1**2+y1**2))) # ripple0
            W -= sqrt(2)*wfe*lam*sin(pi/2*nripple/Rw*(sqrt(xp1**2+yp1**2)))
        elif (atype=='coma5'):
            eps_fac = sqrt( 1 + eps**2 + eps**4 + eps**6 + eps**8 + eps**10 )
            W = 4*sqrt(6/5)*wfe*lam*(x1**2+y1**2)*x1**3/Rw**5/eps_fac # coma5
            W -= 4*sqrt(6/5)*wfe*lam*(xp1**2+yp1**2)*xp1**3/Rw**5/eps_fac
        elif (atype=='tilt'):
            eps_fac = sqrt(1+eps**2)
            W = 2*wfe*lam*x1/Rw/eps_fac # comatilt
            W -= 2*wfe*lam*xp1/Rw/eps_fac
        elif (atype=='comatilt'):
            eps_fac = sqrt( 1 + eps**2 - 7*eps**4 + 9*eps**6 )
            W = 2**1.5*wfe*lam*( 3*(x1**2+y1**2)/Rw**2 - 2 )*x1/Rw/eps_fac # comatilt
            W -= 2**1.5*wfe*lam*( 3*(xp1**2+yp1**2)/Rw**2 - 2 )*xp1/Rw/eps_fac
        elif (atype=='oap'):
            eps_fac = sqrt(1+eps**2+eps**4)
            W = sqrt(6)*wfe*lam*(x1**2-y1**2)/Rw**2/eps_fac # potato-chipping-like astigmatism
            W -= sqrt(6)*wfe*lam*(xp1**2-yp1**2)/Rw**2/eps_fac
        elif (atype=='asti'):
            eps_fac = sqrt(1+eps**4)
            W = 4*wfe*lam*x1**2/Rw**2/eps_fac # astigmatism
            W -= 4*wfe*lam*xp1**2/Rw**2/eps_fac
        elif (atype=='sa'):
            eps_fac = (1-eps**2)*sqrt( 1 + eps**2*(eps**2+7/4) )
            W = 1.5*sqrt(5)*wfe*lam*(x1**2+y1**2)**2/Rw**4/eps_fac # SA
            W -= 1.5*sqrt(5)*wfe*lam*(xp1**2+yp1**2)**2/Rw**4/eps_fac
        else:
            W = 0.

        phase += 2*pi/lam*W

    # work out edges and normalization
    R1,R2 = Dmax/2,Dmax/2*cos(alpha)   #  if R1,R2 set by masking at shear plate, otherwise:
    if (D/2<R1 and w0<=0): R1=D/2
    if (D/2<R2 and w0<=0): R2=D/2

    if (w0>0):
        norm1 = exp(-4*(x/w0)**2-4*(y/w0)**2)
        norm2 = exp(-4*((x-shear)/w0)**2-4*(y/w0)**2)
    else:
        norm1,norm2=0.*x+1,0.*x+1

    # masking
    bad = (x/R1)**2+(y/R2)**2>1
    if (eps>0): bad[(x/(R1*eps))**2+(y/(R2*eps))**2<=1] = True
    norm1[bad] = 0
    bad = ((x-shear)/R1)**2+(y/R2)**2>1
    if (eps>0): bad[((x-shear)/(R1*eps))**2+(y/(R2*eps))**2<=1] = True
    norm2[bad] = 0

    # convert phase to intensity: 
    s = 0.25*(norm1**2+norm2**2) - 0.5*norm1*norm2*cos(phase)

    if (visual):
        s *= 1 + randn(N,N)/10.
        s -= s.min()
        x = (s/s.max())**(1./2.2)
        x[N//2-1:N//2+1,:]=0
        s = zeros((N,N,3),dtype='float32')
        s[:,:,0] = x

    return s
