from numpy import sqrt,cos,sin,pi,arange,newaxis,hstack,array,dot,abs

def dichroic_raytrace(x,y,z,dx,dy,dz,lam=1.0e-6,T=5.e-3,x0=0.,z0=0.):
    """
      fused-silica dichroic
      T is thickness of dichroic
      x,y,z,dx,dy,dz starting ray positions and directions (altered and returned)
    """

    # for fused silica
    a0,a1,a2 = 0.6961663,0.4079426,0.8974794
    b0,b1,b2 = 0.0684043**2,0.1162414**2,9.896161**2
    l2 = (lam*1.e6)**2
    n = sqrt( 1. + a0*l2/(l2-b0) + a1*l2/(l2-b1) + a2*l2/(l2-b2) )

    ca,sa = 1./sqrt(2),1./sqrt(2)   # first surface normal

    # move to first face, ca*(x+t*dx-x0) + sa*(z+t*dz-z0) = 0
    t = -(ca*(x-x0)+sa*(z-z0))/(ca*dx+sa*dz); s0 = 1.*t
    x,y,z = x+t*dx,y+t*dy,z+t*dz

    # now refract 
    mu,nx,nz = 1./n,ca,sa
    n_i = nx*dx+nz*dz
    fac = sqrt(abs(1-mu**2*(1-n_i**2)))-mu*n_i
    dx,dy,dz = mu*dx + nx*fac, mu*dy, mu*dz + nz*fac

    # move to second face, ca*(x+t*dx-x0-T*sqrt(2)) + sa*(z+t*dz-z0) = 0
    t = -(ca*(x-x0-T*sqrt(2))+sa*(z-z0))/(ca*dx+sa*dz); s0 += n*t
    x,y,z = x+t*dx,y+t*dy,z+t*dz

    # now refract 
    mu,nx,nz = n,ca,sa
    n_i = nx*dx+nz*dz
    fac = sqrt(abs(1-mu**2*(1-n_i**2)))-mu*n_i
    dx,dy,dz = mu*dx + nx*fac, mu*dy, mu*dz + nz*fac

    return x,y,z,dx,dy,dz,s0


def prism_raytrace(x,y,z,dx,dy,dz,lam=1.0e-6,b=25.4e-3,lam0=1.5e-6,x0=0.,z0=0.):
    """
      equilateral ca-f2 prism, located at x,y,z=x0,0,z0
      b is length of prism face
      evaluate for wavelength lam; lam0 is symmetric passage wavelength
      x,y,z,dx,dy,dz starting ray positions and directions (altered and returned)
    """
    # for CaF2 prism:
    a0,a1,a2=0.69913,0.11994,4.35181
    b0,b1,b2=0.093742**2,21.182**2,38.46**2
    l2 = (lam0*1.e6)**2
    n0 = sqrt( 1.33973 + a0*l2/(l2-b0) + a1*l2/(l2-b1) + a2*l2/(l2-b2) )
    l2 = (lam*1.e6)**2
    n = sqrt( 1.33973 + a0*l2/(l2-b0) + a1*l2/(l2-b1) + a2*l2/(l2-b2) )

    ca,sa = 0.5,0.5*sqrt(3) # cos(alpha),sin(alpha), alpha=60, first surface normal
    stheta0,ctheta0 = n0/2, sqrt(1-(n0/2)**2) # sin(theta0),cos(theta0), incidence theta0 = arcsin(n0/2)

    # initial on-axis rays move in dx,dy,dz = 0,0,-1
    # rotate by half of deflection angle (2*theta0-alpha) to get to prism frame
    crot,srot = 0.5*(sqrt(3)*ctheta0+stheta0), 0.5*(sqrt(3)*stheta0-ctheta0)

    dx,dz = crot*dx-srot*dz, crot*dz+srot*dx
    x,z = x0 + crot*(x-x0) - srot*(z-z0), z0 + crot*(z-z0) + srot*(x-x0)

    # move to first face, ca*(x+t*dx-x0) + sa*(z+t*dz-z0-0.25*b) = 0
    t = -(ca*(x-x0) + sa*(z-z0-0.25*b))/(ca*dx+sa*dz); s0 = 1.*t
    x,y,z = x+t*dx,y+t*dy,z+t*dz

    # now refract 
    mu,nx,nz = 1./n,-ca,-sa
    n_i = nx*dx+nz*dz
    fac = sqrt(abs(1-mu**2*(1-n_i**2)))-mu*n_i
    dx,dy,dz = mu*dx + nx*fac, mu*dy, mu*dz + nz*fac

    # move to second face, ca*(x+t*dx-x0) - sa*(z+t*dz-z0+0.25*b) = 0
    t = -(ca*(x-x0) - sa*(z-z0+0.25*b))/(ca*dx-sa*dz); s0 += n*t
    x,y,z = x+t*dx,y+t*dy,z+t*dz

    # now refract 
    mu,nx,nz = n,ca,-sa
    n_i = nx*dx+nz*dz
    fac = sqrt(abs(1-mu**2*(1-n_i**2)))-mu*n_i
    dx,dy,dz = mu*dx + nx*fac, mu*dy, mu*dz + nz*fac

    # rotate again by half of deflection angle
    dx,dz = crot*dx-srot*dz, crot*dz+srot*dx
    x,z = x0 + crot*(x-x0) - srot*(z-z0), z0 + crot*(z-z0) + srot*(x-x0)

    return x,y,z,dx,dy,dz,s0


def spot(x00=0., y00=0., z00=0., x11=0., theta_1=45., theta_2=30., Fn=12., efl1=0.1016, efl2=0.27224, lam=1.5e-6, N=203,
         step=8, scl1=1., scl2=1., M1_shift=[], M2_shift=[], diffract_focus=True, through_dichroic=False, verbose=True):
    """
      sending focused light through an 2 OAP relay
        (default units are degrees and inches)

      x00,y00,z00: starting focal plane positions
      x11: detector position
      theta_1 : deflection angle of first OAP (collimator) in degrees
      theta_2 : deflection angle of second OAP (camera) in degrees
      Fn : telescope F-number
      efl1,efl2: (effective) focal lengths of two mirrors
      lam: wavelength of rays (default 1.5e-6)
      N : number of rays to trace (overrides to len(theta) if vector)
      scl1,scl2: allow for differential shrinkage of OAPs relative to spectrograph (default 1,1)
      M1_shift: (list) offset M1 by [a,b,c]
      M2_shift: (list) offset M2 by [a,b,c]
      diffract_focus: return rays positions at diffraction rather than geometric focus (default True)
          Note: an input focused beam goes to the diffraction focus, while a straight laser goes to the geometric focus
              (these have opposite signs in z when z00!=0)
              (changes in M1 shift can do similar things)

      step==8 (default), returns: y,z,wfe
          y,z are positions (same units as efl1,efl2) of the traced rays on the detector
          wfe is wavefront error estimate in units of lam (lambda)
      step<8: returns x,y,z,dx,dy,dz,s0 (ray positions and unit vector and OPD s0)

      Note: Different values of step return earlier ray positions.
        step=-1 initial rays (y-z plane, x-axis toward slit at x=0)
          -> rotation spectrograph to OAPs
        step=0 rays in OAP coordinate system (z-axis in collimation direction)
        step=1 rays on OAP1 surface
        step=2 rays at pupil between OAP1 and OAP2
        step=3 rays before prism
        step=4 rays after prism
        step=5 rays on OAP2
          -> rotation OAPS to spectrograph
        step=6 final rays on OAP2 (x-axis toward detector at x=0)
        step=7 rays after dichroic (if through_dichroic=True, otherwise same as step=6)
        step=8 rays at detector (y-z plane)

      Note: first ray (e.g., x[0],y[0],z[0]) is chief ray
          This ray is blocked, but it defines the beam center.
        So, y[0],z[0] at end is center of spot.  To plot the spot centered at zero do:
        plot ((y[1:]-y[0]),(z[1:]-z[0]),'o')
    """

    #
    # fill the telescope aperture, take rays to x00,y00,z00
    Rs,Rt=0.095,0.25  # radius of secondary, radius of primary
    x = (arange(N)-N//2)*Rt/(N//2)  # generate sample locations
    xx = x[:,newaxis] + 0*x[newaxis,:]
    yy = 0*x[:,newaxis] + x[newaxis,:]
    rad2 = xx**2+yy**2
    j = (rad2>Rs**2)*(rad2<=Rt**2)  # filter out locations not in the telescope aperture
    z = xx[j]; y=yy[j]; x=2*Rt*Fn+0*y  # set locations for position -1 (telescope aperture)
    xx=0;yy=0;rad2=0;j=0  # reset variables
    xs,ys,zs = 1.*x,1.*y,1.*z

    if (verbose): print (f"Initial Beam Size: {z.max()-z.min():.4f}, {y.max()-y.min():.4f} ({2*Rt:.4f})")

    dx,dy,dz=x,-y,-z
    # add in chief ray
    dx=hstack((2*Rt*Fn,dx))
    dy=hstack((0,dy))
    dz=hstack((0,dz))
    norm = sqrt(dx**2+dy**2+dz**2)
    dx/=norm; dy/=norm; dz/=norm

    if (step==-1): return x,y,z,dx,dy,dz,0*dz

    # miror focal lengths parameters
    th1,th2 = theta_1*pi/180, theta_2*pi/180.  # convert to radians
    f1 = 0.5*efl1*(1+cos(th1))  # define focal length of parent parabola
    f2 = 0.5*efl2*(1+cos(th2))

    # to send in the laser directly:
    #z00 = hstack((0,zs))*0.63e-3/(2*Rt) + z00
    #y00 = hstack((0,ys))*0.63e-3/(2*Rt) + y00
    #dx = 0*dx + 1.; dy*=0; dz*=0

    # rotate to oap coordinate system
    dx,dz = dx*sin(th1) + dz*cos(th1) , dz*sin(th1) - dx*cos(th1)
    x0,z0,y0 = x00*sin(th1) + z00*cos(th1) , z00*sin(th1) - x00*cos(th1),1.*y00

    if (step==0): return x0,y0,z0,dx,dy,dz,0*dz

    # path step to oap surface at z = (x^2+y^2)/(4*f1), with z=z0+f1+s*dz
    #s = ( 2*f1*dz-x0*dx-y0*dy + sqrt( (2*f1*dz-x0*dx-y0*dy)**2 + (dx**2+dy**2)*(4*f1*(f1+z0)-x0**2-y0**2) ) )/(dx**2+dy**2)
    if (theta_1==90): scl_fac=-0.6405511811023621   # (h-zp)/f1
    if (theta_1==45): scl_fac=0.008314298505633385
    if (theta_1==30): scl_fac=-0.03338612501963257
    scl=scl1

    if (len(M1_shift)>0):
        x0 -= M1_shift[0]
        y0 -= M1_shift[1]
        z0 -= M1_shift[2]
    
    # old s = ( 2*f1*dz/scl-x0*dx-y0*dy + sqrt( (2*f1*dz/scl-x0*dx-y0*dy)**2 + (dx**2+dy**2)*(4*f1/scl*(f1+z0-scl_fac*f1*(scl-1))-x0**2-y0**2) ) )/(dx**2+dy**2)
    # z  = ( (x-(1-scl)*xp)^2 + y^2 )/(4*f1*scl) - (1-scl)*f1*scl_fac
    xp = efl1*sin(th1)
    s = ( (2*f1*scl*dz-(x0-(1-scl)*xp)*dx-y0*dy) + sqrt( (2*f1*scl*dz-(x0-(1-scl)*xp)*dx-y0*dy)**2+((f1+z0+scl_fac*f1*(1-scl))*4*f1*scl-(x0-(1-scl)*xp)**2-y0**2)*(dx**2+dy**2) ) )/(dx**2+dy**2)

    x,y,z = x0+s*dx,y0+s*dy,z0+f1+s*dz; s0 = 1.*s
    if (step==1): return x,y,z,dx,dy,dz,s0

    if (verbose): print (f"Beam Size on OAP1: {x.max()-x.min():.4f}, {y.max()-y.min():.4f} ({efl1/Fn:.4f})")

    # bounce off the first mirror, normal vector n
    nx,ny,nz = -(x-(1-scl)*xp)/(2*f1*scl),-y/(2*f1*scl),1.
    norm = sqrt(nx*nx+ny*ny+nz*nz) # =sqrt(1+z/f1)
    nx/=norm; ny/=norm; nz/=norm
    vdotn = dx*nx+dy*ny+dz*nz
    # reflection from vector d to r: r = d-2*(d.n)*n
    dx1,dy1,dz1 = dx-2*vdotn*nx, dy-2*vdotn*ny, dz-2*vdotn*nz

    # move up to pupil efl1 from surface of OAP1
    s = (efl1+0.5*efl1*(1-cos(th1))-z)/dz1; s0 += s
    x,y,z = x+s*dx1,y+s*dy1,z+s*dz1

    if (step==2):
        x1,y1,z1=x-x[0],y-y[0],z-z[0]
        # move to pupil that is normal to the chief ray:
        s1=-(x1*dx1[0]+y1*dy1[0]+z1*dz1[0])/(dx1*dx1[0]+dy1*dy1[0]+dz1*dz1[0])
        if (verbose): print ("Wavefront error:",(s1+s0)[1:].std()/lam)
        return x,y,z,dx1,dy1,dz1,s0

    # assume the separation from the pupil to the prism center along ray is efl1
    #  for book-keeping, don't change z
    s = efl1/dz1; s0 += s
    x += s*dx1; y += s*dy1

    # reflect back
    dz1 = -dz1
    z += 0.5*efl2*(1-cos(th2)) - 0.5*efl1*(1-cos(th1))   # z position relative to OAP surfaces

    # center spot on 2nd OAP mirror
    x = efl1*sin(th1) - x - efl2*sin(th2)
    xp = -efl2*sin(th2)

    if (step==3): return x,y,z,dx1,dy1,dz1,s0

    # prism, we assume (approximately true) that the prism center is efl1 away from OAP2
    x,y,z,dx1,dy1,dz1,s1 = prism_raytrace(x,y,z,dx1,dy1,dz1,lam=lam,x0=xp,z0=efl1+0.5*efl2*(1-cos(th2)))
    s0 += s1

    if (step==4):
        x1,y1,z1=x-x[0],y-y[0],z-z[0]
        s1=-(x1*dx1[0]+y1*dy1[0]+z1*dz1[0])/(dx1*dx1[0]+dy1*dy1[0]+dz1*dz1[0])
        ff = (s0.mean()*s1.mean()-(s0*s1).mean())/s1.var()  # account for shift to pupil
        if (verbose): print ("Wavefront error:",(ff*s1+s0)[1:].std()/lam)
        return x,y,z,dx1,dy1,dz1,s0

    # path length to 2nd mirror
    if (theta_2==90): scl_fac=-0.6405511811023621   # (h-zp)/f2
    if (theta_2==45): scl_fac=0.008314298505633385
    if (theta_2==30): scl_fac=-0.03338612501963257
    #s = (x**2+y**2-4*f2*z)/( 2*f2*dz1-x*dx1-y*dy1 - sqrt( (2*f2*dz1-x*dx1-y*dy1)**2 + (dx1**2+dy1**2)*(4*f2*z-x**2-y**2) ) )
    scl=scl2

    # old: s = (x**2+y**2-4*f2/scl*(z-scl_fac*f2*(scl-1)))/( 2*f2*dz1/scl-x*dx1-y*dy1 - sqrt( (2*f2*dz1/scl-x*dx1-y*dy1)**2 + (dx1**2+dy1**2)*(4*f2/scl*(z-scl_fac*f2*(scl-1))-x**2-y**2) ) )
    xp = -efl2*sin(th2)
    s = ( (x-(1-scl)*xp)**2+y**2-4*f2*scl*(z+scl_fac*f2*(1-scl)) )/( 2*f2*scl*dz1-(x-(1-scl)*xp)*dx1-y*dy1 - sqrt( (2*f2*scl*dz1-(x-(1-scl)*xp)*dx1-y*dy1)**2+((z+scl_fac*f2*(1-scl))*4*f2*scl-(x-(1-scl)*xp)**2-y**2)*(dx1**2+dy1**2) ) )

    x += s*dx1; y += s*dy1; z += s*dz1;
    s0 += s

    if (step==5): return x,y,z,dx1,dy1,dz1,s0

    if (verbose): print (f"Beam Size on OAP2: {x.max()-x.min():.4f}, {y.max()-y.min():.4f} ({efl1/Fn:.4f})")

    # focus the spot with mirror 2
    #nx,ny,nz = -x/(2*f2),-y/(2*f2),1.
    nx,ny,nz = -(x-(1-scl)*xp)/(2*f2*scl),-y/(2*f2*scl),1.

    norm = sqrt(nx*nx+ny*ny+nz*nz)
    nx/=norm; ny/=norm; nz/=norm
    vdotn = dx1*nx+dy1*ny+dz1*nz
    dx2 = dx1-2*vdotn*nx
    dy2 = dy1-2*vdotn*ny
    dz2 = dz1-2*vdotn*nz

    if (len(M2_shift)>0):
        x += M2_shift[0]
        y += M2_shift[1]
        z += M2_shift[2]

    # rotate from oap coordinates
    x,z = x*sin(th2) + (z-f2)*cos(th2) , (z-f2)*sin(th2) - x*cos(th2)
    dx2,dz2 = dx2*sin(th2) + dz2*cos(th2) , dz2*sin(th2) - dx2*cos(th2)

    if (step==6): return x,y,z,dx2,dy2,dz2,s0

    if (through_dichroic):
        x,y,z,dx2,dy2,dz2,s1 = dichroic_raytrace(x,y,z,dx2,dy2,dz2,lam=lam,x0=-efl2+2*25.4e-3)
        s0+=s1
        x -= 2195e-6

    if (step==7): return x,y,z,dx2,dy2,dz2,s0

    # move to OAP2 focus
    s = (x11-x)/dx2; s0 += s
    x2,y2,z2 = x+s*dx2, y+s*dy2, z+s*dz2

    # take wavefront to exit pupil
    s0 -= (y2-y2[0])*dy2 + (z2-z2[0])*dz2

    # check for shit in diffraction focus relative to geometric focus
    NN = len(s0[1:])
    v1,v2,v3 = s0[1:].mean(), dot(s0[1:],dy2[1:])/NN, dot(s0[1:],dz2[1:])/NN
    m11,m12,m22 = 1., dy2[1:].mean(), dot(dy2[1:],dy2[1:])/NN
    m13,m23,m33 = dz2[1:].mean(), dot(dy2[1:],dz2[1:])/NN, dot(dz2[1:],dz2[1:])/NN
    mi11,mi12,mi13 = m22*m33-m23**2, m13*m23-m12*m33, m12*m23-m13*m22
    mi22,mi23,mi33 = m11*m33-m13**2, m13*m12-m11*m23, m11*m22-m12**2
    dt = m11*(m22*m33-m23**2) - m12*(m12*m33-m13*m23) + m13*(m12*m23-m13*m22)
    s00     = (mi11*v1+mi12*v2+mi13*v3)/dt
    dy0_fac = (mi12*v1+mi22*v2+mi23*v3)/dt
    dz0_fac = (mi13*v1+mi23*v2+mi33*v3)/dt

    s0 -= dy0_fac*dy2 + dz0_fac*dz2
    if (verbose): print (f"Diffraction focus shift: 0.0,{-dy0_fac:.6f},{-dz0_fac:.6f}")
    wfe = s0.std()/lam
    if (verbose): print (f"Wavefront Error: {wfe:.3f}")

    # correct rays to diffraction focus as opposed to geometric focus
    if (diffract_focus):
        y2 -= dy0_fac
        z2 -= dz0_fac

    if (step==8): return y2,z2,wfe
