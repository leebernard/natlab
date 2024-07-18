from numpy import sqrt,cos,sin,pi,arange,newaxis,hstack,array,dot,abs,zeros,histogram2d,exp,array,loadtxt,where,ceil,floor,log
from scipy.fftpack import fft2,fftshift,fft,ifft
from scipy.special import gammaincc,sinc
from scipy.interpolate import RectBivariateSpline,interp1d

from matplotlib.pyplot import plot,xlabel,ylabel,title,imshow

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


class Spot:

    def __init__(self, theta_1=45., theta_2=30., Fn=12., efl1=0.1016, efl2=0.27224, Rt=0.25, Rs=0.095, Rc=1.5, spider=0.05, N=64):
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
        self.x00=0.
        self.y00=0.
        self.z00=0.
        self.x11=0.
        self.lam0=1.5e-6 # reference wavelength
        self.lam=self.lam0
        self.scl1=1.
        self.scl2=1.
        self.M1_shift=[]
        self.M2_shift=[]
        self.verbose=False
        self.theta_1 = theta_1
        self.theta_2 = theta_2
        self.Fn = Fn
        self.efl1 = efl1
        self.efl2 = efl2
        self.Fnd = Fn*efl2/efl1  # F# at detector
        self.Rt = Rt
        self.Rs = Rs
        self.Rc = Rc # cold stop diameter in units of beam diameter
        self.N = N
        self.s0 = 0.
        self.s00 = 0.
        self.sw = 100.e-6 # slit width
        self.sl = 3.e-3 # slit length
        self.slit_x00 = 0. # slit displacement along optical axis from focus
        self.spider = spider  #  width of support arms relative to telescope radius
        self.ps = 18.e-6 # detector pixel size
        self.mask = [] # to hold telescope aperture mask
        self.dmask = [] # to hold instrument aperture mask (cold pupil)

        #
        # fill the telescope aperture, take rays to x00,y00,z00
        if (Rc>1): N1 = 2*int(N*Rc/2)
        else: N1=N
        self.N1 = N1
        x = (arange(1,N)-N//2)*Rt/(N//2)
        x1 = (arange(1,N1)-N1//2)*Rt/(N//2)
        xx = (x[:,newaxis] + 0.*x1[newaxis,:])
        yy = (0*x[:,newaxis] + x1[newaxis,:])
        rad2 = ((x**2)[:,newaxis] + (x1**2)[newaxis,:])
        j0 = rad2<=(Rt*N1/N)**2
        y,z = xx[j0],yy[j0]
        # put chief ray at beginning
        i0 = where((y==0)*(z==0))[0][0]
        y = hstack((y[i0],y[:i0],y[i0+1:]))
        z = hstack((z[i0],z[:i0],z[i0+1:]))
        x = 2*Rt*Fn+0*y

        rad2 = y**2+z**2
        self.j0 = rad2<=Rt**2   #  primary beam
        self.obscured = ~self.j0
        self.obscured[rad2<=Rs**2] = True
        self.obscured[abs(y)<=spider*Rt/2] = True
        self.obscured[abs(z)<=spider*Rt/2] = True
        self.zobscured = abs(y)<=spider*Rt/2
        self.not_obscured = (~self.obscured).sum()/self.j0.sum()

        if (self.verbose): print (f"Initial Beam Size: {z[self.j0].max()-z[self.j0].min():.4f}, {y[self.j0].max()-y[self.j0].min():.4f} ({2*self.Rt:.4f})")

        # create vectors and normalize
        dx,dy,dz=x,-y,-z
        norm = sqrt(dx**2+dy**2+dz**2)
        dx/=norm; dy/=norm; dz/=norm

        self.x0, self.y0, self.z0  = x,y,z
        self.dx0, self.dy0, self.dz0  = dx,dy,dz

        # now run it
        self.raytrace()

    def oap1(self):
        """ reflect off of OAP1 """
        th1 = self.theta_1*pi/180
        f1 = 0.5*self.efl1*(1+cos(th1))

        # rotate to oap coordinate system
        dx,dz,dy = self.dx0*sin(th1) + self.dz0*cos(th1) , self.dz0*sin(th1) - self.dx0*cos(th1),1.*self.dy0
        x0,z0,y0 = self.x00*sin(th1) + self.z00*cos(th1) , self.z00*sin(th1) - self.x00*cos(th1),1.*self.y00

        # path step to oap surface at z = (x^2+y^2)/(4*f1), with z=z0+f1+s*dz
        if (self.theta_1==90): scl_fac=-0.6405511811023621   # (h-zp)/f1
        if (self.theta_1==45): scl_fac=0.008314298505633385
        if (self.theta_1==30): scl_fac=-0.03338612501963257
        scl=self.scl1

        if (len(self.M1_shift)>0):
            x0 -= self.M1_shift[0]
            y0 -= self.M1_shift[1]
            z0 -= self.M1_shift[2]
    
        # z  = ( (x-(1-scl)*xp)^2 + y^2 )/(4*f1*scl) - (1-scl)*f1*scl_fac
        xp = self.efl1*sin(th1)
        s = ( (2*f1*scl*dz-(x0-(1-scl)*xp)*dx-y0*dy) + sqrt( (2*f1*scl*dz-(x0-(1-scl)*xp)*dx-y0*dy)**2+((f1+z0+scl_fac*f1*(1-scl))*4*f1*scl-(x0-(1-scl)*xp)**2-y0**2)*(dx**2+dy**2) ) )/(dx**2+dy**2)

        x,y,z = x0+s*dx,y0+s*dy,z0+f1+s*dz;
        self.s00 = s
        if (self.verbose): print (f"Beam Size on OAP1: {x[self.j0].max()-x[self.j0].min():.4f}, {y[self.j0].max()-y[self.j0].min():.4f} ({self.efl1/self.Fn:.4f})")

        # bounce off the first mirror, normal vector n
        nx,ny,nz = -(x-(1-scl)*xp)/(2*f1*scl),-y/(2*f1*scl),1.
        norm = sqrt(nx*nx+ny*ny+nz*nz) # =sqrt(1+z/f1)
        nx/=norm; ny/=norm; nz/=norm
        vdotn = dx*nx+dy*ny+dz*nz
        # reflection from vector d to r: r = d-2*(d.n)*n
        dx1,dy1,dz1 = dx-2*vdotn*nx, dy-2*vdotn*ny, dz-2*vdotn*nz

        self.x1, self.y1, self.z1  = x,y,z
        self.dx1, self.dy1, self.dz1  = dx1,dy1,dz1

    def move2pupil(self): 
        """ move up to pupil efl1 from surface of OAP1 """
        x,y,z = self.x1,self.y1,self.z1
        dx1,dy1,dz1 = self.dx1,self.dy1,self.dz1

        s = (self.efl1+0.5*self.efl1*(1-cos(self.theta_1*pi/180))-z)/dz1
        self.s00 += s
        x,y,z = x+s*dx1,y+s*dy1,z+s*dz1

        x1,y1,z1=x-x[0],y-y[0],z-z[0]
        # move to pupil that is normal to the chief ray:
        s1=-(x1*dx1[0]+y1*dy1[0]+z1*dz1[0])/(dx1*dx1[0]+dy1*dy1[0]+dz1*dz1[0])
        self.s1 = self.s00 + s1
        self.wfe0 = self.s1[~self.obscured].std()/self.lam0
        if (self.verbose): print ("Wavefront error at pupil (@1.5-micron):",self.wfe0)

        # assume the separation from the pupil to the prism center along ray is efl1
        #  for book-keeping, don't change z
        s = self.efl1/dz1
        self.s00 += s
        x += s*dx1; y += s*dy1

        # reflect back
        dz1 = -dz1
        z += 0.5*self.efl2*(1-cos(self.theta_2*pi/180)) - 0.5*self.efl1*(1-cos(self.theta_1*pi/180))   # z position relative to OAP surfaces

        # center spot on 2nd OAP mirror
        x = self.efl1*sin(self.theta_1*pi/180) - x - self.efl2*sin(self.theta_2*pi/180)

        self.x2, self.y2, self.z2  = x,y,z
        self.dx2, self.dy2, self.dz2  = dx1,dy1,dz1

    def disperse(self):
        """ pass through prism """
        xp = -self.efl2*sin(self.theta_2*pi/180)

        x,y,z = self.x2,self.y2,self.z2
        dx1,dy1,dz1 = self.dx2,self.dy2,self.dz2
        # prism, we assume (approximately true) that the prism center is efl1 away from OAP2
        x,y,z,dx1,dy1,dz1,s1 = prism_raytrace(x,y,z,dx1,dy1,dz1,lam=self.lam,lam0=self.lam0,x0=xp,z0=self.efl1+0.5*self.efl2*(1-cos(self.theta_2*pi/180)))
        self.s0 = self.s00 + s1

        if (self.verbose):
            x1,y1,z1=x-x[0],y-y[0],z-z[0]
            s1=-(x1*dx1[0]+y1*dy1[0]+z1*dz1[0])/(dx1*dx1[0]+dy1*dy1[0]+dz1*dz1[0])
            ff = (self.s0[~self.obscured].mean()*s1[~self.obscured].mean()-(self.s0*s1)[~self.obscured].mean())/s1[~self.obscured].var()  # account for shift to pupil
            print ("Wavefront error after prism:",(ff*s1+self.s0)[1:].std()/self.lam)

        self.x3, self.y3, self.z3  = x,y,z
        self.dx3, self.dy3, self.dz3  = dx1,dy1,dz1

    def oap2(self):
        """ reflect off of OAP2 """
        th2 = self.theta_2*pi/180.
        f2 = 0.5*self.efl2*(1+cos(th2))

        x,y,z = self.x3,self.y3,self.z3
        dx1,dy1,dz1 = self.dx3,self.dy3,self.dz3

        # path length to 2nd mirror
        if (self.theta_2==90): scl_fac=-0.6405511811023621   # (h-zp)/f2
        if (self.theta_2==45): scl_fac=0.008314298505633385
        if (self.theta_2==30): scl_fac=-0.03338612501963257
        #s = (x**2+y**2-4*f2*z)/( 2*f2*dz1-x*dx1-y*dy1 - sqrt( (2*f2*dz1-x*dx1-y*dy1)**2 + (dx1**2+dy1**2)*(4*f2*z-x**2-y**2) ) )
        scl=self.scl2

        xp = -self.efl2*sin(th2)
        s = ( (x-(1-scl)*xp)**2+y**2-4*f2*scl*(z+scl_fac*f2*(1-scl)) )/( 2*f2*scl*dz1-(x-(1-scl)*xp)*dx1-y*dy1 - sqrt( (2*f2*scl*dz1-(x-(1-scl)*xp)*dx1-y*dy1)**2+((z+scl_fac*f2*(1-scl))*4*f2*scl-(x-(1-scl)*xp)**2-y**2)*(dx1**2+dy1**2) ) )

        x += s*dx1; y += s*dy1; z += s*dz1;
        self.s0 += s

        if (self.verbose): print (f"Beam Size on OAP2: {x[self.j0].max()-x[self.j0].min():.4f}, {y[self.j0].max()-y[self.j0].min():.4f} ({self.efl1/self.Fn:.4f})")

        # focus the spot with mirror 2
        nx,ny,nz = -(x-(1-scl)*xp)/(2*f2*scl),-y/(2*f2*scl),1.

        norm = sqrt(nx*nx+ny*ny+nz*nz)
        nx/=norm; ny/=norm; nz/=norm
        vdotn = dx1*nx+dy1*ny+dz1*nz
        dx2 = dx1-2*vdotn*nx
        dy2 = dy1-2*vdotn*ny
        dz2 = dz1-2*vdotn*nz

        if (len(self.M2_shift)>0):
            x += self.M2_shift[0]
            y += self.M2_shift[1]
            z += self.M2_shift[2]

        # rotate from oap coordinates
        x,z = x*sin(th2) + (z-f2)*cos(th2) , (z-f2)*sin(th2) - x*cos(th2)
        dx2,dz2 = dx2*sin(th2) + dz2*cos(th2) , dz2*sin(th2) - dx2*cos(th2)

        self.x4,self.y4,self.z4 = x,y,z
        self.dx4,self.dy4,self.dz4 = dx2,dy2,dz2

    def dichroic(self,through=False,sep=2*25.4e-3):
        """ send to dichroic, go through if through=True """       
        x,y,z = self.x4,self.y4,self.z4
        dx2,dy2,dz2 = self.dx4,self.dy4,self.dz4

        x0 = -self.efl2+sep
        if (through):
            x,y,z,dx2,dy2,dz2,s1 = dichroic_raytrace(x,y,z,dx2,dy2,dz2,lam=self.lam,x0=x0)
            x -= 2195e-6
        else:
            # just move to dichroic front surface z = x0-x
            s1 = (x0-x-z)/(dx2+dz2)
            x,y,z = x+s1*dx2,y+s1*dy2,z+s1*dz2

        self.s0+=s1
        self.x5,self.y5,self.z5 = x,y,z
        self.dx5,self.dy5,self.dz5 = dx2,dy2,dz2


    def focus(self):
        """ move to OAP2 focus """
        x,y,z = self.x5,self.y5,self.z5
        dx2,dy2,dz2 = self.dx5,self.dy5,self.dz5

        s = (self.x11-x)/dx2
        self.s0 += s
        x2,y2,z2 = x+s*dx2, y+s*dy2, z+s*dz2

        # take wavefront to exit pupil
        self.s0 -= (y2-y2[0])*dy2 + (z2-z2[0])*dz2

        # check for shit in diffraction focus relative to geometric focus
        NN = (~self.obscured).sum()
        v1,v2,v3 = self.s0[~self.obscured].mean(), dot(self.s0[~self.obscured],dy2[~self.obscured])/NN, dot(self.s0[~self.obscured],dz2[~self.obscured])/NN
        m11,m12,m22 = 1., dy2[~self.obscured].mean(), dot(dy2[~self.obscured],dy2[~self.obscured])/NN
        m13,m23,m33 = dz2[~self.obscured].mean(), dot(dy2[~self.obscured],dz2[~self.obscured])/NN, dot(dz2[~self.obscured],dz2[~self.obscured])/NN
        mi11,mi12,mi13 = m22*m33-m23**2, m13*m23-m12*m33, m12*m23-m13*m22
        mi22,mi23,mi33 = m11*m33-m13**2, m13*m12-m11*m23, m11*m22-m12**2
        dt = m11*(m22*m33-m23**2) - m12*(m12*m33-m13*m23) + m13*(m12*m23-m13*m22)
        s00     = (mi11*v1+mi12*v2+mi13*v3)/dt
        dy0_fac = (mi12*v1+mi22*v2+mi23*v3)/dt
        dz0_fac = (mi13*v1+mi23*v2+mi33*v3)/dt

        self.s0 -= dy0_fac*dy2 + dz0_fac*dz2
        if (self.verbose): print (f"Diffraction focus shift: 0.0,{-dy0_fac:.6f},{-dz0_fac:.6f}")
        self.wfe = self.s0[~self.obscured].std()/self.lam
        if (self.verbose): print (f"Wavefront Error at focus: {self.wfe:.3f}")

        # correct rays to diffraction focus as opposed to geometric focus
        y2 -= dy0_fac
        z2 -= dz0_fac

        self.dy0_fac, self.dz0_fac = dy0_fac,dz0_fac
        self.x6,self.y6,self.z6 = x2,y2,z2
        self.dx6,self.dy6,self.dz6 = dx2,dy2,dz2

    def raytrace(self,lam=1.5e-6):   
        """ start with telescope rays and move to focus """
        self.lam = lam
        self.mask = []
        self.oap1()
        self.move2pupil()
        self.disperse()
        self.oap2()
        self.dichroic(through=lam>2.5e-6)
        self.focus()

    def changelambda(self,lam=1.5e-6,gothrough=-1):   
        """ start after pupil """
        self.lam = lam
        self.disperse()
        self.oap2()
        if (gothrough==-1): self.dichroic(through=lam>2.5e-6)
        elif (gothrough==0): self.dichroic(through=False)
        else: self.dichroic(through=True)
        self.focus()

    def make_spot(self,f=4,create_image=True):
        """
        f: oversampling
        uses pupil image
        """
        # have N0 samples of the telescope aperture
        # pupil is shape N1,N1 (with N1>N0)
        N0,N1 = self.N,self.N1; N=int(f)*N0
        while (N<N1):
            f *= 2
            N = int(f)*N0

        # (exit) pupil definition
        apert = zeros((N,N),dtype='complex64') # aperture
        a0,b0 = N//2 - N0//2,N//2 - N0//2 + N0
        a,b = N//2 - N1//2,N//2 - N1//2 + N1
        a1,b1 = N1//2 - N0//2,N1//2 - N0//2 + N0

        self.strehl0 = abs(self.mask[a1:b1,a1:b1]/self.ng)**2
        self.spotsum = f**2 / ( pi/4*( 1 - (self.Rs/self.Rt)**2 ) - self.spider*(1-self.Rs/self.Rt) )

        # apply the slit
        if (self.sw>0):
            f1 = 1
            if (f<32): f1 = 32//f   # over-sample to get sinc ringing of aperture
            N2 = N*f1
            ap1 = zeros((N0,N2),dtype='complex64') # aperture
            ap1[:,N2//2-N0//2:N2//2-N0//2+N0] = self.mask[a1:b1,a1:b1]
            phi0 = [0.]
            xx = (arange(N)-N//2)*f/(N//2)
            if (self.slit_x00!=0): phi0 = pi*self.slit_x00*xx**2/(4*self.lam*self.Fn**2)
            if (self.z00!=0): phi0 -= pi*xx*self.z00/(self.lam*self.Fn)
            if (len(phi0)>1): ap1[:,N2//2-N//2:N2//2-N//2+N] *= exp(1j*phi0)
            xx = (arange(N2)-N2//2)*f*f1/(N2//2)
            yy = fftshift( sinc(xx*self.sw/(2*self.Fn*self.lam)) )
            ap1 = ifft( fft(ap1)*fft(yy) )[:,N2//2-N//2:N2//2-N//2+N] * self.sw/(self.Fn*self.N*self.lam)
            if (len(phi0)>1): ap1 *= exp(-1j*phi0)
            #self.slitloss = 1-(abs(ap1)**2).sum()/self.ng
            self.slitloss = 1 - (abs(self.mask[a1:b1,a1:b1]*ap1[:,N//2-N0//2:N//2-N0//2+N0])).sum()/self.ng
            apert[a0:b0] = ap1
        else:
            apert[a0:b0,a0:b0] = self.mask[a1:b1,a1:b1]
            self.slitloss = 1-(abs(apert[a:b])**2).sum()/self.ng

        # apply the cold stop
        apert[a:b,a:b] *= self.dmask*exp(1j*2*pi*self.pupil/self.lam)
        apert[:,:a]=0; apert[:,b:]=0

        self.strehl = abs(apert[a:b,a:b].sum()/self.ng)**2
        self.loss = 1-(abs(apert[a:b,a:b])**2).sum()/self.ng

        # create the focal plane image
        if (create_image):
            im = fftshift(abs(fft2(apert)/self.ng)**2)
            if (self.sl>0):
                xx = (arange(N)-(N-1)/2)*self.lam*self.Fn/f + self.y00
                im[abs(xx)>self.sl/2] = 0

            #im *= (1-self.loss)*self.spotsum/im.sum()
            return im
        else:
            return 0

    def large_mask(self,f=8):
        """
          create a more-finely-sampled pupil mask for the telescope
        """
        N = self.N
        x = (arange(f*N)-(f*N)//2)/((f*N)//2)
        xx = (x[:,newaxis] + 0.*x[newaxis,:])
        yy = (0*x[:,newaxis] + x[newaxis,:])
        rad2 = ((x**2)[:,newaxis] + (x**2)[newaxis,:])

        mask = (rad2<=1)*(rad2>(self.Rs/self.Rt)**2)
        mask[abs(xx)<self.spider/2] = False
        mask[abs(yy)<self.spider/2] = False

        # defocus in units of lam*F#^2
        pupil = mask.astype('complex64')
        #if (self.slit_x00!=0):
        #   defoc = self.slit_x00/(self.lam*self.Fn**2)
        #   pupil *= exp(1j*0.25*pi*defoc*rad2)

        if (self.z00!=0):
            z00 = self.z00/(self.lam*self.Fn)
            pupil *= exp(-1j*pi*x*z00)

        return mask,pupil

    def makepupil(self):
        """
        draw the pupil
        """
        iy = (self.N1//2 + self.y0*(self.N//2)/self.Rt).astype('int16')
        iz = (self.N1//2 + self.z0*(self.N//2)/self.Rt).astype('int16')
        self.pupil = zeros((self.N1,self.N1),dtype='float32')
        self.pupil[iy,iz] = self.s0-self.s0[0]

        if (len(self.mask)==0):
            self.mask = zeros((self.N1,self.N1),dtype='bool')
            self.dmask = zeros((self.N1,self.N1),dtype='bool')
            self.mask[iy[~self.obscured],iz[~self.obscured]] = True
            if (self.Rc>=1): self.dmask[iy,iz] = True
            else:
                rad2 = ( (iy-self.N1//2)**2+(iz-self.N1//2)**2 )/(self.N//2)**2
                h = rad2<=self.Rc**2
                self.dmask[iy[h],iz[h]] = True

            self.dmask[iy[self.zobscured],iz[self.zobscured]] = False
            self.ng = self.mask.sum()

        self.pupil -= self.pupil[self.mask].mean()

    def makeimage(self,oversamp=4,pixellate=False,dzfac=0.):
        """ 
            produce the diffracted spot, sampling is lambda*F#/oversamp
        """
        self.makepupil()
        im = self.make_spot(f=oversamp)

        if (pixellate):
            # image pixel grid:
            pix = self.lam*self.Fnd/oversamp/self.ps
            N = len(im)
            x = pix*(arange(N)-N//2+0.5)
            # new pixel grid
            N1 = 2*(int(x.max()-x.min())//2)
            dy = self.y6[0]/self.ps; dy -= int(dy)
            #dz = (self.z6[0]+self.z00*self.efl2/self.efl1)/self.ps; dz -= int(dz)
            dz = self.z6[0]/self.ps; dz -= int(dz)
            x1 = arange(N1)-(N1-1)/2

            imc = (im/self.spotsum).cumsum(axis=0).cumsum(axis=1)
            deltaz = dzfac*self.lam*self.Fnd/self.ps
            res = RectBivariateSpline(x,x,imc,kx=2,ky=2)
            im2 = res(x1-dy,x1[1:]+0.5*deltaz-dz) - res(x1-dy,x1[:-1]-0.5*deltaz-dz)
            return (im2[1:]-im2[:-1]).clip(0)/(1+deltaz)
        else:
            return im.clip(0)

    def lamgrid(self,oversamp=1,lam1=0.8e-6,lam2=3.5e-6):
        """ define a wavelength grid oversampled with respect to diffraction resolution
            return the grid and associated detector position (wavelength solution)
        """
        start_lam = self.lam
        lam0=lam1
        lam_grid=[lam0]

        while (lam0<lam2):
            self.lam=lam0
            lam0 *= 1 + 1./(oversamp*self.getR())
            lam_grid.append(lam0)

        self.lam = start_lam
        lam_grid = array(lam_grid)
        zsoln = zeros(len(lam_grid),dtype='float32')
        zsoln[1:] = lam_grid[:-1].cumsum()*self.Fnd/oversamp

        # set z to zero at self.lam0 (1.5e-6)
        i0 = (lam_grid<self.lam0).sum()-1
        zsoln -= zsoln[i0] - (lam_grid[i0]-self.lam0)*(zsoln[i0+1]-zsoln[i0])/(lam_grid[i0+1]-lam_grid[i0])

        return lam_grid,zsoln
 

    def getR(self,lam_grid=[]):
        """ get spectral resolution R at lambda assuming lam*F# spatial resolution set by diffraction """
        l0 = self.lam0*1.e6; l20 = l0*l0
        if (len(lam_grid)>0):
            l = lam_grid*1.e6; l2 = l*l
        else:
            l = self.lam*1.e6; l2 = l*l
        # for CaF2 prism:
        a0,a1,a2=0.69913,0.11994,4.35181
        b0,b1,b2=0.093742**2,21.182**2,38.46**2
        n0 = sqrt( 1.33973 + a0*l20/(l20-b0) + a1*l20/(l20-b1) + a2*l20/(l20-b2) )
        n = sqrt( 1.33973 + a0*l2/(l2-b0) + a1*l2/(l2-b1) + a2*l2/(l2-b2) )
        dn_dl = 1./n*abs( a0*b0*l/(l2-b0)**2 + a1*b1*l/(l2-b1)**2 + a2*b2*l/(l2-b2)**2 ) 
        s1p,s0p = 0.5*sqrt(3)*sqrt(n**2-(0.5*n0)**2)-0.25*n0, 0.5*n0/n                                       
        dtheta_dn = 0.5*sqrt(3)/sqrt(1-s0p**2)/sqrt(1-s1p**2) 
        return dtheta_dn*dn_dl*self.efl1/self.Fn*1.e6

    def source_spectrum(self,lamgrid,mag_K=6.,eff=1.):
        """
            eff is efficiency
            return f_source, multiply by dlam/lam to get counts/s
        """
        lam_K = 2.16e-6
        A = eff*pi*self.not_obscured*(self.Rt)**2 # effective area
        f_source = A*(666.7e-26/6.6260701e-34)*10**(-0.4*mag_K)*(lam_K/lamgrid)**2

        #add a line
        #N=len(f_source)
        #f_source[N//2] *= 2

        return f_source

    def background_spectrum(self,lamgrid,Tw=253.,epsw=0.15,epsc=0.1,Tc=120.,lam_max=3.5e-6,epsm=0.95,Tm=120.,Omega=0.1):
        """
            Tw is telescope temperature
            Tc is cold optics temperature
            Tm is wall temperature
            epsw is telescope emissivity
            epsc is emissivity of cold optics
            epsm is emissivity of walls
            Omega solid angle for wall emission
        """
        R = self.getR(lamgrid)
        # warm optics: blackbody from primary, per pixel
        #
        c,hck=2.9979246e+08,0.014387768762454801
        f_tel = epsw * pi/2 * (self.ps/(lamgrid*self.Fnd))**3 * (self.sw/self.ps)*(self.efl2/self.efl1)*(c/lamgrid) / ( exp(hck/lamgrid/Tw) - 1 ) / R

        ########### non-dispersed signal ###########
        # cold optics: blackbody within spectrograph, A.Omega = pi/4*(ps/F#)^2
        #
        f_cold = epsc* 4*c*Tc/hck*(self.ps*Tc/hck/self.Fnd)**2*gammaincc(3,hck/lam_max/Tc)

        # instrument walls
        # solid angle due to collimation
        #
        f_walls = epsm*Omega* 4*c*Tm/hck*(self.ps*Tm/hck)**2*gammaincc(3,hck/lam_max/Tm)

        return f_tel+f_cold+f_walls

    def draw_spectrum(self,channel=1,oversamp=4,oversamp1=2,dzfac=0.,contamination=0.01):
        """
          create an image of the spectrum between lam1 and lam2
          oversamp is oversampling of PSF in units of lambda*F#
          oversamp1 is oversampling of lambda grid relative to spectral FWHM
          returns wavelength grid, source spectrum, background, wavelength grid1, slitloss, all_loss, strehl ratio
          dzfac: blur out the PSF at a given wavelength by dzfac*lambda*F# in the dispersion direction
        """
        l0,ch1,ch2 = loadtxt('ch1ch2.txt',unpack=True)
        if (contamination>0):
            ch1[ch1<contamination] = contamination
            ch2[ch2<contamination] = contamination

        lam1,lam2 = l0.min(),l0.max()

        lamgrid,zsoln=self.lamgrid(oversamp=oversamp1,lam1=0.5e-6,lam2=5.e-6)
        eff1 = interp1d(l0,ch1,bounds_error=False,fill_value=0)(lamgrid)
        eff2 = interp1d(l0,ch2,bounds_error=False,fill_value=0)(lamgrid)

        if (channel==1):
            ch = eff1>0
            eff = eff1[ch]
            through=0
        elif (channel==2):
            ch = eff2>0
            eff = eff2[ch]
            through=1
        else:
            ch = (lamgrid>=0.8e-6)*(lamgrid<=3.5e-6)
            eff = 0.*lamgrid[ch]+1.
            through=0

        lg = lamgrid[ch]
        M = len(lg)
        lam1,lam2 = lg.min(),lg.max()

        flx = self.source_spectrum(lg)
        R = self.getR(lg)
        flx *= eff/(R*oversamp1)

        #y0,z0 = 0., self.z00*self.efl2/self.efl1
        y0,z0 = 0.,0.
        self.changelambda(lam2,gothrough=through)
        yb,zb=(self.y6[0]+y0)/self.ps,(self.z6[0]+z0)/self.ps
        self.changelambda(lam1,gothrough=through)
        ya,za=(self.y6[0]+y0)/self.ps,(self.z6[0]+z0)/self.ps
        Ly,Lz = int(ya)-int(yb),int(za)-int(zb)

        im0 = self.makeimage(oversamp=oversamp,pixellate=True,dzfac=dzfac)*flx[0]
        N0 = len(im0)
        N1,N2 = abs(Ly)+N0,abs(Lz)+N0
        im = zeros((N1,N2),dtype='float32')
        dim = zeros((N1,N2),dtype='float32')
        im[Ly:,Lz:] += im0

        slitloss = zeros(M-1,dtype='float32'); slitloss[0] = self.slitloss
        allloss = zeros(M-1,dtype='float32'); allloss[0] = self.loss
        strehl = zeros(M-1,dtype='float32'); strehl[0] = self.strehl

        for i in range(1,M-1):
            self.changelambda(lg[i],gothrough=through)
            im0 = self.makeimage(oversamp=oversamp,pixellate=True,dzfac=dzfac)
            slitloss[i] = self.slitloss
            allloss[i] = self.loss
            strehl[i] = self.strehl
            N = len(im0)
            a,b = N//2-N0//2,N//2+N0//2+1
            dim[:] = 0
            dim[Ly:,Lz:] = im0[a:b,a:b]*flx[i]
            dy,dz=int((self.y6[0]+y0)/self.ps)-int(ya),int((self.z6[0]+z0)/self.ps)-int(za)
            im[:N1+dy,:N2+dz] += dim[-dy:,-dz:]

        if (channel==1):
            z = arange(N2,dtype='float32')
            z = -(z+int(za)-z[Lz+N0//2])*self.ps
        else:
            self.changelambda(1.5e-6,gothrough=through)
            z = arange(N2,dtype='float32')[::-1]
            z = (z + (self.z6[0]+z0)/self.ps-za - z[Lz+N0//2])*self.ps

        lam1 = interp1d(zsoln,lamgrid,bounds_error=False,fill_value=lamgrid.min())(z)
        bg = self.background_spectrum(lam1) 
        if (channel==1): bg *= interp1d(lamgrid,eff1)(lam1)
        elif (channel==2): bg *= interp1d(lamgrid,eff2)(lam1)

        return lam1,im,bg,lg[:-1],slitloss,allloss,strehl

    def sumOverFwhm(self,x,y,ofac=1.):
        """ sum output of draw spectrum over fwhm along trace """
        l = ofac*x*self.Fnd/self.ps
        l0 = floor(l).astype('int16')
        cy=y.cumsum(axis=0)
        a,b = y.shape
        N = a//2
        ii = arange(b)
        s = cy[N+l0,ii]*(1-l+l0) + cy[N+l0+1,ii]*(l-l0)
        s -= cy[N-l0,ii]*(1-l+l0) + cy[N-l0-1,ii]*(l-l0)
        mask = zeros((a,b),dtype='bool')
        for i in range(b):
            mask[N-l0[i]:N+l0[i]+1,i] = True

        return s,mask

    def spotplot(self):
        """
        """
        y = (self.y6[~self.obscured]-self.y6[0])/self.ps
        z = (self.z6[~self.obscured]-self.z6[0])/self.ps
        airy = 1.22*self.lam*self.Fnd/self.ps 
        phi = 2*pi*arange(100)/99
        plot (y,z,'o',alpha=0.1)
        plot (airy*cos(phi),airy*sin(phi),'k:',alpha=0.5)
        xlabel("X Pixel",fontsize=14); ylabel("Y Pixel",fontsize=14)
        strl = exp(-(2*pi*self.wfe)**2)
        title(r"""RMS Wavefront Error: %.3f/$\lambda$ (Strehl$\approx$%.2f)""" % (self.wfe,strl))

    def implot(self,image,scl=1.e-6,logscale=True):
        """   
        """   
        if (logscale): imshow(-log(scl+image),cmap='gray')
        else: imshow(-image,cmap='gray')
        title("""Strehl: %.2f, Slit-loss: %.1f%%, All-loss: %.1f%%""" % (min(self.strehl,1),self.slitloss*100,self.loss*100))
