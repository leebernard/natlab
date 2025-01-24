from scipy.fftpack import fft2,fftshift
from numpy import zeros,arange,exp,newaxis,pi

def make_spot(defoc=2.,tilt=0.,N0=128,f=8,ap_ratio=0.38,spider=0.05,phi=[],coma=False):
    """
    f: oversampling
       spot pixels size is lambda*F#/f, or  ~about 1 FWHM/f
    N0: defines number of aperture samples: N0 samples across aperture
       size of final image spans N0 FWHM's of the PSF
    ap_ratio: central obscuration ratio
    spider: thickness of spiders / radius
    coma: treat coma instead of defocus
    defoc: units of lambda*F#^2
    tilt: units of lambda*F#
    normalize: =False (central pixel -> 1) = True (sum -> 1)
    """
    if (len(phi)>0): N0 = len(phi)
    N=int(f)*N0

    if (N0<=defoc): print (f"Warning: large defocus, need N0>{defoc:.0f}")

    # (exit) pupil definition
    apert = zeros((N,N),dtype='complex128') # aperture
    a,b = N//2 - N0//2,N//2 - N0//2 + N0

    if (len(phi)>0):
        good = phi>0
        apert[a:b,a:b] = good*exp(1j*2*pi*phi)
    else:
        #x = (arange(N0)-(N0-1)/2)/(N0/2)
        x = (arange(N0)-N0//2)/(N0/2)
        rad2 = x[:,newaxis]**2+x[newaxis,:]**2
        good = (rad2>ap_ratio**2)*(rad2<=1.)
        if (spider>0):
            good[abs(x)<spider/2] = False
            good[:,abs(x)<spider/2] = False

        if (coma):
            apert[a:b,a:b] = good * exp(1j*pi*defoc*rad2*x[:,newaxis]/4) # coma
        else:
            apert[a:b,a:b] = good * exp(1j*pi*(tilt*x[:,newaxis]+defoc*rad2/4))  # tilt and/or defocus

    ng = good.sum()
    return fftshift(abs(fft2(apert)/ng)**2)
