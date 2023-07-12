
# pull the size of the fine wavelength grid
fine_size = fine_wl.shape[0]
# need to figure out whether I pull right side or left side
# see searchsorted docstring

# find the nearest wavelengths corresponding to the wavelength grid boundaries
ii = np.searchsorted(fine_wl, wlgrid)

# make a cumlative sum of the fine wavelengths
Fpc = np.zeros(fine_size + 1, dtype='float64')
Fpc[1:] = Fp.cumsum()

# calculate the mean of values within the given wavelength bins
Fint = (Fpc[ii[1:] + 1] - Fpc[ii[:-1]]) / (ii[1:] - ii[:-1] + 1)

# calculate the mean wavelength of each bin
mean_wavelengths = (wlgrid[:-1] + wlgrid[1:]) / 2



# test which method is correct


fpc_total_signal = np.sum( Fpc[ii[1:] + 1] - Fpc[ii[:-1]])
fp_total_signal = np.sum( Fp[:ii[-1]])

print('With +1:', fpc_total_signal - fp_total_signal)

fpc_total_signal = np.sum( Fpc[ii[1:]] - Fpc[ii[:-1]])
fp_total_signal = np.sum( Fp[:ii[-1]])


print('Without +1:', fpc_total_signal - fp_total_signal)
# got a better result without the +1





import numpy as np
import matplotlib.pyplot as plt

from numpy.random import rand
from numpy.random import randn

'''
Nat's code
'''
#
# simulate some data for testing
#   fine data "dat" are on grid lam
#   coarse data "dat0" to be on grid lam0
#   notice that we always have len(lam) = len(dat)+1 since we are working with bin edges
#
go=True
lam = rand(300); lam.sort(); lam[0]=0; lam[-1]=1.
dat = randn(299) + 20
while (go):
    lam = rand(300); lam.sort(); lam[0]=0; lam[-1]=1.
    dat = randn(299) + 20
    lam0 = rand(30); lam0.sort(); lam0[0]=0; lam0[-1]=1
    ii = lam.searchsorted(lam0)
    # make sure there are no duplicate indices
    if (len(ii)==len(np.unique(ii))): go=False


#
# now the algorithm 2 ways, ("dat0") with cumsum and edge mop up
#    and linear interpolation ("dat1") which doesn't need edge mop up
#  (however get the normalization to go from sum to average will be trickier)
#


cdat = np.zeros(len(lam) + 1,dtype='float64')
cdat[1:] = dat.cumsum()
from scipy.interpolate import interp1d
res = interp1d(lam,cdat)
dat1 = res(lam0[1:]) - res(lam0[:-1])

dat0 = cdat[ii[1:]] - cdat[ii[:-1]]
norm = 1. * (ii[1:] - ii[:-1])
# mop up the bin edges by shuffling some signal between bins:
delta1 = ( lam[ii[1:]]-lam0[1:] )/( lam[ii[1:]]-lam[ii[1:]-1] )
dat0 -= delta1*dat[ii[1:]-1]
dat0[1:] += delta1[:-1]*dat[ii[1:-1]-1]
norm -= delta1
norm[1:] += delta1[:-1]

if (ii[0]>0):
    delta1 = ( lam[ii[0]]-lam0[0] )/( lam[ii[0]]-lam[ii[0]-1] )
    dat0[0] += delta1*dat[ii[0]-1]
    norm[0] += delta1

# note: the average is now dat0/norm

# to validate the binning between the two algorithms, check that dat0.sum() == dat1.sum()
print (dat1.sum()/dat0.sum()-1)

# you can also make a plot to check the first algorithm
plt.plot (0.5*(lam[1:]+lam[:-1]),dat)
plt.plot (0.5*(lam0[1:]+lam0[:-1]),dat0/norm)


# check that the signal is preserved
fpc_total_signal = np.sum(dat0)
fp_total_signal = np.sum( dat[:ii[-1]])


print('Differnece in signal summation:', fp_total_signal - fpc_total_signal)

plt.plot(dat0)
plt.plot(dat[:ii[-1]])


output = dat0 / norm

# compare to the function version
from toolkit import improved_non_uniform_tophat

test_dat, _ = improved_non_uniform_tophat(lam0, lam, dat)




