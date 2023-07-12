import numpy as np
path = '/home/lee/natlab/excite_optic_data/'
file = 'short_pass_filter_sp0.txt'
sp0_data = np.loadtxt(path+file, skiprows=1)
print('shape', sp0_data.shape)

wl = sp0_data[:,0]/1000  # convert to um
trans = sp0_data[:,1]/100  # convert to fractional
refl = sp0_data[:,2]/100  # convert to fractional

np.savetxt(path+'scratch_wl.txt', wl)
np.savetxt(path+'scratch_trans.txt', trans)
np.savetxt(path+'scratch_refl.txt', refl)

