import numpy as np
import matplotlib.pyplot as plt

test_data_file = 'excite/DATA_1_warmup_20230531.csv'
warmup_data = np.loadtxt(test_data_file)

fname = 'excite/ccc_warmupdata_run1/ccc_ls_ktemp_d4'
temperature_data = np.fromfile(fname, dtype='float32')

fname = 'excite/ccc_warmupdata_run1/ccc_time'
time_data = np.fromfile(fname, dtype='float32')

# cvs file index 11 corresponds to KST index 368285
cvs_sample_start = 11
kst_start = 368285
kst_end = 1237856  # This is where we lost telemetry
time_slice = time_data[kst_start:kst_end]
time_slice = time_slice - time_slice[0]
minutes = time_slice/60

plt.plot(warmup_data[11:2951])
plt.ylabel('mV')

cvs_sample_start = 11
kst_start = 368285



