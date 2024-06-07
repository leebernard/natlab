import numpy as np
import matplotlib.pyplot as plt

test_data_file = 'excite/DATA_1_warmup_20230531.csv'
voltage_data = np.loadtxt(test_data_file)  # units are mV

fname = 'excite/ccc_warmupdata_run1/ccc_ls_ktemp_d4'
temperature_data = np.fromfile(fname, dtype='float32')

fname = 'excite/ccc_warmupdata_run1/ccc_time'
time_data = np.fromfile(fname, dtype='float32')

# cvs file index 11 corresponds to KST index 368285

kst_start = 362850
kst_end = 1237879  # Presumably this is where we lost telemetry
time_slice = time_data[kst_start:kst_end]
time_slice = time_slice - time_slice[0]
kst_minutes = time_slice / 60  # assuming time data is in seconds

plate_temp = temperature_data[kst_start:kst_end]

cvs_sample_start = 11
cvs_sample_end = int(kst_minutes[-1])
voltage_slice = voltage_data[cvs_sample_start:cvs_sample_end]  # convert to mV

voltage_time = np.arange(voltage_slice.size)
# np.searchsorted()  # finds nearest neighbor

# map the plate temperature data to the voltage data, using time to interpolate
plate_temp_interp = np.interp(voltage_time, kst_minutes, plate_temp)

# check if the mapping was done correctly
# plt.scatter(voltage_time, plate_temp_interp)
# plt.ylabel('Kelvin')
# plt.xlabel('Seconds')

fig0, ax0 = plt.subplots(figsize=(8, 6), tight_layout=True)
ax0.plot(plate_temp_interp, voltage_slice, label='Photodiode voltage curve')
ax0.set_xlabel('Spectrometer plate temperature (K)')
ax0.set_ylabel('mV')
ax0.legend()

fig1, ax1 = plt.subplots(figsize=(8, 6), tight_layout=True)
color1 = 'tab:blue'
ax1.plot(voltage_time, voltage_slice, label='Photodiode signal warmup curve', color=color1)
ax1.set_xlabel('Minutes')
ax1.set_ylabel('mV', color=color1)
# ax1.legend()

ax2 = ax1.twinx()
color2='tab:orange'
ax2.plot(voltage_time, plate_temp_interp,  label='Temperature warmup curve', color=color2)
ax2.set_xlabel('Minutes')
ax2.set_ylabel('Kelvin', color=color2)
# ax2.legend(loc=0)

fig1.legend(loc=2, bbox_to_anchor=[.15, .95])  # best doesn't work




