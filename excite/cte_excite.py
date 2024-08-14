import matplotlib
matplotlib.use("qt5agg")
import matplotlib.pyplot as plt
import numpy as np

# Stefan-Boltzmann constant
sigma = 6.570373e-8  # W/(m^2 K^4)

def ti_6_4_cte(T):
    ### titanium 6al 4V
    # returns units of um/(M K)
    a = -1.711e2
    b = -2.140e-1
    c = 4.807e-3
    d = -7.111e-6
    return a + b*T + c*T**2 + d*T**3


def al6061_cte(T):
    ### titanium 6al 4V
    a = -4.1277e2
    b = -3.0389E-1
    c = 8.7696E-3
    d = -9.9821E-6
    return a + b*T + c*T**2 + d*T**3


def al5083_cte(T):
    ### titanium 6al 4V
    a = -4.1277E2
    b = -3.0389E-1
    c = 8.7696E-3
    d = -9.9821E-6
    return a + b*T + c*T**2 + d*T**3



T = np.linspace(300, 100, num=500)

fig, ax = plt.subplots()
ax.plot(T, ti_6_4_cte(T), label='Ti 6-4')
ax.plot(T, al6061_cte(T), label='Al 6061')
ax.plot(T, al5083_cte(T), label='Al 5083')

ax.legend()

### emission

# cooling power

# cooling area is approx a 11 inch flat circle.
# the circle has an edge, but it also has sections removed.
# any areas on the circle that doen't have coating are compensated for by the coating on the mounts
# the mounts overall have more emitting surface area than their footprints
# so the flat circle approx most likely underestimates cooling power

bench_r = 11 * .0254 / 2  # bench radius in meters

emission_area = np.pi * bench_r**2

# stephan-boltzmann law
# power = sigma * T^4 * area (emissivity ~ 1)
emission_power = sigma * T**4 * emission_area

T_steady = 120  # kelvin
emission_steady = sigma * T_steady**4 * emission_area

shell_power_5k = sigma * (T - 5)**4 * emission_area
shell_power_10k = sigma * (T - 10)**4 * emission_area
shell_power_3k = sigma * (T - 3)**4 * emission_area


cooling_power_5k = emission_power - shell_power_5k
cooling_power_10k = emission_power - shell_power_10k
cooling_power_3k = emission_power - shell_power_3k

T_shell = 295.8
T_plate = 296.
plate_emission = sigma * T_plate**4 * emission_area
shell_emission = sigma * T_shell**4 * emission_area
radiative_cooling = plate_emission - shell_emission
print(f'Radiative cooling: {radiative_cooling: .2f} W')

fig, (ax1, ax2) = plt.subplots(2, figsize=(6,8),  tight_layout=True)
ax1.plot(T, emission_power, label='Emissive power')
ax1.plot(T_steady, emission_steady, label='Operating Temperature', linestyle='None', marker='o', color='tab:red')
ax1.annotate(f'emissive power: {emission_steady:.2f} W', xy=(T_steady, emission_steady), xytext=(165, 5),
            arrowprops=dict(facecolor='black'))
ax1.set_xlabel('Bench Temperature (K)')
ax1.set_ylabel('Emissive power of the spectrograph coating (W)')
ax1.xaxis.set_inverted(True)

ax1.legend()

ax2.plot(T, cooling_power_5k, label='delta = 5 K')
ax2.plot(T, cooling_power_10k, label='delta = 10 K')
ax2.plot(T, cooling_power_3k, label='delta = 3 K')
ax2.set_xlabel('Bench Temperature (K)')
ax2.set_ylabel('Radiative cooling power  (W)')
ax2.xaxis.set_inverted(True)
ax2.legend()




