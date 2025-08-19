import postgkyl as pg
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as multip
from concurrent.futures import ProcessPoolExecutor
import time
from matplotlib.colors import LogNorm

matplotlib.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.size': 12,
    'axes.titlesize': 18,
    'axes.labelsize': 22,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10
})

last_frame = 100 # MUST UPDATE THIS TO 100 LATER

# Universal constants (matching input_file.c)
eps0 = 8.8541878128e-12  # F/m (permittivity of free space)
mu0 = 4 * np.pi * 1e-7  # H/m (permeability of free space)
eV = 1.602176634e-19  # J (elementary charge / electron volt)
mp = 1.67262192369e-27  # kg (proton mass)
me = 9.1093837015e-31  # kg (electron mass)

# Plasma parameters (from input_file.c)
mi = 2.014 * mp  # kg (deuteron mass)
qi = eV  # C (ion charge)
qe = -eV  # C (electron charge)
Te0 = 940 * eV  # J (electron temperature)
n0 = 3e19  # m^-3 (number density)
B_p = 0.53  # T (magnetic field strength)
beta = 0.4  # plasma beta
tau = (B_p**2) * beta / (2.0 * mu0 * n0 * Te0) - 1.0  # temperature ratio
Ti0 = tau * Te0  # J (ion temperature)
kperpRhos = 0.1  # normalized perpendicular wavenumber

# Derived parameters
vti = np.sqrt(Ti0 / mi)  # m/s (ion thermal velocity)
vte = np.sqrt(Te0 / me)  # m/s (electron thermal velocity)
c_s = np.sqrt(Te0 / mi)  # m/s (sound speed)
omega_ci = eV * B_p / mi  # rad/s (ion cyclotron frequency)
rho_s = c_s / omega_ci  # m (sound gyroradius)
kperp = kperpRhos / rho_s  # m^-1 (perpendicular wavenumber)

# Collision frequencies
nuFrac = 1.0
elc_nuFrac = 1/5.489216862238348
logLambdaElc = 6.6 - 0.5 * np.log(n0 / 1e20) + 1.5 * np.log(Te0 / eV)
nuElc = elc_nuFrac * nuFrac * logLambdaElc * (eV**4) * n0 / (6 * np.sqrt(2) * (np.pi**(3/2)) * (eps0**2) * np.sqrt(me) * (Te0**(3/2)))
logLambdaIon = 6.6 - 0.5 * np.log(n0 / 1e20) + 1.5 * np.log(Ti0 / eV)
nuIon = nuFrac * logLambdaIon * (eV**4) * n0 / (12 * (np.pi**(3/2)) * (eps0**2) * np.sqrt(mi) * (Ti0**(3/2)))

# Geometry parameters
z_min = -2.0
z_max = 2.0
psi_min = 1e-6
psi_eval = 1e-3
psi_max = 3e-3

# Velocity space parameters
vpar_max_elc = 30 * vte
mu_max_elc = me * (3 * vte)**2 / (2 * B_p)
vpar_max_ion = 30 * vti
mu_max_ion = mi * (3 * vti)**2 / (2 * B_p)

# Source parameters
ion_source_amplitude = 1e20  # m^-3/s
ion_source_sigma = 0.5
ion_source_temp = 5000 * eV  # J


# Plot 0th BiMaxwellian Moments
data_BiMax0 = pg.GData('gk_wham-ion_BiMaxwellianMoments_0.gkyl')
interp_BiMax0 = pg.GInterpModal(data_BiMax0, 1, 'ms')
x, dens0 = interp_BiMax0.interpolate(0)
x, upar0 = interp_BiMax0.interpolate(1)
x, Tpar_div_m0 = interp_BiMax0.interpolate(2)
x, Tperp_div_m0 = interp_BiMax0.interpolate(3)

dens0 = np.squeeze(dens0)
upar0 = np.squeeze(upar0)/c_s
Tpar_div_m0 = np.squeeze(Tpar_div_m0)
Tperp_div_m0 = np.squeeze(Tperp_div_m0)

Tpar0 = Tpar_div_m0 * mi / eV  # Convert to eV
Tperp0 = Tperp_div_m0 * mi / eV  # Convert to eV

avg_density = np.mean(np.abs(dens0))
avg_upar = np.mean(np.abs(upar0))
avg_Tpar = np.mean(np.abs(Tpar0))
avg_Tperp = np.mean(np.abs(Tperp0))

rms_density = np.zeros(last_frame)
rms_upar = np.zeros(last_frame)
rms_Tpar = np.zeros(last_frame)
rms_Tperp = np.zeros(last_frame)

for i in range(last_frame):
  data_BiMax = pg.GData('gk_wham-ion_BiMaxwellianMoments_'+str(i)+'.gkyl')
  interp_BiMax = pg.GInterpModal(data_BiMax, 1, 'ms')
  x, dens = interp_BiMax.interpolate(0)
  x, upar = interp_BiMax.interpolate(1)
  x, Tpar_div_m = interp_BiMax.interpolate(2)
  x, Tperp_div_m = interp_BiMax.interpolate(3)

  dens = np.squeeze(dens)
  upar = np.squeeze(upar)/c_s
  Tpar_div_m = np.squeeze(Tpar_div_m)
  Tperp_div_m = np.squeeze(Tperp_div_m)

  Tpar = Tpar_div_m * mi / eV  # Convert to eV
  Tperp = Tperp_div_m * mi / eV  # Convert to

  rms_density[i] = np.std(dens - dens0) / avg_density
  rms_upar[i] = np.std(upar - upar0) / avg_upar
  rms_Tpar[i] = np.std(Tpar - Tpar0) / avg_Tpar
  rms_Tperp[i] = np.std(Tperp - Tperp0) / avg_Tperp

# Calculate RMS error for the field

data_field0 = pg.GData('gk_wham-field_0.gkyl')
interp_field0 = pg.GInterpModal(data_field0, 1, 'ms')
x, phi0 = interp_field0.interpolate(0)
phi0 = np.squeeze(phi0) * eV / Te0
avg_phi0 = np.mean(np.abs(phi0))
rms_phi = np.zeros(last_frame)
for i in range(last_frame):
  data_field = pg.GData('gk_wham-field_'+str(i)+'.gkyl')
  interp_field = pg.GInterpModal(data_field, 1, 'ms')
  x, phi = interp_field.interpolate(0)
  phi = np.squeeze(phi) * eV / Te0
  rms_phi[i] = np.std(phi - phi0) / avg_phi0

fig, ax = plt.subplots( 2, 3, figsize=(10,6))

xaxis = np.arange(1, last_frame)*1e-6 * vti / 2

ax[0, 0].plot(xaxis, rms_density[1:], label='σ Ion Density', color="#0080D5")  # Use a strong orange from ColorBrewer/ggplot
ax[0, 0].set_ylabel(r'$CV(n_i)$')

ax[0, 1].plot(xaxis, rms_upar[1:], label='σ Ion Parallel Velocity', color="#0080D5") 
ax[0, 1].set_ylabel(r'$CV(u_{||}/c_s)$')

ax[1, 0].plot(xaxis, rms_Tpar[1:], label='σ Ion Parallel Temperature', color="#0080D5")
ax[1, 0].set_ylabel(r'$CV(T_{||})$')

ax[1, 1].plot(xaxis, rms_Tperp[1:], label='σ Ion Perpendicular Temperature', color="#0080D5")
ax[1, 1].set_ylabel(r'$CV(T_{\perp})$')

ax[1, 2].plot(xaxis, rms_phi[1:], label='σ Field', color="#0080D5")
ax[1, 2].set_ylabel(r'$CV(e \phi / T_e)$')

# Remove ax[0,2] since it is not used
fig.delaxes(ax[0, 2])

# Add a grid to each subplot
for a in ax.flat:
  a.grid(True)
  a.set_xlim(0, xaxis[-1])  # Set x-axis limits for all subplots

# Set the x-axis label for the bottom row
ax[1, 0].set_xlabel(r'$t / \tau_{\rm transit}$')
ax[1, 1].set_xlabel(r'$t / \tau_{\rm transit}$')
ax[1, 2].set_xlabel(r'$t / \tau_{\rm transit}$')

plt.tight_layout()
# plt.show()
plt.savefig('bimaxwellian_moments_and_field_std.pdf')
plt.close()  # Close figure to free memory

