import numpy as np
import matplotlib.pyplot as plt


def vel_map(vpc):
  b = 1.45
  vel_range = 16
  vp = np.zeros_like(vpc)
  linear_velocity_threshold = 1./3.
  vpar_threshold = 1/b*np.arctan(linear_velocity_threshold * np.tan(b))
  frac_linear = vpar_threshold
  print(f'vpar_threshold = {vpar_threshold}')
  func_frac = np.tan(frac_linear*b) / np.tan(b)
  mask = abs(vpc) < frac_linear
  vp[mask]  = vel_range * vpc[mask] * func_frac / frac_linear
  vp[~mask] = vel_range * np.tan(vpc[~mask]*b) / np.tan(b)
  return vp
  # vel_range = 30
  # return vel_range * np.tan(vpc*b) / np.tan(b) * s + s * vpc

vpc = np.linspace(-1, 1, 32)
vp = vel_map(vpc)
vp_half = vp[len(vp)//2:]

fig, axs = plt.subplots(1, 2, figsize=(6, 4))


# Plot the scatter plot
axs[0].plot(vp, 'b.')
axs[0].set_xlabel('Cell number')
axs[0].set_ylabel('$v_{||}$ ($v_{th}$ units)')
axs[0].tick_params(axis='both', which='major')

# Plot the histogram
axs[1].hist(vp, bins=60)
axs[1].set_xlabel('$v_{||}$ ($v_{th}$ units)')
axs[1].set_ylabel('Frequency')
axs[1].tick_params(axis='both', which='major')


# Put a title on the whole figure
fig.suptitle('Non-uniform velocity grid for $v_{||}$')
fig.tight_layout()
# Set all font sizes to 32
plt.savefig('velocity-grid.png', dpi=1000)
plt.show()