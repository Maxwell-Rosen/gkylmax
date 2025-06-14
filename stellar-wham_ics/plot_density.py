import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# Open the HDF5 file
hf = h5py.File('cql3d_f_RZ.h5', 'r')
# Read the datasets
v_norm = hf['v_norm'][()] # Velocity normalization coefficient. cm/s
psiGrid = hf['psiGrid'][:]
zGrid = hf['zGrid'][:]
uGrid = hf['uGrid'][:] # Normalized spherical velocity
vGrid = hf['vGrid'][:]
theta = hf['theta'][:] # Pitch angle. Radians
charge = hf['charge'][:]
mass = hf['mass'][:]
BdB0 = hf['BdB0'][:]
B0 = hf['B0'][:]
phi = hf['phi'][:]
f_dist = hf['f_dist'][:]

print(BdB0.shape)
tile_B0 = np.transpose(np.tile(B0, (len(zGrid), 1)))
BGrid = np.multiply(BdB0, tile_B0)

species = 0
psi = 16
print(v_norm)
print(uGrid)

# f is read f[species, psi, z, v, theta]
density = np.zeros((len(psiGrid), len(zGrid)))
for j in range(len(psiGrid)):
  for i in range(len(zGrid)):
    f_slice_0 = f_dist[1, j, i, :, :]
    integrate_theta = np.trapz(f_slice_0, theta, axis=1)
    v_integral_theta = np.multiply(integrate_theta, uGrid*v_norm)
    density[j,i] = np.trapz(v_integral_theta, uGrid*v_norm) * 2 * np.pi

# density[0] = density[1]
density[:,0] = density[:,1]
psiGrid *= -1e-8
zGrid *= 1e-2

# plt.pcolormesh(zGrid, psiGrid, density)
# plt.title("Electron Density vs Z and Psi")
# plt.ylabel("Psi")
# plt.xlabel("Z, cm")
# plt.colorbar()
# plt.show()

plt.pcolormesh(zGrid, psiGrid, phi)
plt.title("Phi vs Z and Psi")
plt.ylabel("Psi")
plt.xlabel("Z, cm")
plt.colorbar()
plt.show()

# plt.plot(zGrid, density[psi,:]*1e8)
# plt.ylabel("Density, cm^-3")
# plt.xlabel("Z, cm")
# plt.title("Density vs Z")
# plt.show()