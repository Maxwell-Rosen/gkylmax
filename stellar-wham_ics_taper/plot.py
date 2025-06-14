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

species = 0
psi = 16
lineLength = 1
u = 5
theta_ix = 5

f_eval = f_dist[species, psi, lineLength, u, theta_ix]
psi_eval = psiGrid[psi]
z_eval = zGrid[lineLength]
v_eval = uGrid[u]
theta_eval = theta[theta_ix]

# print(np.shape(phi))
# print(np.shape(psiGrid))
# print(np.shape(zGrid))
# print(psiGrid)
# print(zGrid)
# plt.pcolormesh(zGrid, psiGrid, phi[1:, 1:])
# # plt.pcolormesh(phi)
# plt.colorbar()
# plt.title("Phi")
# plt.set_cmap("plasma")
# plt.ylabel("psi")
# plt.xlabel("z")
# plt.savefig("python-plots/phi.png")
# plt.show()

# zGrid = zGrid * 1e-2
# phi = phi * 1e3
# psiGrid = psiGrid * -1e-8
# print(psiGrid)
# plt.plot(zGrid, phi[-1, :])
# plt.title("Phi at psi = " + str(psiGrid[-1]))
# plt.ylabel("phi, V")
# plt.xlabel("z, m")
# plt.savefig("python-plots/phi_98.png")
# plt.show()

# plt.plot(zGrid, phi[0, :])
# plt.title("Phi at psi = " + str(psiGrid[0]))
# plt.ylabel("phi, V")
# plt.xlabel("z, m")
# plt.savefig("python-plots/phi_0.png")
# plt.show()

# Integrate to get density
print("making the pl")
fig, axs = plt.subplots(1, 2)

cax0 = axs[0].pcolormesh(f_dist[1, 0, 0, :, :])
fig.colorbar(cax0, ax=axs[0])
axs[0].set_title('Subfigure 1')
axs[0].set_xlabel('x')
axs[0].set_ylabel('y')

cax1 = axs[1].pcolormesh(f_dist[1, 0, 1, :, :])
fig.colorbar(cax1, ax=axs[1])
axs[1].set_title('Subfigure 2')
axs[1].set_xlabel('x')
axs[1].set_ylabel('y')

plt.show()

# f_slice_0 = f_dist[0, psi, 0, :, :]
# plt.pcolormesh(theta, uGrid*v_norm, f_slice_0[1:, 1:])
# plt.ylim(0, 0.01*v_norm)
# plt.colorbar()
# plt.title("Ion distribution function at lineLength = 0")
# plt.set_cmap("plasma")
# plt.ylabel("u")
# plt.xlabel("theta")
# plt.savefig("python-plots/distf_ion_0.png")
# plt.close()

# f_slice_1 = f_dist[0, psi, -1, :, :]
# plt.pcolormesh(theta, uGrid*v_norm, f_slice_1[1:, 1:])
# plt.ylim(0, 0.01*v_norm)
# plt.colorbar()
# plt.title("Ion distribution function at lineLength = " + str(zGrid[-1]))
# plt.set_cmap("plasma")
# plt.ylabel("u")
# plt.xlabel("theta")
# plt.savefig("python-plots/distf_ion_98.png")
# plt.close()

# f_slice_2 = f_dist[0, psi, 80, :, :]
# plt.pcolormesh(theta, uGrid*v_norm, f_slice_2[1:, 1:])
# plt.ylim(0, 0.01*v_norm)
# plt.colorbar()
# plt.title("Ion distribution function at lineLength = " + str(zGrid[80]))
# plt.set_cmap("plasma")
# plt.ylabel("u")
# plt.xlabel("theta")
# plt.savefig("python-plots/distf_ion_60.png")
# plt.close()

# f_slice_3 = f_dist[1, psi, 0, :, :]
# plt.pcolormesh(theta, uGrid*v_norm, f_slice_3[1:, 1:])
# plt.ylim(0, v_norm)
# plt.colorbar()
# plt.title("Elc distribution function at lineLength = 0")
# plt.set_cmap("plasma")
# plt.ylabel("u")
# plt.xlabel("theta")
# plt.savefig("python-plots/distf_elc_0.png")
# plt.close()

# f_slice_4 = f_dist[1, psi, -1, :, :]
# plt.pcolormesh(theta, uGrid*v_norm, f_slice_4[1:, 1:])
# plt.ylim(0, v_norm)
# plt.colorbar()
# plt.title("Elc distribution function at lineLength = " + str(zGrid[-1]))
# plt.set_cmap("plasma")
# plt.ylabel("u")
# plt.xlabel("theta")
# plt.savefig("python-plots/distf_elc_98.png")
# plt.close()

# f_slice_5 = f_dist[1, psi, 80, :, :]
# plt.pcolormesh(theta, uGrid*v_norm, f_slice_5[1:, 1:])
# plt.ylim(0, v_norm)
# plt.colorbar()
# plt.title("Elc distribution function at lineLength = " + str(zGrid[80]))
# plt.set_cmap("plasma")
# plt.xlabel("theta")
# plt.ylabel("u")
# plt.savefig("python-plots/distf_elc_60.png")
# plt.close()

# # Convert from polar to cartesian
# v = uGrid*v_norm
# theta = theta
# vperp = np.outer(v,np.sin(theta))**2
# vpar = np.outer(v,np.cos(theta))

# f_slice_6 = f_dist[1, psi, -1, :, :]
# plt.pcolormesh(vpar, vperp, f_slice_6[1:, 1:])
# plt.ylim(0, 5e19)
# plt.colorbar()
# plt.title("Elc distribution function at lineLength = " + str(zGrid[-1]))
# plt.set_cmap("plasma")
# plt.ylabel("v_perp**2")
# plt.xlabel("v_par")
# plt.savefig("python-plots/distf_elc_98_cartesian.png")
# plt.close()

# f_slice_7 = f_dist[1, psi, 80, :, :]
# plt.pcolormesh(vpar, vperp, f_slice_7[1:, 1:])
# plt.ylim(0, 5e19)
# plt.colorbar()
# plt.title("Elc distribution function at lineLength = " + str(zGrid[80]))
# plt.set_cmap("plasma")
# plt.ylabel("v_perp**2")
# plt.xlabel("v_par")
# plt.savefig("python-plots/distf_elc_60_cartesian.png")
# plt.close()

# f_slice_8 = f_dist[1, psi, 0, :, :]
# plt.pcolormesh(vpar, vperp, f_slice_8[1:, 1:])
# plt.ylim(0, 5e19)
# plt.colorbar()
# plt.title("Elc distribution function at lineLength = 0")
# plt.set_cmap("plasma")
# plt.ylabel("v_perp**2")
# plt.xlabel("v_par")
# plt.savefig("python-plots/distf_elc_0_cartesian.png")
# plt.close()

# f_slice_9 = f_dist[0, psi, -1, :, :]
# plt.pcolormesh(vpar, vperp, f_slice_9[1:, 1:])
# plt.ylim(0, 5e16)
# plt.xlim(-1e8, 1e8)
# plt.colorbar()
# plt.title("Ion distribution function at lineLength = " + str(zGrid[-1]))
# plt.set_cmap("plasma")
# plt.ylabel("v_perp**2")
# plt.xlabel("v_par")
# plt.savefig("python-plots/distf_ion_98_cartesian.png")
# plt.close()

# f_slice_10 = f_dist[0, psi, 80, :, :]
# plt.pcolormesh(vpar, vperp, f_slice_10[1:, 1:])
# plt.ylim(0, 5e16)
# plt.xlim(-1e8, 1e8)
# plt.colorbar()
# plt.title("Ion distribution function at lineLength = " + str(zGrid[80]))
# plt.set_cmap("plasma")
# plt.ylabel("v_perp**2")
# plt.xlabel("v_par")
# plt.savefig("python-plots/distf_ion_60_cartesian.png")
# plt.close()

# f_slice_11 = f_dist[0, psi, 0, :, :]
# plt.pcolormesh(vpar, vperp, f_slice_11[1:, 1:])
# plt.ylim(0, 5e16)
# plt.xlim(-1e8, 1e8)
# plt.colorbar()
# plt.title("Ion distribution function at lineLength = 0")
# plt.set_cmap("plasma")
# plt.ylabel("v_perp**2")
# plt.xlabel("v_par")
# plt.savefig("python-plots/distf_ion_0_cartesian.png")
# plt.close()
