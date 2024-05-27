import numpy as np
import matplotlib.pyplot as plt

# with open("zGrid.bin", "rb") as binary_file:
#     zGrid = binary_file.read()
#     numbers = struct.unpack('f' * (len(zGrid) // 8), zGrid)
#     print(numbers)
B0 = np.fromfile("B0.bin")
BdB0 = np.fromfile("BdB0.bin")
charge = np.fromfile("charge.bin")
f_dist_elc = np.fromfile("f_dist_elc.bin")
f_dist_ion = np.fromfile("f_dist_ion.bin")
mass = np.fromfile("mass.bin")
phi = np.fromfile("phi.bin")
psiGrid = np.fromfile("psiGrid.bin")
theta = np.fromfile("theta.bin")
uGrid = np.fromfile("uGrid.bin")
vGrid = np.fromfile("vGrid.bin")
v_norm = np.fromfile("v_norm.bin")
zGrid = np.fromfile("zGrid.bin")

f_dist_elc = f_dist_elc.reshape((len(psiGrid), len(zGrid), len(uGrid), len(theta)))
f_dist_ion = f_dist_ion.reshape((len(psiGrid), len(zGrid), len(uGrid), len(theta)))

# plt.pcolormesh(f_dist_elc[:,:,0,0])
# plt.colorbar()
# plt.show()


fig, axs = plt.subplots(1, 2)

cax0 = axs[0].pcolormesh(f_dist_elc[0, 0, :, :])
fig.colorbar(cax0, ax=axs[0])
axs[0].set_title('Subfigure 1')
axs[0].set_xlabel('x')
axs[0].set_ylabel('y')

cax1 = axs[1].pcolormesh(f_dist_elc[0, 1, :, :])
fig.colorbar(cax1, ax=axs[1])
axs[1].set_title('Subfigure 2')
axs[1].set_xlabel('x')
axs[1].set_ylabel('y')

plt.show()

# f is read f[species, psi, z, v, theta]
density = np.zeros(len(zGrid))
psi = 0
for i in range(len(zGrid)):
  f_slice_0 = f_dist_elc[psi, i, :, :]
  integrate_theta = np.trapz(f_slice_0, theta, axis=1)
  v_integral_theta = np.multiply(integrate_theta, uGrid*v_norm)
  density[i] = np.trapz(v_integral_theta, uGrid*v_norm) * 2 * np.pi

# density[0] = density[1]

plt.plot(zGrid, density*1e8)
plt.ylabel("Density, cm^-3")
plt.xlabel("Z, cm")
plt.title("Density vs Z")
plt.show()