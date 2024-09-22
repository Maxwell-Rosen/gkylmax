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
phi = phi.reshape((len(psiGrid), len(zGrid)))
BdB0 = BdB0.reshape((len(psiGrid), len(zGrid)))
# plt.pcolormesh(f_dist_elc[:,:,0,0])
# plt.colorbar()
# plt.show()


# fig, axs = plt.subplots(1, 2)

# cax0 = axs[0].pcolormesh(f_dist_elc[0, 0, :, :])
# fig.colorbar(cax0, ax=axs[0])
# axs[0].set_title('Subfigure 1')
# axs[0].set_xlabel('x')
# axs[0].set_ylabel('y')

# cax1 = axs[1].pcolormesh(f_dist_elc[0, 1, :, :])
# fig.colorbar(cax1, ax=axs[1])
# axs[1].set_title('Subfigure 2')
# axs[1].set_xlabel('x')
# axs[1].set_ylabel('y')

# plt.show()

# f is read f[species, psi, z, v, theta]
density = np.zeros(len(zGrid))
psi = 0
for i in range(len(zGrid)):
  f_slice_0 = f_dist_elc[psi, i, :, :]
  integrate_theta = np.trapz(f_slice_0, theta, axis=1)
  v_integral_theta = np.multiply(integrate_theta, uGrid*v_norm)
  density[i] = np.trapz(v_integral_theta, uGrid*v_norm) * 2 * np.pi

# density[0] = density[1]

# plt.plot(zGrid, density*1e8)
# plt.ylabel("Density, cm^-3")
# plt.xlabel("Z, cm")
# plt.title("Density vs Z")
# plt.show()

# plot phi
psiGrid_interp = np.zeros(len(psiGrid)+1)
psiGrid_interp[0] = psiGrid[0]
psiGrid_interp[-1] = psiGrid[-1]
for i in range(1, len(psiGrid)):
  psiGrid_interp[i] = (psiGrid[i] + psiGrid[i-1]) / 2

zGrid_interp = np.zeros(len(zGrid)+1)
zGrid_interp[0] = zGrid[0]
zGrid_interp[-1] = zGrid[-1]
for i in range(1, len(zGrid)):
  zGrid_interp[i] = (zGrid[i] + zGrid[i-1]) / 2
  

# plt.pcolormesh(phi) 
# plt.colorbar()
# plt.show()

# plt.plot(phi[:,-1],'r',label='index -1')
# plt.plot(phi[:,-2],'b',label='index -2')
# plt.plot(phi[:,-3],'g',label='index -3')
# plt.plot(phi[:,-4],'y',label='index -4')
# phi_shape = phi.shape
# plt.plot(psiGrid,phi[:, phi_shape[1]//2],'r',label='index 0')
# plt.xscale('log')
# plt.legend()
# plt.show()


plt.pcolormesh(BdB0)
plt.colorbar()
plt.show()