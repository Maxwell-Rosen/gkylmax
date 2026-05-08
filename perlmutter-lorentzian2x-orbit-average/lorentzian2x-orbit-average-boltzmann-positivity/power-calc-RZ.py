import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from postgkyl.data import GData, GInterpModal

data = GData("gk_lorentzian_mirror-ion_source_M2_0.gkyl", mapc2p_name="gk_lorentzian_mirror-mapc2p_deflated.gkyl")
grid, vals = GInterpModal(data).interpolate()

cells = np.shape(grid)
grid_out = []
for i in range((cells[0])):
  grid_out.append(0.5 * (grid[i][:-1, :-1] + grid[i][1:, 1:]))

R = np.squeeze(grid_out[0])
Z = np.squeeze(grid_out[1])
f = np.squeeze(vals)

# Set values of f to zero outside the range Z = -0.98:0.98
f[(Z < -0.98) | (Z > 0.98)] = 0.0

# Integrate along the R direction first
f_integrated = np.trapz(f, R, axis=0)
# Integrate along Z second
Z_1d  = Z[0,:]
f_integrated = np.trapz(f_integrated, Z_1d)

mass = 1.6726219e-27 * 2.014

power_1d = 7.682495e+32 * 1/2 * mass
power_2d = f_integrated * 1/2 * mass
print(f"Integrated M2: {power_2d}")
print(f"1D Power: {power_1d}")
print(f"Ratio: {power_1d / power_2d}")

# print(f_integrated)
# plt.plot(R, f_integrated)
# plt.xlabel('R')
# plt.ylabel('Integrated M2')
# plt.title('Integrated M2 along Z vs R')
# plt.grid()
# plt.show()

# plt.pcolormesh(Z, R, f)
# plt.show()

# print(np.shape(R))
# print(np.shape(Z))
# print(np.shape(f))