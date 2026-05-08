import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from postgkyl.data import GData, GInterpModal

data = GData("gk_lorentzian_mirror-ion_source_M2_0.gkyl", mapc2p_name="gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl")
grid, vals = GInterpModal(data).interpolate()

cells = np.shape(grid)
grid_out = []
for i in range((cells[0])):
  grid_out.append(0.5 * (grid[i][:-1, :-1] + grid[i][1:, 1:]))

psi = np.squeeze(grid_out[0])
z = np.squeeze(grid_out[1])
f = np.squeeze(vals)

# Set values of f to zero outside the range z = -0.98:0.98
f[(z < -0.98) | (z > 0.98)] = 0.0

# Integrate along the z direction (axis 1)
f_integrated = np.trapz(f, z, axis=1)

plt.plot(psi, f_integrated)
plt.xlabel('psi')
plt.ylabel('Integrated M2')
plt.title('Integrated M2 along z vs psi')
plt.grid()
plt.show()

# plt.pcolormesh(z, psi, f)
# plt.show()

print(np.shape(psi))
print(np.shape(z))
print(np.shape(f))