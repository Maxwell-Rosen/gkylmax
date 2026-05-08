import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from postgkyl.data import GData, GInterpModal

data = GData("gk_lorentzian_mirror-ion_source_M2_0.gkyl")
grid, vals = GInterpModal(data).interpolate()

jacob = GData("gk_lorentzian_mirror-jacobgeo.gkyl")
J_grid, J_vals = GInterpModal(jacob).interpolate()

psi = np.squeeze(grid[0])
Z = np.squeeze(grid[1])

psi = 0.5 * (psi[:-1] + psi[1:])
Z = 0.5 * (Z[:-1] + Z[1:])

m2 = np.squeeze(vals)

jm2 = np.squeeze(J_vals) * m2

power = np.trapz(np.trapz(jm2, Z, axis=1), psi)
print(f"Integrated M2: {power}")

mass = 1.6726219e-27 * 2.014

power_3d = np.pi * mass * power
print(f"Integrated M2: {power_3d}")
