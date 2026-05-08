import numpy as np

integjm2 = 7.527231e+32 # pgkyl gk_lorentzian_mirror-ion_source_M2_0.gkyl gk_lorentzian_mirror-jacobgeo.gkyl interp ev "f[0] f[1] *" integ 0 info

psi_max = 0.00265309
mass = 1.6726219e-27 * 2.014

power_3d = np.pi * mass * integjm2 * psi_max
print(f"Integrated M2: {power_3d}")
