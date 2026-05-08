import numpy as np

integjm2 = 1.117152e+29 # pgkyl gk_lorentzian_mirror-ion_source_M2_0.gkyl gk_lorentzian_mirror-jacobgeo.gkyl interp ev "f[0] f[1] *" integ 1,0 info

mass = 1.6726219e-27 * 2.014

power_3d = np.pi * mass * integjm2
print(f"Integrated M2: {power_3d}")
