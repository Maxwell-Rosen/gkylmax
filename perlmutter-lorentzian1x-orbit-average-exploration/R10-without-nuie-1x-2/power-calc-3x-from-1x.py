import numpy as np

integjm2 = 8.258368e+30   # pgkyl gk_lorentzian_mirror-ion_source_M2_0.gkyl gk_lorentzian_mirror-jacobgeo.gkyl interp ev "f[0] f[1] *" integ 0 info
integjm0 = 3.415077e+18 # pgkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl gk_lorentzian_mirror-jacobgeo.gkyl interp ev "f[0] f[1] *" integ 0 info

psi_max = 0.00262927
mass = 1.6726219e-27 * 2.014

power_3d = np.pi * mass * integjm2 * psi_max
flux_3d = 2 * np.pi * integjm0 * psi_max
print(f"Integrated flux: {flux_3d}") # 5.636485845501093e+16
print(f"Integrated power: {power_3d}") # 229.82984265577193
