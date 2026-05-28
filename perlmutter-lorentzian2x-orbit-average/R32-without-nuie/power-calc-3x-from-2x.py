import numpy as np

integjm2 = 2.171348e+28 # pgkyl gk_lorentzian_mirror-ion_source_M2_0.gkyl gk_lorentzian_mirror-jacobgeo.gkyl interp ev "f[0] f[1] *" integ 1,0 info
integjm0 = 8.982309e+15 # pgkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl gk_lorentzian_mirror-jacobgeo.gkyl interp ev "f[0] f[1] *" integ 1,0 info

mass = 1.6726219e-27 * 2.014

flux_3d = 2 * np.pi * integjm0
power_3d = np.pi * mass * integjm2
print(f"Integrated flux: {flux_3d}") # 5.6416312466120536e+16
print(f"Integrated power: {power_3d}") # 229.79287075807366

#(1.187386e31 - 6.6125e30) / 0.001 = 5.26186e+33
