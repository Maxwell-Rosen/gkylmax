import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg

output_dir = "64-cells/"

filename_ion_distf = "64-cells/gk_wham-ion_300.gkyl"
filename_mc2nu_pos = "64-cells/gk_wham-mc2nu_pos.gkyl"
filename_mc2nu_vel = "64-cells/gk_wham-ion_mapc2p_vel.gkyl"

poly_order = 1
ions = pg.GData(filename_ion_distf)
ion_interp = pg.GInterpModal(ions, poly_order, 'gkhyb')
x, f = ion_interp.interpolate()
f = f[:,:,:,0]

mc2nu_pos = pg.GData(filename_mc2nu_pos)
mc2nu_pos_interp = pg.GInterpModal(mc2nu_pos, poly_order, 'ms')
x, mc2nu = mc2nu_pos_interp.interpolate()
z = mc2nu[:,0]

mc2nu_vel = pg.GData(filename_mc2nu_vel)
num_interp_points = 3 # Need to do this because of the gkhybrid basis, so p2 in vpar
mc2nu_vel_interp = pg.GInterpModal(mc2nu_vel, poly_order, 'ms', num_interp_points)
x, vel_grid = mc2nu_vel_interp.interpolate()
vpar = vel_grid[:,0,0]

mc2nu_vel_interp = pg.GInterpModal(mc2nu_vel, poly_order, 'ms')
x, vel_grid = mc2nu_vel_interp.interpolate()
mu = vel_grid[0,:,0]

# Write these arrays as binary files
output_file = f"{output_dir}/f_dist_ion.bin"
with open(output_file, 'wb') as binary_file:
    binary_file.write(f.tobytes())

output_file = f"{output_dir}/zGrid.bin"
with open(output_file, 'wb') as binary_file:
    binary_file.write(z.tobytes())

output_file = f"{output_dir}/vparGrid.bin"
with open(output_file, 'wb') as binary_file:
    binary_file.write(vpar.tobytes())

output_file = f"{output_dir}/muGrid.bin"
with open(output_file, 'wb') as binary_file:
    binary_file.write(mu.tobytes())