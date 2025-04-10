#[ ........................................................... ]#
#[
#[ Maxwell Rosen
#[ plot.py
#[ 2021-04-20
#  For plotting comparison of uniform and nonuniform grid runs.
#[
#[ ........................................................... ]#

import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from matplotlib.colors import LogNorm
import multiprocessing
from scipy.integrate import cumulative_trapezoid as cumtrapz
import imageio.v2 as imageio
from scipy.optimize import curve_fit

frame = 50

uniformBiMax_name = './data-unif/gk_wham-ion_BiMaxwellianMoments_'+str(frame)+'.gkyl'
nonunifBiMax_name = './data-nonunif/gk_wham-ion_BiMaxwellianMoments_'+str(frame)+'.gkyl'

uniformField_name = './data-unif/gk_wham-field_'+str(frame)+'.gkyl'
nonunifField_name = './data-nonunif/gk_wham-field_'+str(frame)+'.gkyl'

uniformPosMap_name = './data-unif/gk_wham-mc2nu_pos.gkyl'
nonunifPosMap_name = './data-nonunif/gk_wham-mc2nu_pos.gkyl'

uniformBiMax_pgdata = pg.GData(uniformBiMax_name)
nonunifBiMax_pgdata = pg.GData(nonunifBiMax_name)
uniformField_pgdata = pg.GData(uniformField_name)
nonunifField_pgdata = pg.GData(nonunifField_name)
uniformPosMap_pgdata = pg.GData(uniformPosMap_name)
nonunifPosMap_pgdata = pg.GData(nonunifPosMap_name)

uniformBiMax_pginterp = pg.GInterpModal(uniformBiMax_pgdata, 1, 'ms')
nonunifBiMax_pginterp = pg.GInterpModal(nonunifBiMax_pgdata, 1, 'ms')
uniformField_pginterp = pg.GInterpModal(uniformField_pgdata, 1, 'ms')
nonunifField_pginterp = pg.GInterpModal(nonunifField_pgdata, 1, 'ms')
uniformPosMap_pginterp = pg.GInterpModal(uniformPosMap_pgdata, 1, 'ms')
nonunifPosMap_pginterp = pg.GInterpModal(nonunifPosMap_pgdata, 1, 'ms')

coords, uniformDens = uniformBiMax_pginterp.interpolate(0)
coords, nonunifDens = nonunifBiMax_pginterp.interpolate(0)

coords, uniformField = uniformField_pginterp.interpolate()
coords, nonunifField = nonunifField_pginterp.interpolate()
coords, uniformPosMap = uniformPosMap_pginterp.interpolate(2)
coords, nonunifPosMap = nonunifPosMap_pginterp.interpolate(2)

plt.figure(figsize=(8, 6))
plt.plot(uniformPosMap[:,0], uniformDens[:,0], label='Uniform grid')
plt.plot(nonunifPosMap[:,0], nonunifDens[:,0], label='Nonuniform grid', linestyle='--')
plt.xlabel("Field line length, normalized")
plt.ylabel("Density, $m^{-3}$")
plt.title("Comparison of ion density on uniform and nonuniform grids")
plt.legend()
plt.tight_layout()
plt.yscale('log')
plt.savefig('./plots/unif_vs_nonunif_density.png', dpi=600)
plt.show()

# Plot the electric potential on uniform and nonuniform grids
plt.figure(figsize=(8, 6))
plt.plot(uniformPosMap[:,0], uniformField[:,0], label='Uniform grid')
plt.plot(nonunifPosMap[:,0], nonunifField[:,0], label='Nonuniform grid', linestyle='--')
plt.xlabel("Field line length, normalized")
plt.ylabel("Electric potential, V")
plt.title("Comparison of electric potential on uniform and nonuniform grids")
plt.legend()
plt.tight_layout()
plt.savefig('./plots/unif_vs_nonunif_field.png', dpi=600)
plt.show()


# if plot_bimax_moms:
#   def make_moms(frame_number):
#     print("Getting moments for frame ", frame_number)

#     filename_ion = str(dataDir+unifFile+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
#     pgData_ion = pg.GData(filename_ion)
#     pgInterp_ion = pg.GInterpModal(pgData_ion, polyOrder, 'ms')
#     coords, n_ion = pgInterp_ion.interpolate(0)
#     coords, u_ion = pgInterp_ion.interpolate(1)
#     coords, Tpar_ion = pgInterp_ion.interpolate(2)
#     coords, Tperp_ion = pgInterp_ion.interpolate(3)

#     filename_field = str(dataDir+unifFile+'-field_'+str(frame_number)+'.gkyl')
#     pgData_field = pg.GData(filename_field)
#     pgInterp_field = pg.GInterpModal(pgData_field, polyOrder, 'ms')
#     coords, phi = pgInterp_field.interpolate()

#     filename_elc_mod = str(dataDir+modifiedFile+'-elc_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
#     pgData_elc_mod = pg.GData(filename_elc_mod)
#     pgInterp_elc_mod = pg.GInterpModal(pgData_elc_mod, polyOrder, 'ms')
#     coords, n_elc_mod = pgInterp_elc_mod.interpolate(0)
#     coords, u_elc_mod = pgInterp_elc_mod.interpolate(1)
#     coords, Tpar_elc_mod = pgInterp_elc_mod.interpolate(2)
#     coords, Tperp_elc_mod = pgInterp_elc_mod.interpolate(3)
    
#     filename_ion_mod = str(dataDir+modifiedFile+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
#     pgData_ion_mod = pg.GData(filename_ion_mod)
#     pgInterp_ion_mod = pg.GInterpModal(pgData_ion_mod, polyOrder, 'ms')
#     coords, n_ion_mod = pgInterp_ion_mod.interpolate(0)
#     coords, u_ion_mod = pgInterp_ion_mod.interpolate(1)
#     coords, Tpar_ion_mod = pgInterp_ion_mod.interpolate(2)
#     coords, Tperp_ion_mod = pgInterp_ion_mod.interpolate(3)
    
#     filename_field_mod = str(dataDir+modifiedFile+'-field_'+str(frame_number)+'.gkyl')
#     pgData_field_mod = pg.GData(filename_field_mod)
#     pgInterp_field_mod = pg.GInterpModal(pgData_field_mod, polyOrder, 'ms')
#     coords, phi_mod = pgInterp_field_mod.interpolate()

#     data = pg.GData(str(dataDir+unifFile+"-nodes.gkyl"))
#     vals = data.get_values()
#     nodes_R = vals[:,0]
#     nodes_Z = vals[:,1]
#     nodes_phi = vals[:,2]

#     shape_R = np.shape(nodes_R)
#     midplane_R_min = nodes_R[shape_R[0]//2]
#     midplane_R_max = nodes_R[shape_R[0]//2]
#     throat_R_min = nodes_R[shape_R[0]//4]
#     throat_R_max = nodes_R[shape_R[0]//4]

#     def expand_1D_array(original_array):
#       new_length = 2 * len(original_array) - 1
#       new_array = np.zeros(new_length)
#       new_array[0] = original_array[0]
#       for i in range(1, len(original_array)):
#           new_array[2*i - 1] = (original_array[i - 1] + original_array[i]) / 2
#           new_array[2*i] = original_array[i]
#       return new_array
    
#     def expand_2D_array(original_array):
#       original_shape = np.shape(original_array)
#       new_shape = (2*original_shape[0]-1, 2*original_shape[1]-1)
#       new_array = np.zeros(new_shape)
#       for i in range(1, original_shape[0]-1):
#         for j in range(1, original_shape[1]-1):
#           new_array[2*i, 2*j] = original_array[i, j]
#           new_array[2*i, 2*j+1] = (original_array[i, j] + original_array[i, j+1]) / 2
#           new_array[2*i, 2*j-1] = (original_array[i, j] + original_array[i, j-1]) / 2
#           new_array[2*i+1, 2*j] = (original_array[i, j] + original_array[i+1, j]) / 2
#           new_array[2*i-1, 2*j] = (original_array[i, j] + original_array[i-1, j]) / 2
#           new_array[2*i+1, 2*j+1] = (original_array[i, j] + original_array[i+1, j+1]) / 2
#           new_array[2*i-1, 2*j-1] = (original_array[i, j] + original_array[i-1, j-1]) / 2
#           new_array[2*i+1, 2*j-1] = (original_array[i, j] + original_array[i+1, j-1]) / 2
#           new_array[2*i-1, 2*j+1] = (original_array[i, j] + original_array[i-1, j+1]) / 2
#       new_array[:,0] = expand_1D_array(original_array[:,0])
#       new_array[:,-1] = expand_1D_array(original_array[:,-1])
#       new_array[0,:] = expand_1D_array(original_array[0,:])
#       new_array[-1,:] = expand_1D_array(original_array[-1,:])
#       return new_array
    
#     def expand_3D_array(original_array):
#       original_shape = np.shape(original_array)
#       new_shape = (2*original_shape[0]-1, 2*original_shape[1]-1, 2*original_shape[2]-1)
#       new_array = np.zeros(new_shape)
#       for i in range(1, original_shape[0]-1):
#         for j in range(1, original_shape[1]-1):
#           for k in range(1, original_shape[2]-1):
#             new_array[2*i, 2*j, 2*k] = original_array[i, j, k]

#             new_array[2*i, 2*j, 2*k+1] = (original_array[i, j, k] + original_array[i, j, k+1]) / 2
#             new_array[2*i, 2*j, 2*k-1] = (original_array[i, j, k] + original_array[i, j, k-1]) / 2
#             new_array[2*i, 2*j+1, 2*k] = (original_array[i, j, k] + original_array[i, j+1, k]) / 2
#             new_array[2*i, 2*j-1, 2*k] = (original_array[i, j, k] + original_array[i, j-1, k]) / 2
#             new_array[2*i+1, 2*j, 2*k] = (original_array[i, j, k] + original_array[i+1, j, k]) / 2
#             new_array[2*i-1, 2*j, 2*k] = (original_array[i, j, k] + original_array[i-1, j, k]) / 2

#             new_array[2*i+1, 2*j+1, 2*k] = (original_array[i, j, k] + original_array[i+1, j+1, k]) / 2
#             new_array[2*i-1, 2*j-1, 2*k] = (original_array[i, j, k] + original_array[i-1, j-1, k]) / 2
#             new_array[2*i+1, 2*j-1, 2*k] = (original_array[i, j, k] + original_array[i+1, j-1, k]) / 2
#             new_array[2*i-1, 2*j+1, 2*k] = (original_array[i, j, k] + original_array[i-1, j+1, k]) / 2

#             new_array[2*i+1, 2*j, 2*k+1] = (original_array[i, j, k] + original_array[i+1, j, k+1]) / 2
#             new_array[2*i-1, 2*j, 2*k-1] = (original_array[i, j, k] + original_array[i-1, j, k-1]) / 2
#             new_array[2*i+1, 2*j, 2*k-1] = (original_array[i, j, k] + original_array[i+1, j, k-1]) / 2
#             new_array[2*i-1, 2*j, 2*k+1] = (original_array[i, j, k] + original_array[i-1, j, k+1]) / 2

#             new_array[2*i, 2*j+1, 2*k+1] = (original_array[i, j, k] + original_array[i, j+1, k+1]) / 2
#             new_array[2*i, 2*j-1, 2*k-1] = (original_array[i, j, k] + original_array[i, j-1, k-1]) / 2
#             new_array[2*i, 2*j+1, 2*k-1] = (original_array[i, j, k] + original_array[i, j+1, k-1]) / 2
#             new_array[2*i, 2*j-1, 2*k+1] = (original_array[i, j, k] + original_array[i, j-1, k+1]) / 2
#       new_array[:,0,:] = expand_2D_array(original_array[:,0,:])
#       new_array[:,-1,:] = expand_2D_array(original_array[:,-1,:])
#       new_array[:,:,0] = expand_2D_array(original_array[:,:,0])
#       new_array[:,:,-1] = expand_2D_array(original_array[:,:,-1])
#       new_array[0,:,:] = expand_2D_array(original_array[0,:,:])
#       new_array[-1,:,:] = expand_2D_array(original_array[-1,:,:])
#       return new_array

    
#     nodes_Z = expand_1D_array(nodes_Z)
#     nodes_Z = nodes_Z[1:]
#     # nodes_R = expand_1D_array(nodes_R)


#     n_elc = n_elc[:,0]
#     u_elc = u_elc[:,0]
#     Tpar_elc = Tpar_elc[:,0] * me / eV
#     Tperp_elc = Tperp_elc[:,0] * me / eV
#     T_elc = (Tpar_elc + 2*Tperp_elc)/3
#     n_ion = n_ion[:,0]
#     u_ion = u_ion[:,0]
#     Tpar_ion = Tpar_ion[:,0] * mi / eV
#     Tperp_ion = Tperp_ion[:,0] * mi / eV
#     T_ion = (Tpar_ion + 2*Tperp_ion)/3
#     phi = phi[:,0]
#     midplane_Te = T_elc[T_elc.shape[0]//2]
#     ephioTe =  phi / midplane_Te

#     n_elc_mod = n_elc_mod[:,0]
#     u_elc_mod = u_elc_mod[:,0]
#     Tpar_elc_mod = Tpar_elc_mod[:,0] * me / eV
#     Tperp_elc_mod = Tperp_elc_mod[:,0] * me / eV
#     T_elc_mod = (Tpar_elc_mod + 2*Tperp_elc_mod)/3
#     n_ion_mod = n_ion_mod[:,0]
#     u_ion_mod = u_ion_mod[:,0]
#     Tpar_ion_mod = Tpar_ion_mod[:,0] * mi / eV
#     Tperp_ion_mod = Tperp_ion_mod[:,0] * mi / eV
#     T_ion_mod = (Tpar_ion_mod + 2*Tperp_ion_mod)/3
#     phi_mod = phi_mod[:,0]
#     midplane_Te_mod = T_elc_mod[T_elc_mod.shape[0]//2]
#     ephioTe_mod =  phi_mod / midplane_Te_mod

#     # # Compute polarization density for ions
#     # # Read in the magnetic field
#     # filename_bmag = str(dataDir+unifFile+'-bmag.gkyl')
#     # pgData_bmag = pg.GData(filename_bmag)
#     # pgInterp_bmag = pg.GInterpModal(pgData_bmag, polyOrder, 'ms')
#     # coords, bmag = pgInterp_bmag.interpolate()
#     # bmag_shape = bmag.shape

#     # # Checked in gkyl_gk_geometry_bmag_mid that this is how it's done
#     # epsilon_i = mi * n_pol / bmag[bmag_shape[0]//2,bmag_shape[1]//2,0]**2

#     # filename_jacobgeo = str(dataDir+unifFile+'-jacobgeo.gkyl')
#     # pgData_jacobgeo = pg.GData(filename_jacobgeo)
#     # pgInterp_jacobgeo = pg.GInterpModal(pgData_jacobgeo, polyOrder, 'ms')
#     # coords, jacobgeo = pgInterp_jacobgeo.interpolate()
#     # jacobgeo = jacobgeo[:,0]

#     # filename_gxx = str(dataDir+unifFile+'-gxxj.gkyl')
#     # pgData_gxx = pg.GData(filename_gxx)
#     # pgInterp_gxx = pg.GInterpModal(pgData_gxx, polyOrder, 'ms')
#     # coords, gxxj = pgInterp_gxx.interpolate()
#     # gxxj = gxxj[:,0]

#     # D =  gxxj * epsilon_i
#     # dpsi = coords[1][0] - coords[0][0]
#     # ni_pol = np.zeros(D.shape)
#     # for i in range (D.shape[0]-1):
#     #   ni_pol[i,:] = -1/jacobgeo[i,:] / dpsi * (((D[i+1,:]+D[i,:])/2) * ((phi[i+1,:]-phi[i,:])/dpsi) - \
#     #                                            ((D[i-1,:]+D[i,:])/2) * ((phi[i,:]-phi[i-1,:])/dpsi))
#     # ni_pol[0,:] = 0.0

#     # make an array grid that is the size of coords

#     X = nodes_Z

#     # X = nodes_Z[:,:]
#     # Y = nodes_R[:,:]

#     # Print where n_ion is nan
#     print(np.argwhere(np.isnan(n_ion)))

    
#     fig, ax = plt.subplots(5, 3, figsize=(12,12))
#     fig.suptitle(str(frame_number*time_per_frame)+' seconds', fontsize=20)

#     def plot_moment_data(data, data_mod, ax, fig, title, locx, locy):
#       ax[locx,locy].plot(X, data, label='Standard collisions')
#       ax[locx,locy].plot(X, data_mod, label='Modified collisions', linestyle='--')
#       ax[locx,locy].set_xlabel('Z cylindrical axis, m')
#       ax[locx,locy].set_ylabel(title)
#       ax[locx,locy].set_title(title, fontsize=16)

#     plot_moment_data(n_elc, n_elc_mod, ax, fig, '$n_e$, $m^{-3}$', 0, 0)
#     plot_moment_data(u_elc, u_elc_mod, ax, fig, '$U_{e,||}$, $m/s$', 0, 2)
#     plot_moment_data(Tpar_elc, Tpar_elc_mod, ax, fig, '$T_{e,||}$, $eV$', 1, 0)
#     plot_moment_data(Tperp_elc, Tperp_elc_mod, ax, fig, '$T_{e,\perp}$, $eV$', 1, 1)
#     plot_moment_data(T_elc, T_elc_mod, ax, fig, '$T_e$, $eV$', 1, 2)
#     plot_moment_data(n_ion, n_ion_mod, ax, fig, '$n_i$, $m^{-3}$', 2, 0)
#     plot_moment_data(u_ion, u_ion_mod, ax, fig, '$U_{i,||}$, $m/s$', 2, 2)
#     plot_moment_data(Tpar_ion, Tpar_ion_mod, ax, fig, '$T_{i,||}$, $eV$', 3, 0)
#     plot_moment_data(Tperp_ion, Tperp_ion_mod, ax, fig, '$T_{i,\perp}$, $eV$', 3, 1)
#     plot_moment_data(T_ion, T_ion_mod, ax, fig, '$T_i$, $eV$', 3, 2)

#     ax[0,2].legend()

#     # Plot electron density on a log scale
#     ax[0,1].plot(X,n_elc, label='Standard collisions')
#     ax[0,1].plot(X,n_elc_mod, label='Modified collisions', linestyle='--')
#     ax[0,1].set_yscale('log')
#     ax[0,1].set_xlabel('Z cylindrical axis, m')
#     ax[0,1].set_ylabel('$n_e$')
#     ax[0,1].set_title('$n_e$ (log scale) $m^{-3}$', fontsize=16)

#     # Plot the ion density on a log scale
#     ax[2,1].plot(X,n_ion, label='Standard collisions')
#     ax[2,1].plot(X,n_ion_mod, label='Modified collisions', linestyle='--')
#     ax[2,1].set_yscale('log')
#     ax[2,1].set_xlabel('Z cylindrical axis, m')
#     ax[2,1].set_ylabel('$n_i$')
#     ax[2,1].set_title('$n_i$ (log scale) $m^{-3}$', fontsize=16)

#     plot_moment_data(phi, phi_mod, ax, fig, '$\phi$, V', 4, 0)
#     plot_moment_data(ephioTe, ephioTe_mod, ax, fig, '$e \phi / T_e$', 4, 1)
#     ax[4,2].remove()

#     plt.tight_layout()
#     plt.savefig(outDir+'moments_'+str(frame_number)+figureFileFormat, dpi=600)
#     plt.close()


