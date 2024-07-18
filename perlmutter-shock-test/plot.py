import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from matplotlib.colors import LogNorm
import multiprocessing
from scipy.integrate import cumulative_trapezoid as cumtrapz

dataDir = './outputs/'
unifFile = 'gk_bgk_periodic_sod_shock_1x2v_p1'
nonunifFile = 'gk_bgk_periodic_sod_shock_1x2v_p1_nonunif'
outDir = './python-plots/'
polyOrder = 1

for frame_number in range(0, 2):
  filename_elc = str(dataDir+unifFile+'-neut_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
  pgData_elc = pg.GData(filename_elc)
  pgInterp_elc = pg.GInterpModal(pgData_elc, polyOrder, 'ms')
  coords, n_unif = pgInterp_elc.interpolate(0)
  coords, u_unif = pgInterp_elc.interpolate(1)
  coords, Tpar_unif = pgInterp_elc.interpolate(2)
  coords, Tperp_unif = pgInterp_elc.interpolate(3)

  filename_nonunif = str(dataDir+nonunifFile+'-neut_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
  pgData_nonunif = pg.GData(filename_nonunif)
  pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'ms')
  coords, n_nonunif = pgInterp_nonunif.interpolate(0)
  coords, u_nonunif = pgInterp_nonunif.interpolate(1)
  coords, Tpar_nonunif = pgInterp_nonunif.interpolate(2)
  coords, Tperp_nonunif = pgInterp_nonunif.interpolate(3)

  filename_mc2p = str(dataDir+nonunifFile+'-mapc2p.gkyl')
  pgData_mc2p = pg.GData(filename_mc2p)
  pgInterp_mc2p = pg.GInterpModal(pgData_mc2p, polyOrder, 'ms')
  x_mc2p, mc2p_z = pgInterp_mc2p.interpolate(2) ## Cylindrical axis Z
  y_mc2p, mc2p_x = pgInterp_mc2p.interpolate(0) ## Cylindrical axis Y


  filename_mc2p = str(dataDir+unifFile+'-mapc2p.gkyl')
  pgData_mc2p = pg.GData(filename_mc2p)
  pgInterp_mc2p = pg.GInterpModal(pgData_mc2p, polyOrder, 'ms')
  x_mc2p, mc2p_z_unif = pgInterp_mc2p.interpolate(2) ## Cylindrical axis Z

  mc2p_z = mc2p_z[:,0]
  unif_z = mc2p_z_unif[:,0]


  plt.plot(mc2p_z, label='Non-uniform')
  plt.plot(unif_z, label='Uniform')
  plt.savefig(outDir+'mc2p_z_'+str(frame_number)+'.png', dpi=300)
  plt.close()

  n_unif = n_unif[:,0]
  n_nonunif = n_nonunif[:,0]
  u_unif = u_unif[:,0]
  u_nonunif = u_nonunif[:,0]
  Tpar_unif = Tpar_unif[:,0]
  Tpar_nonunif = Tpar_nonunif[:,0]
  Tperp_unif = Tperp_unif[:,0]
  Tperp_nonunif = Tperp_nonunif[:,0]

  fig, ax = plt.subplots(2, 2, figsize=(8,8))
  fig.suptitle('Frame ' + str(frame_number), fontsize=14)

  def plot_moment_data(X, data, ax, fig, title, legend, locx, locy, linestye='-'):
    ax[locx,locy].plot(X, data, label=legend, linestyle=linestye)
    ax[locx,locy].set_xlabel('Z cylindrical axis, m')
    ax[locx,locy].set_title(title, fontsize=12)


  plot_moment_data(unif_z, n_unif, ax, fig, 'Density', 'Uniform', 0, 0)
  plot_moment_data(mc2p_z, n_nonunif, ax, fig, 'Density', 'Non-uniform', 0, 0, '--')
  plot_moment_data(unif_z, u_unif, ax, fig, 'Velocity', 'Uniform', 0, 1)
  plot_moment_data(mc2p_z, u_nonunif, ax, fig, 'Velocity', 'Non-uniform', 0, 1, '--')
  plot_moment_data(unif_z, Tpar_unif, ax, fig, 'Tpar', 'Uniform', 1, 0)
  plot_moment_data(mc2p_z, Tpar_nonunif, ax, fig, 'Tpar', 'Non-uniform', 1, 0, '--')
  plot_moment_data(unif_z, Tperp_unif, ax, fig, 'Tperp', 'Uniform', 1, 1)
  plot_moment_data(mc2p_z, Tperp_nonunif, ax, fig, 'Tperp', 'Non-uniform', 1, 1, '--')
  ax[0,0].legend()
  ax[0,1].legend()
  ax[1,0].legend()
  ax[1,1].legend()

  plt.tight_layout()
  plt.savefig(outDir+'moments_'+str(frame_number)+'.png', dpi=300)
  plt.close()