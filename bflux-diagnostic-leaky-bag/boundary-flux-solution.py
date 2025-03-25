import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from scipy.special import erf

is_organized = 1

plot_uniform = 1
plot_edge = 1
plot_central = 1

if plot_uniform:
  # Load data
  if is_organized:
    data = pg.GData("gkyldata/gk_boundary_flux_1x2v_p1-ion_integrated_moms.gkyl")
  else:
    data = pg.GData("gk_boundary_flux_1x2v_p1-ion_integrated_moms.gkyl")
  vals = data.get_values()
  time = np.squeeze(np.array(data.get_grid()))

  n0 = 1.0
  vt = 1.0
  L = 1.0

  # Analytic solution
  N = 2*n0/(np.sqrt(2*np.pi)*vt)
  N *= np.sqrt(np.pi/2) * L * vt * erf(L/(np.sqrt(2)*time*vt)) + time*vt*(np.exp(-L**2/(2*time**2*vt**2)) - 1)

  # Plot
  plt.figure(figsize=(10,5))
  plt.suptitle("Uniform grid")
  plt.subplot(1,2,1)
  plt.plot(time, vals[:,0], 'r-')
  plt.plot(time, N, 'k--')
  plt.legend(['Gkeyll', 'Analytic'])
  plt.xlabel('Time')
  plt.ylabel('Integrated density')
  plt.title('Integrated density vs time')
  plt.xlim([np.min(time), np.max(time)])
  plt.ylim([np.min(vals[:,0]), np.max(vals[:,0])])

  plt.subplot(1,2,2)
  plt.plot(time, vals[:,0] - N, 'r-')
  plt.xlabel('Time')
  plt.ylabel('Error')
  plt.title('Error in integrated density vs time')
  plt.xlim([np.min(time), np.max(time)])
  plt.tight_layout()
  plt.savefig('boundary-flux-solution-uniform.png')
  # plt.show()

if plot_edge:
    # Load data
  if is_organized:
    data = pg.GData("gkyldata/gk_boundary_flux_1x2v_p1_nux_edge-ion_integrated_moms.gkyl")
  else:
    data = pg.GData("gk_boundary_flux_1x2v_p1_nux_edge-ion_integrated_moms.gkyl")
  vals = data.get_values()
  time = np.squeeze(np.array(data.get_grid()))

  n0 = 1.0
  vt = 1.0
  L = 1.0

  # Analytic solution
  N = 2*n0/(np.sqrt(2*np.pi)*vt)
  N *= np.sqrt(np.pi/2) * L * vt * erf(L/(np.sqrt(2)*time*vt)) + time*vt*(np.exp(-L**2/(2*time**2*vt**2)) - 1)

  # Plot
  plt.figure(figsize=(10,5))
  plt.suptitle("Cells concentrated at the edge")
  plt.subplot(1,2,1)
  plt.plot(time, vals[:,0], 'r-')
  plt.plot(time, N, 'k--')
  plt.legend(['Gkeyll', 'Analytic'])
  plt.xlabel('Time')
  plt.ylabel('Integrated density')
  plt.title('Integrated density vs time')
  plt.xlim([np.min(time), np.max(time)])
  plt.ylim([np.min(vals[:,0]), np.max(vals[:,0])])

  plt.subplot(1,2,2)
  plt.plot(time, vals[:,0] - N, 'r-')
  plt.xlabel('Time')
  plt.ylabel('Error')
  plt.title('Error in integrated density vs time')
  plt.xlim([np.min(time), np.max(time)])
  plt.tight_layout()
  plt.savefig('boundary-flux-solution-edge.png')
  # plt.show()

if plot_central:
    # Load data
  if is_organized:
    data = pg.GData("gkyldata/gk_boundary_flux_1x2v_p1_nux_center-ion_integrated_moms.gkyl")
  else:
    data = pg.GData("gk_boundary_flux_1x2v_p1_nux_center-ion_integrated_moms.gkyl")
  vals = data.get_values()
  time = np.squeeze(np.array(data.get_grid()))

  n0 = 1.0
  vt = 1.0
  L = 1.0

  # Analytic solution
  N = 2*n0/(np.sqrt(2*np.pi)*vt)
  N *= np.sqrt(np.pi/2) * L * vt * erf(L/(np.sqrt(2)*time*vt)) + time*vt*(np.exp(-L**2/(2*time**2*vt**2)) - 1)

  # Plot
  plt.figure(figsize=(10,5))
  plt.suptitle("Cells concentrated at the center")
  plt.subplot(1,2,1)
  plt.plot(time, vals[:,0], 'r-')
  plt.plot(time, N, 'k--')
  plt.legend(['Gkeyll', 'Analytic'])
  plt.xlabel('Time')
  plt.ylabel('Integrated density')
  plt.title('Integrated density vs time')
  plt.xlim([np.min(time), np.max(time)])
  plt.ylim([np.min(vals[:,0]), np.max(vals[:,0])])

  plt.subplot(1,2,2)
  plt.plot(time, vals[:,0] - N, 'r-')
  plt.xlabel('Time')
  plt.ylabel('Error')
  plt.title('Error in integrated density vs time')
  plt.xlim([np.min(time), np.max(time)])
  plt.tight_layout()
  plt.savefig('boundary-flux-solution-central.png')
  # plt.show()