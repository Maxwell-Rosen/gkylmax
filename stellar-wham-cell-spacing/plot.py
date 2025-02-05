
import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from matplotlib.colors import LogNorm
import multiprocessing
from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.optimize import root_scalar
import imageio.v2 as imageio


polyOrder = 1
basisType = 'gkhyb'

eps0, mu0 = 8.8541878176204e-12, 1.2566370614359e-06
eV        = 1.602176487e-19
qe, qi    = -1.602176487e-19, 1.602176487e-19
me, mp    = 9.10938215e-31, 1.672621637e-27

mi        = 2.014*mp                         #[ Deuterium ion mass.
Te0       = 940*eV
n0        = 3.e19
B_p       = 0.53
beta      = 0.4                              #[ Ratio of plasma to magnetic pressure.
tau       = (B_p**2)*beta/(2*mu0*n0*Te0)-1    #[ Ti/Te ratio.
Ti0       = tau*Te0
n_pol     = 3e19

timestep = 8.42244e-12

#[ Thermal speeds.
vti = np.sqrt(Ti0/mi)
vte = np.sqrt(Te0/me)
mui0 = 0.5*mi*(vti**2)/B_p
c_s = np.sqrt(Te0/mi)

#[ Gyrofrequencies and gyroradii.
omega_ci = eV*B_p/mi
rho_s    = c_s/omega_ci

filename_bmag = str('gk_wham-bmag.gkyl')
pgData_bmag = pg.GData(filename_bmag)
pgInterp_bmag = pg.GInterpModal(pgData_bmag, polyOrder, 'ms')
x_bmag, dataOut_bmag = pgInterp_bmag.interpolate()

B = dataOut_bmag[47:145]
z = x_bmag[0][47:145]
shp_B = np.shape(B)
B = B[shp_B[0]//2:,0]
z = z[shp_B[0]//2:]
mid_B = B[0]

filename_elc_vel_grid = str('gk_wham-elc_mapc2p_vel.gkyl')
pgData_elc_vel_grid = pg.GData(filename_elc_vel_grid)
dg_coeffs = pgData_elc_vel_grid._values
c0 = dg_coeffs[:,0,0]
c1 = dg_coeffs[:,0,1]
c2 = dg_coeffs[0,:,2]
c3 = dg_coeffs[0,:,3]

vpar = c0
vpar_shp = np.shape(vpar)
vpar = vpar[vpar_shp[0]//2:]
mu = c2

arr_num_cells = np.zeros([len(mu), len(vpar)])

for im in range(len(mu)):
  mu_max = mu[im]

  Energies = 1/2*me*vpar**2 + mu_max*mid_B
  trapped = Energies - mu_max*B[-1]

  zt = np.zeros(len(vpar))
  for i in range(len(vpar)):
      vp_i = np.abs(vpar[i])
      def minimize_energy(zt_val):
          return 1/2*me*(vp_i**2) + mu_max*mid_B - mu_max*np.interp(zt_val, z, B) 
      try:
          sol = root_scalar(minimize_energy, bracket=[z[0], z[-1]])
          zt[i] = sol.root
      except ValueError:
          zt[i] = np.nan

  dBdz = np.gradient(B, z)
  dBdz_tp = np.interp(zt, z, dBdz)

  dVpar = np.diff(vpar)
  dVpar = np.append(dVpar[0], dVpar)

  Delta_theta = me * vpar * dVpar / (mu_max * dBdz_tp)
  Delta_theta = Delta_theta[~np.isnan(Delta_theta)]
  Lz = 2*np.pi
  min_delta_theta = np.min(Delta_theta)
  n_cells = Lz/min_delta_theta 

  Delta_theta = me * vpar * dVpar / (mu_max * dBdz_tp)
  arr_num_cells[im,:] = Lz/Delta_theta

  print(n_cells)
  
plt.pcolor(vpar, mu, arr_num_cells)
plt.colorbar()
plt.xlabel('$v_{\\parallel}$')
plt.ylabel('$\\mu$')
plt.title('Number of cells required in simulation \n based on phase space compression')
plt.show()


  