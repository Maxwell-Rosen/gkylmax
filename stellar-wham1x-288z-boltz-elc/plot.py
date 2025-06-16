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
from scipy.optimize import curve_fit, minimize
from scipy import stats

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression


# dataDir = '/home/mr1884/scratch/Link to scratch_traverse/gkylmax/traverse-wham1x-compare_unif_vs_nonunif/outputs/'
# dataDir = './data-hires-lorad/'
dataDir = './'
simName = 'gk_wham'
origFile = '../stellar-wham1x-288z-run/misc/gk_wham-ion_integrated_moms.gkyl'
frame_max_plus1 = 341
time_per_frame = 10e-6
in_folders_orgnaized = 0 # 1 if it's in the folders, 0 if data is in present directory

plot_potential_trace = 0
plot_bimax_moms = 0
plot_bimax_moms_2D_time_trace = 0
plot_intEnergy_trace = 0
plot_integrate_positivity = 0
plot_stitched_integrated_moms = 0
calculate_linear_interpolation = 0
plot_source_derivative = 0
plot_steady_state_comparison = 1

# frame_arr = np.arange(0,11)
# frame_arr = np.array([1:4])
# save_figure_as_file= 1     #[ If True, save figure to file. If False, display figure on screen.

# plot_moments       = 1  # Plot density, potential, upar, tperp, tpar.
# plot_distvpar      = 0  # plot distribution function in vpar.
# plot_distmu        = 0  # plot distribution function in mu.
# plot_distf_at_z    = 0
# z_loctions = np.array([0, 0.3, 0.98, 2.4])

# print(frame_arr)
# def process_frame(frameNum):
outDir = './python-plots/'

figureFileFormat = '.png'    #[ Can be .png, .pdf, .ps, .eps, .svg.

#   #[ ............... End of user inputs (MAYBE) ..................... ]#

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

#[ Some RGB colors. These are MATLAB-like.
defaultBlue    = [0, 0.4470, 0.7410]
defaultOrange  = [0.8500, 0.3250, 0.0980]
defaultGreen   = [0.4660, 0.6740, 0.1880]
defaultPurple  = [0.4940, 0.1840, 0.5560]
defaultRed     = [0.6350, 0.0780, 0.1840]
defaultSkyBlue = [0.3010, 0.7450, 0.9330]
grey           = [0.5, 0.5, 0.5]
#[ Colors in a single array.
defaultColors = [defaultBlue,defaultOrange,defaultGreen,defaultPurple,defaultRed,defaultSkyBlue,grey,'black']

#[ LineStyles in a single array.
lineStyles = ['-','--',':','-.','None','None','None','None']
markers    = ['None','None','None','None','o','d','s','+']

#[ Some fontsizes used in plots.
xyLabelFontSize       = 17
titleFontSize         = 17
colorBarLabelFontSize = 17
tickFontSize          = 14
legendFontSize        = 14
textFontSize          = 16

#.Set the font size of the ticks to a given size.
def setTickFontSize(axIn,fontSizeIn):
  axIn.tick_params(labelsize = fontSizeIn)

#.Plot vertical lines at +/- given x location.
def plot_verticalLinesPM(xIn, axIn):
  ymin, ymax = axIn.get_ylim()
  eps = 0.5*ymax,
  axIn.plot([xIn, xIn], [ymin-eps, ymax+eps], linestyle=":", color='grey')
  axIn.plot([-xIn, -xIn], [ymin-eps, ymax+eps], linestyle=":", color='grey')
  axIn.set_ylim(ymin, ymax)

#   #................................................................................#

if plot_potential_trace:
  print("Plotting potential trace")
  if in_folders_orgnaized:
    filename_bmag = str(dataDir+'Geometry/'+simName+'-bmag.gkyl')
  else:
    filename_bmag = str(dataDir+simName+'-bmag.gkyl')
  pgData_bmag = pg.GData(filename_bmag)
  pgInterp_bmag = pg.GInterpModal(pgData_bmag, polyOrder, 'ms')
  x_bmag, dataOut_bmag = pgInterp_bmag.interpolate()

  bmag_shape = dataOut_bmag.shape
  midpoint = int(bmag_shape[0]/2)
  upperhalf = dataOut_bmag[midpoint:]
  peak = np.argmax(upperhalf)
  peak_idx = midpoint+peak

  def loadphi(frame_number, filename):
    if in_folders_orgnaized:
      filename_phi = str(dataDir+'Field/'+filename+'-field_'+str(frame_number)+'.gkyl')
    else:
      filename_phi = str(dataDir+filename+'-field_'+str(frame_number)+'.gkyl')
    pgData_phi = pg.GData(filename_phi)
    pgInterp_phi = pg.GInterpModal(pgData_phi, polyOrder, 'ms')
    x_phi, dataOut_phi = pgInterp_phi.interpolate()
    return dataOut_phi
  
  def get_temp(frame_number, filename):
    # filename_elc = str(dataDir+filename+'-elc_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
    # pgData_elc = pg.GData(filename_elc)
    # pgInterp_elc = pg.GInterpModal(pgData_elc, polyOrder, 'ms')
    # coords, Tpar_elc = pgInterp_elc.interpolate(2)
    # coords, Tperp_elc = pgInterp_elc.interpolate(3)
    # Temp = (Tpar_elc[midpoint,0] + 2*Tperp_elc[midpoint,0])/3 * me / eV
    return Te0

  potential = np.zeros(frame_max_plus1)
  expander_potential = np.zeros(frame_max_plus1)
  Temp = np.zeros(frame_max_plus1)
  Temp_mod = np.zeros(frame_max_plus1)  
  for i in range(frame_max_plus1):
    dataOut_phi = loadphi(i, simName)
    Temp[i] = get_temp(i, simName)
    midphi = dataOut_phi[midpoint]
    phi_peak = dataOut_phi[peak_idx]
    potential[i] = (midphi[0] - phi_peak[0]) / Temp[i]
    wall_potential = (dataOut_phi[0][0] + dataOut_phi[-1][0]) / 2
    expander_potential[i] = (phi_peak[0] - wall_potential) / Temp[i]

  Temp = Temp*eV
  potential *= eV
  expander_potential *= eV


  plt.plot(np.arange(frame_max_plus1)*1e-6, potential)
  plt.xlabel('Time, seconds')
  plt.ylabel('Potential difference, $e \phi / T_e(\psi_{min},z=0)$')
  plt.title('Potential difference between midplane and peak magnetic field')
  plt.xscale('log')
  plt.savefig(outDir+'potential_trace'+figureFileFormat)
  plt.close()

  plt.plot(np.arange(frame_max_plus1)*1e-6, expander_potential)
  plt.xlabel('Time, seconds')
  plt.ylabel('Potential difference, $e \phi / T_e(\psi_{min},z=0)$')
  plt.title('Potential difference between mirror throat and wall')
  plt.xscale('log')
  plt.savefig(outDir+'potential_expander_trace'+figureFileFormat)
  plt.close()

  # starting_fit_from_frame = frame_max_plus1 - 20
  # x = np.arange(frame_max_plus1-starting_fit_from_frame)*1e-6
  # y_mod = potential_mod[starting_fit_from_frame:]
  # y_van = potential[starting_fit_from_frame:]
  # # Fit an exponential of the form a - b*exp(-c*x)
  # def fit_func(x, a, b, c):
  #   return a - b*np.exp(-c*x)
  
  # popt_mod, pcov_mod = curve_fit(fit_func, x, y_mod, p0=[5.0, 0.2, 0.0])
  # popt_van, pcov_van = curve_fit(fit_func, x, y_van, p0=[5.0, 0.2, 0.0])

  # print("Fitted parameters for modified collisions: ")
  # print(" a = {:.4f} ± {:.4f} e phi / Te".format(popt_mod[0], np.sqrt(pcov_mod[0,0])))
  # print(" b = {:.4f} ± {:.4f} e phi / Te".format(popt_mod[1], np.sqrt(pcov_mod[1,1])))
  # pct_err = np.sqrt(pcov_mod[2,2]) / popt_mod[2]
  # print(" 1/c = {:.1f} ± {:.1f} microseconds".format(1/(popt_mod[2])*1e6, 1/(popt_mod[2])*pct_err*1e6))

  # print("Fitted parameters for standard collisions: ")
  # print(" a = {:.4f} ± {:.4f} e phi / Te".format(popt_van[0], np.sqrt(pcov_van[0,0])))
  # print(" b = {:.4f} ± {:.4f} e phi / Te".format(popt_van[1], np.sqrt(pcov_van[1,1])))
  # pct_err = np.sqrt(pcov_van[2,2]) / popt_van[2]
  # print(" 1/c = {:.1f} ± {:.1f} microseconds".format(1/(popt_van[2])*1e6, 1/(popt_van[2])*pct_err*1e6))

  # plt.plot(x, y_mod, label='Data Modified Collisions', color='orange', linestyle='-')
  # plt.plot(x, fit_func(x, *popt_mod), label='Fit Modified Collisions: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt_mod), color='goldenrod', linestyle='--')
  # plt.plot(x, y_van, label='Data Standard Collisions', color='blue', linestyle='-')
  # plt.plot(x, fit_func(x, *popt_van), label='Fit Standard Collisions: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt_van), color='skyblue', linestyle='--')
  # plt.xlabel('Time, seconds')
  # plt.ylabel('Potential difference, $e \phi / T_e(\psi_{min},z=0)$')
  # plt.title('Potential difference between midplane and peak magnetic field')
  # plt.legend()

  # plt.savefig(outDir+'potential_trace_fit'+figureFileFormat)
  # plt.close()

  # # Plot the error in this model
  # plt.plot(x, y_mod - fit_func(x, *popt_mod), label='Modified Collisions', color='orange')
  # plt.plot(x, y_van - fit_func(x, *popt_van), label='Standard Collisions', color='blue')
  # plt.xlabel('Time, seconds')
  # plt.ylabel('Error in fit')
  # plt.title('Error in fit of potential difference between midplane and peak magnetic field')
  # plt.legend()
  # plt.savefig(outDir+'potential_trace_fit_error'+figureFileFormat)
  # plt.close()

  # plt.plot(np.arange(frame_max_plus1)*1e-6, Temp/eV, label = 'Standard collisions')
  # plt.plot(np.arange(frame_max_plus1)*1e-6, Temp_mod/eV, linestyle='--', label = 'Modified collisions')
  # plt.xlabel('Time, seconds')
  # plt.ylabel('Temperature, $T_e(\psi_{min},z=0)$')
  # plt.title('Temperature at midplane')
  # plt.legend()
  # plt.savefig(outDir+'temperature_trace'+figureFileFormat)
  # plt.close()

if plot_bimax_moms:
  def make_moms(frame_number):
    print("Getting moments for frame ", frame_number)

    if in_folders_orgnaized:
      filename_ion = str(dataDir+'BiMaxwellianMoments/'+simName+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
    else:
      filename_ion = str(dataDir+simName+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
    pgData_ion = pg.GData(filename_ion)
    pgInterp_ion = pg.GInterpModal(pgData_ion, polyOrder, 'ms')
    coords, n_ion = pgInterp_ion.interpolate(0)
    coords, u_ion = pgInterp_ion.interpolate(1)
    coords, Tpar_ion = pgInterp_ion.interpolate(2)
    coords, Tperp_ion = pgInterp_ion.interpolate(3)

    if in_folders_orgnaized:
      filename_field = str(dataDir+'Field/'+simName+'-field_'+str(frame_number)+'.gkyl')
    else:
      filename_field = str(dataDir+simName+'-field_'+str(frame_number)+'.gkyl')
    pgData_field = pg.GData(filename_field)
    pgInterp_field = pg.GInterpModal(pgData_field, polyOrder, 'ms')
    coords, phi = pgInterp_field.interpolate()
    
    if in_folders_orgnaized:
      filename_mc2nu_pos = pg.GData(str(dataDir+'Geometry/'+simName+'-mc2nu_pos.gkyl'))
    else:
      filename_mc2nu_pos = pg.GData(str(dataDir+simName+'-mc2nu_pos.gkyl'))
    interp = pg.GInterpModal(filename_mc2nu_pos, 1, 'ms')
    nodes_Z = interp.interpolate(2)[1]
    nodes_Z = np.squeeze(nodes_Z)

    n_ion = n_ion[:,0]
    u_ion = u_ion[:,0]
    Tpar_ion = Tpar_ion[:,0] * mi / eV
    Tperp_ion = Tperp_ion[:,0] * mi / eV
    T_ion = (Tpar_ion + 2*Tperp_ion)/3
    phi = phi[:,0]
    midplane_Te = Te0
    ephioTe =  eV * phi / midplane_Te

    X = nodes_Z

    # Print where n_ion is na
    print(np.argwhere(np.isnan(n_ion)))
    
    fig, ax = plt.subplots(3, 3, figsize=(12,12))
    fig.suptitle(str(frame_number*time_per_frame)+' seconds', fontsize=20)

    def plot_moment_data(data, ax, fig, title, locx, locy):
      ax[locx,locy].plot(X, data)
      ax[locx,locy].set_xlabel('Field line length, radians')
      ax[locx,locy].set_ylabel(title)
      ax[locx,locy].set_title(title, fontsize=16)

    plot_moment_data(n_ion, ax, fig, '$n_i$, $m^{-3}$', 0, 0)
    plot_moment_data(u_ion, ax, fig, '$U_{i,||}$, $m/s$', 0, 2)
    plot_moment_data(Tpar_ion, ax, fig, '$T_{i,||}$, $eV$', 1, 0)
    plot_moment_data(Tperp_ion, ax, fig, '$T_{i,\perp}$, $eV$', 1, 1)
    plot_moment_data(T_ion, ax, fig, '$T_i$, $eV$', 1, 2)

    # Plot the ion density on a log scale
    ax[0,1].plot(X, n_ion)
    ax[0,1].set_yscale('log')
    ax[0,1].set_xlabel('Field line length, radians')
    ax[0,1].set_ylabel('$n_i$')
    ax[0,1].set_title('$n_i$ (log scale) $m^{-3}$', fontsize=16)

    plot_moment_data(phi, ax, fig, '$\phi$, V', 2, 0)
    plot_moment_data(ephioTe, ax, fig, '$e \phi / T_e$', 2, 1)
    ax[2,2].remove()

    plt.tight_layout()
    plt.savefig(outDir+'moments_'+str(frame_number)+figureFileFormat, dpi=600)
    plt.close()

  # Number of processes to run in parallel
  # make_moms(0)
  frame_arr = np.arange(0,frame_max_plus1)
  num_processes = multiprocessing.cpu_count()
  print('Number of processes: ', num_processes)
  pool = multiprocessing.Pool(processes=num_processes)
  pool.map(make_moms, frame_arr)
  pool.close()
  pool.join()


  # # Define the filenames in order
  # filenames = [f'moments_{i}.png' for i in range(0, frame_max_plus1)]
  # filenames = [outDir+f'moments_{i}.png' for i in range(0, frame_max_plus1)]

  # # Create a writer object specifying the output file name and frame rate
  # with imageio.get_writer(outDir+'moments_movie.mp4', mode='I', fps=5) as writer:
  #     for filename in filenames:
  #         image = imageio.imread(filename)
  #         writer.append_data(image)
  # print("Movie created successfully.")

if plot_bimax_moms_2D_time_trace:
    def get_moms(frame_number):
      if in_folders_orgnaized:
        filename_ion = str(dataDir+'BiMaxwellianMoments/'+simName+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
      else:
        filename_ion = str(dataDir+simName+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
      pgData_ion = pg.GData(filename_ion)
      pgInterp_ion = pg.GInterpModal(pgData_ion, polyOrder, 'ms')
      _, n = pgInterp_ion.interpolate(0)
      _, upar = pgInterp_ion.interpolate(1)
      _, Tpar = pgInterp_ion.interpolate(2)
      _, Tperp = pgInterp_ion.interpolate(3)

      if in_folders_orgnaized:
        filename_field = str(dataDir+'Field/'+simName+'-field_'+str(frame_number)+'.gkyl')
      else:
        filename_field = str(dataDir+simName+'-field_'+str(frame_number)+'.gkyl')
      pgData_field = pg.GData(filename_field)
      pgInterp_field = pg.GInterpModal(pgData_field, polyOrder, 'ms')
      _, phi = pgInterp_field.interpolate()

      n = n[:,0]
      upar = upar[:,0]
      Tpar = Tpar[:,0]
      Tperp = Tperp[:,0]
      phi = phi[:,0]
      return n, upar, Tpar, Tperp, phi
    
    if in_folders_orgnaized:
      filename_mc2nu_pos = pg.GData(str(dataDir+'Geometry/'+simName+'-mc2nu_pos.gkyl'))
    else:
      filename_mc2nu_pos = pg.GData(str(dataDir+simName+'-mc2nu_pos.gkyl'))
    interp = pg.GInterpModal(filename_mc2nu_pos, 1, 'ms')
    nodes_Z = interp.interpolate(2)[1]
    nodes_Z = np.squeeze(nodes_Z)

    time = np.arange(frame_max_plus1+1)*time_per_frame


    n_ion = np.zeros((frame_max_plus1, nodes_Z.shape[0]))
    u_ion = np.zeros((frame_max_plus1, nodes_Z.shape[0]))
    Tpar_ion = np.zeros((frame_max_plus1, nodes_Z.shape[0]))
    Tperp_ion = np.zeros((frame_max_plus1, nodes_Z.shape[0]))
    phi = np.zeros((frame_max_plus1, nodes_Z.shape[0]))

    def process_frame(frame_number):
      n, upar, Tpar, Tperp, phi_frame = get_moms(frame_number)
      return frame_number, n, upar, Tpar, Tperp, phi_frame

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
      results = pool.map(process_frame, range(frame_max_plus1))

    for frame_number, n, upar, Tpar, Tperp, phi_frame in results:
      n_ion[frame_number, :] = n
      u_ion[frame_number, :] = upar
      Tpar_ion[frame_number, :] = Tpar
      Tperp_ion[frame_number, :] = Tperp
      phi[frame_number, :] = phi_frame

    Tpar_ion = Tpar_ion * mi / eV
    Tperp_ion = Tperp_ion * mi / eV
    T_ion = (Tpar_ion + 2*Tperp_ion)/3
    Energy = n_ion * T_ion
    ephioTe = eV * phi / Te0

    z_lower = -np.pi + 0.1
    z_upper = np.pi - 0.1
    z_lefts = nodes_Z[0:-1]
    z_rights = nodes_Z[1:]
    z_avg = (z_lefts + z_rights) / 2
    z_edges = np.append(z_lower, z_avg)
    z_edges = np.append(z_edges, z_upper)

    fig, ax = plt.subplots(3, 3, figsize=(20,12))
    fig.suptitle('Time trace of moments', fontsize=20)

    ax[0,0].pcolormesh(time, z_edges, n_ion.T, shading='auto', cmap='inferno')
    ax[0,0].set_xlabel('Time, seconds')
    ax[0,0].set_ylabel('Field line length, radians')
    cbar = plt.colorbar(ax[0,0].pcolormesh(time, z_edges, n_ion.T, shading='auto', cmap='inferno'))
    cbar.set_label('$n_i$, $m^{-3}$')
    ax[0,0].set_title('$n_i$, $m^{-3}$')

    ax[0,1].pcolormesh(time, z_edges, n_ion.T, shading='auto', norm=LogNorm(), cmap='inferno')
    ax[0,1].set_xlabel('Time, seconds')
    ax[0,1].set_ylabel('Field line length, radians')
    cbar = plt.colorbar(ax[0,1].pcolormesh(time, z_edges, n_ion.T, shading='auto', norm=LogNorm(), cmap='inferno'))
    cbar.set_label('$n_i$, $m^{-3}$')
    ax[0,1].set_title('$n_i$, $m^{-3}$')

    ax[0,2].pcolormesh(time, z_edges, u_ion.T, shading='auto', cmap='inferno')
    ax[0,2].set_xlabel('Time, seconds')
    ax[0,2].set_ylabel('Field line length, radians')
    cbar = plt.colorbar(ax[0,2].pcolormesh(time, z_edges, u_ion.T, shading='auto', cmap='inferno'))
    cbar.set_label('$U_{i,||}$, $m/s$')
    ax[0,2].set_title('$U_{i,||}$, $m/s$')

    ax[1,0].pcolormesh(time, z_edges, Tpar_ion.T, shading='auto', cmap='inferno')
    ax[1,0].set_xlabel('Time, seconds')
    ax[1,0].set_ylabel('Field line length, radians')
    cbar = plt.colorbar(ax[1,0].pcolormesh(time, z_edges, Tpar_ion.T, shading='auto', cmap='inferno'))
    cbar.set_label('$T_{i,||}$, $eV$')
    ax[1,0].set_title('$T_{i,||}$, $eV$')

    ax[1,1].pcolormesh(time, z_edges, Tperp_ion.T, shading='auto', cmap='inferno')
    ax[1,1].set_xlabel('Time, seconds')
    ax[1,1].set_ylabel('Field line length, radians')
    cbar = plt.colorbar(ax[1,1].pcolormesh(time, z_edges, Tperp_ion.T, shading='auto', cmap='inferno'))
    cbar.set_label('$T_{i,\perp}$, $eV$')
    ax[1,1].set_title('$T_{i,\perp}$, $eV$')
    
    ax[1,2].pcolormesh(time, z_edges, T_ion.T, shading='auto', cmap='inferno')
    ax[1,2].set_xlabel('Time, seconds')
    ax[1,2].set_ylabel('Field line length, radians')
    cbar = plt.colorbar(ax[1,2].pcolormesh(time, z_edges, T_ion.T, shading='auto', cmap='inferno'))
    cbar.set_label('$T_i$, $eV$')
    ax[1,2].set_title('$T_i$, $eV$')

    ax[2,0].pcolormesh(time, z_edges, phi.T, shading='auto', cmap='inferno')
    ax[2,0].set_xlabel('Time, seconds')
    ax[2,0].set_ylabel('Field line length, radians')
    cbar = plt.colorbar(ax[2,0].pcolormesh(time, z_edges, phi.T, shading='auto', cmap='inferno'))
    cbar.set_label('$\phi$, V')
    ax[2,0].set_title('$\phi$, V')

    ax[2,1].pcolormesh(time, z_edges, ephioTe.T, shading='auto', cmap='inferno')
    ax[2,1].set_xlabel('Time, seconds')
    ax[2,1].set_ylabel('Field line length, radians')
    cbar = plt.colorbar(ax[2,1].pcolormesh(time, z_edges, ephioTe.T, shading='auto', cmap='inferno'))
    cbar.set_label('$e \phi / T_e$')
    ax[2,1].set_title('$e \phi / T_e$')

    ax[2,2].pcolor(time, z_edges, Energy.T, shading='auto', cmap='inferno')
    ax[2,2].set_xlabel('Time, seconds')
    ax[2,2].set_ylabel('Field line length, radians')
    cbar = plt.colorbar(ax[2,2].pcolor(time, z_edges, Energy.T, shading='auto', cmap='inferno'))
    cbar.set_label('Energy, $eV m^{-3}$')
    ax[2,2].set_title('Energy, $eV m^{-3}$')

    plt.tight_layout()
    plt.savefig(outDir+'time_trace_moments_n'+figureFileFormat, dpi=600)
    plt.close()


if plot_intEnergy_trace:
    print("Getting integrated energy trace")
    def get_moms(frame_number):

      if in_folders_orgnaized:
        filename_ion = str(dataDir+'BiMaxwellianMoments/'+simName+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
      else:
        filename_ion = str(dataDir+simName+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
      pgData_ion = pg.GData(filename_ion)
      pgInterp_ion = pg.GInterpModal(pgData_ion, polyOrder, 'ms')
      _, n = pgInterp_ion.interpolate(0)
      _, upar = pgInterp_ion.interpolate(1)
      _, Tpar = pgInterp_ion.interpolate(2)
      _, Tperp = pgInterp_ion.interpolate(3)

      if in_folders_orgnaized:
        filename_field = str(dataDir+'Field/'+simName+'-field_'+str(frame_number)+'.gkyl')
      else:
        filename_field = str(dataDir+simName+'-field_'+str(frame_number)+'.gkyl')
      pgData_field = pg.GData(filename_field)
      pgInterp_field = pg.GInterpModal(pgData_field, polyOrder, 'ms')
      _, phi = pgInterp_field.interpolate()

      n = n[:,0]
      upar = upar[:,0]
      Tpar = Tpar[:,0]
      Tperp = Tperp[:,0]
      phi = phi[:,0]
      return n, upar, Tpar, Tperp, phi
    
    if in_folders_orgnaized:
      filename_mc2nu_pos = pg.GData(str(dataDir+'Geometry/'+simName+'-mc2nu_pos.gkyl'))
    else:
      filename_mc2nu_pos = pg.GData(str(dataDir+simName+'-mc2nu_pos.gkyl'))
    interp = pg.GInterpModal(filename_mc2nu_pos, 1, 'ms')
    nodes_Z = interp.interpolate(2)[1]
    nodes_Z = np.squeeze(nodes_Z)

    time = np.arange(frame_max_plus1+1)*time_per_frame


    n_ion = np.zeros((frame_max_plus1, nodes_Z.shape[0]))
    u_ion = np.zeros((frame_max_plus1, nodes_Z.shape[0]))
    Tpar_ion = np.zeros((frame_max_plus1, nodes_Z.shape[0]))
    Tperp_ion = np.zeros((frame_max_plus1, nodes_Z.shape[0]))
    phi = np.zeros((frame_max_plus1, nodes_Z.shape[0]))

    def process_frame(frame_number):
      n, upar, Tpar, Tperp, phi_frame = get_moms(frame_number)
      return frame_number, n, upar, Tpar, Tperp, phi_frame

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
      results = pool.map(process_frame, range(frame_max_plus1))

    for frame_number, n, upar, Tpar, Tperp, phi_frame in results:
      n_ion[frame_number, :] = n
      u_ion[frame_number, :] = upar
      Tpar_ion[frame_number, :] = Tpar
      Tperp_ion[frame_number, :] = Tperp
      phi[frame_number, :] = phi_frame

    Tpar_ion = Tpar_ion * mi / eV
    Tperp_ion = Tperp_ion * mi / eV
    T_ion = (Tpar_ion + 2*Tperp_ion)/3
    Energy = n_ion * T_ion
    ephioTe = eV * phi / Te0

    intEnergy = np.trapz(Energy, x=nodes_Z, axis=1)
    plt.figure(figsize=(12,8))
    plt.plot(time[:-1]*1e3, intEnergy)
    plt.xlabel('Time, ms')
    plt.ylabel('Integrated energy $\int n T dx$, $eV m^{-3}$')
    plt.title('Integrated energy $\int n T dx$')
    plt.savefig(outDir+'integrated_energy_trace'+figureFileFormat)
    plt.close()

  
if plot_integrate_positivity:
    print("Getting integrated moments")
    
#  pgkyl "$name"-"$species"_integrated_moms.gkyl -t f\
#  "$name"-"$species"_positivity_shift_integrated_moms.gkyl -t p \
#  activate -t f,p ev -t poverf 'p f /' \
#  activate -t poverf pl --title "Mp/Mf" --saveas "$saveLoc-positivity-moms-over-f.png" --no-show&

    if in_folders_orgnaized:
      filename_ion = str(dataDir+'misc/'+simName+'-ion_integrated_moms.gkyl')
    else:
      filename_ion = str(dataDir + simName + '-ion_integrated_moms.gkyl')
    pgData_ion = pg.GData(filename_ion)
    M_ion = pgData_ion.get_values()
    M0_ion = np.array(M_ion[:,0])
    M1_ion = np.array(M_ion[:,1])
    M2par_ion = np.array(M_ion[:,2])
    M2perp_ion = np.array(M_ion[:,3])
    time = np.squeeze(np.array(pgData_ion.get_grid()))

    n_ion = M0_ion
    u_ion = M1_ion / M0_ion
    Tpar_ion = (M2par_ion - u_ion / M0_ion) * mi / eV / M0_ion
    Tperp_ion = M2perp_ion / M0_ion * mi / eV / 2.0
    T_ion = (Tpar_ion + 2*Tperp_ion)/3

    if in_folders_orgnaized:
      filename_ion_positivity = str(dataDir+'misc/'+simName+'-ion_positivity_shift_integrated_moms.gkyl')
    else:
      filename_ion_positivity = str(dataDir+simName+'-ion_positivity_shift_integrated_moms.gkyl')
    pgData_ion_positivity = pg.GData(filename_ion_positivity)
    M_ion_positivity = pgData_ion_positivity.get_values()
    M0_ion_positivity = np.array(M_ion_positivity[:,0])
    M1_ion_positivity = np.array(M_ion_positivity[:,1])
    M2par_ion_positivity = np.array(M_ion_positivity[:,2])
    M2perp_ion_positivity = np.array(M_ion_positivity[:,3])

    n_ion_positivity = M0_ion_positivity
    u_ion_positivity = np.divide(M1_ion_positivity, M0_ion_positivity, where=M0_ion_positivity!=0)
    Tpar_ion_positivity = (M2par_ion_positivity - np.divide(u_ion_positivity, M0_ion_positivity, where=M0_ion_positivity!=0)) * mi / eV / M0_ion
    Tperp_ion_positivity = M2perp_ion_positivity / M0_ion * mi / eV / 2.0
    T_ion_positivity = (Tpar_ion_positivity + 2*Tperp_ion_positivity)/3

    fig, ax = plt.subplots(4, 3, figsize=(12,20))
    fig.suptitle('Integrated moments', fontsize=20)

    def plot_moment_data(data, ax, fig, title, locx, locy):
      ax[locx,locy].plot(time, data)
      ax[locx,locy].set_xlabel('Time, seconds')
      ax[locx,locy].set_title(title)

    plot_moment_data(n_ion, ax, fig, '$n_i$, $m^{-3}$', 0, 0)
    plot_moment_data(u_ion, ax, fig, '$U_{i,||}$, $m/s$', 0, 2)
    plot_moment_data(Tpar_ion, ax, fig, '$T_{i,||}$, $eV$', 1, 0)
    plot_moment_data(Tperp_ion, ax, fig, '$T_{i,\perp}$, $eV$', 1, 1)
    plot_moment_data(T_ion, ax, fig, '$T_i$, $eV$', 1, 2)

    plot_moment_data(n_ion_positivity, ax, fig, 'Positivity $n_i$, $m^{-3}$', 2, 0)
    plot_moment_data(u_ion_positivity, ax, fig, 'Positivity $U_{i,||}$, $m/s$', 2, 2)
    plot_moment_data(Tpar_ion_positivity, ax, fig, 'Positivity $T_{i,||}$, $eV$', 3, 0)
    plot_moment_data(Tperp_ion_positivity, ax, fig, 'Positivity $T_{i,\perp}$, $eV$', 3, 1)
    plot_moment_data(T_ion_positivity, ax, fig, 'Positivity $T_i$, $eV$', 3, 2)

    plt.tight_layout()
    plt.savefig(outDir+'integrated_moments'+figureFileFormat, dpi=600)
    plt.close()

    ##########################################################################################

    fig, ax = plt.subplots(4, 2, figsize=(12,20))
    fig.suptitle('Integrated Ms', fontsize=20)

    plot_moment_data(M0_ion, ax, fig, 'M0 ion', 0, 0)
    plot_moment_data(M1_ion, ax, fig, 'M1 ion', 0, 1)
    plot_moment_data(M2par_ion, ax, fig, 'M2par ion', 1, 0)
    plot_moment_data(M2perp_ion, ax, fig, 'M2perp ion', 1, 1)

    plot_moment_data(M0_ion_positivity, ax, fig, 'M0 ion positivity', 2, 0)
    plot_moment_data(M1_ion_positivity, ax, fig, 'M1 ion positivity', 2, 1)
    plot_moment_data(M2par_ion_positivity, ax, fig, 'M2par ion positivity', 3, 0)
    plot_moment_data(M2perp_ion_positivity, ax, fig, 'M2perp ion positivity', 3, 1)

    plt.tight_layout()
    plt.savefig(outDir+'integrated_Ms'+figureFileFormat, dpi=600)
    plt.close()

    M0_ion_ratio = M0_ion_positivity / M0_ion
    M1_ion_ratio = M1_ion_positivity / M1_ion
    M2par_ion_ratio = M2par_ion_positivity / M2par_ion
    M2perp_ion_ratio = M2perp_ion_positivity / M2perp_ion

    fig, ax = plt.subplots(2, 2, figsize=(12,10))
    fig.suptitle('Ratios of $M_{i,positivity} / M_i$', fontsize=20)

    plot_moment_data(M0_ion_ratio, ax, fig, 'M0 ion ratio', 0, 0)
    plot_moment_data(M1_ion_ratio, ax, fig, 'M1 ion ratio', 0, 1)
    plot_moment_data(M2par_ion_ratio, ax, fig, 'M2par ion ratio', 1, 0)
    plot_moment_data(M2perp_ion_ratio, ax, fig, 'M2perp ion ratio', 1, 1)

    plt.tight_layout()
    plt.savefig(outDir+'integrated_Ms_ratios'+figureFileFormat, dpi=600)
    plt.close()

    M0_ion_ratio_total = cumtrapz(M0_ion_ratio, time, initial=0) / timestep
    M1_ion_ratio_total = cumtrapz(M1_ion_ratio, time, initial=0) / timestep
    M2par_ion_ratio_total = cumtrapz(M2par_ion_ratio, time, initial=0) / timestep
    M2perp_ion_ratio_total = cumtrapz(M2perp_ion_ratio, time, initial=0) / timestep

    fig, ax = plt.subplots(2, 2, figsize=(12,10))
    fig.suptitle('Time integrated positivity ratios', fontsize=20)

    plot_moment_data(M0_ion_ratio_total, ax, fig, 'M0 ion ratio', 0, 0)
    plot_moment_data(M1_ion_ratio_total, ax, fig, 'M1 ion ratio', 0, 1)
    plot_moment_data(M2par_ion_ratio_total, ax, fig, 'M2par ion ratio', 1, 0)
    plot_moment_data(M2perp_ion_ratio_total, ax, fig, 'M2perp ion ratio', 1, 1)

    plt.tight_layout()
    plt.savefig(outDir+'integrated_Ms_ratios_in_time'+figureFileFormat, dpi=600)
    plt.close()

    positivity_M0_total = cumtrapz(M0_ion_positivity, time, initial=0) / timestep
    positivity_M1_total = cumtrapz(M1_ion_positivity, time, initial=0) / timestep
    positivity_M2par_total = cumtrapz(M2par_ion_positivity, time, initial=0) / timestep
    positivity_M2perp_total = cumtrapz(M2perp_ion_positivity, time, initial=0) / timestep

    fig, ax = plt.subplots(2, 2, figsize=(12,10))
    fig.suptitle('Time integrated positivity', fontsize=20)
    plot_moment_data(positivity_M0_total, ax, fig, 'M0 ion positivity', 0, 0)
    # plot_moment_data(M0_ion, ax, fig, 'M0 ion', 0, 0)
    plot_moment_data(positivity_M1_total, ax, fig, 'M1 ion positivity', 0, 1)
    # plot_moment_data(M1_ion, ax, fig, 'M1 ion', 0, 1)
    plot_moment_data(positivity_M2par_total, ax, fig, 'M2par ion positivity', 1, 0)
    # plot_moment_data(M2par_ion, ax, fig, 'M2par ion', 1, 0)
    plot_moment_data(positivity_M2perp_total, ax, fig, 'M2perp ion positivity', 1, 1)
    # plot_moment_data(M2perp_ion, ax, fig, 'M2perp ion', 1, 1)

    plt.tight_layout()
    plt.savefig(outDir+'integrated_Ms_positivity_in_time'+figureFileFormat, dpi=600)
    plt.close()

if plot_stitched_integrated_moms:
  print("Getting integrated moments")
  if in_folders_orgnaized:
    filename_ion = str(dataDir+'misc/'+simName+'-ion_integrated_moms.gkyl')
  else:
    filename_ion = str(dataDir + simName + '-ion_integrated_moms.gkyl')
  pgData_ion = pg.GData(filename_ion)
  M_ion = pgData_ion.get_values()
  M0_ion = np.array(M_ion[:,0])
  M1_ion = np.array(M_ion[:,1])
  M2par_ion = np.array(M_ion[:,2])
  M2perp_ion = np.array(M_ion[:,3])
  time = np.squeeze(np.array(pgData_ion.get_grid()))

  orig_integ_moms_data = pg.GData(origFile)
  orig_integ_moms = orig_integ_moms_data.get_values()
  orig_time = np.squeeze(np.array(orig_integ_moms_data.get_grid()))
  orig_M0_ion = np.array(orig_integ_moms[:,0])
  orig_M1_ion = np.array(orig_integ_moms[:,1])
  orig_M2par_ion = np.array(orig_integ_moms[:,2])
  orig_M2perp_ion = np.array(orig_integ_moms[:,3])

  continued_time = time + orig_time[-1]

  fig, ax = plt.subplots(2,2, figsize=(12,12))
  fig.suptitle('Integrated moments', fontsize=20)

  ax[0,0].plot(orig_time, orig_M0_ion)
  ax[0,0].plot(continued_time, M0_ion)
  ax[0,0].set_title('M0')
  ax[0,0].set_xlabel('Time')
  ax[0,0].set_ylabel('Integrated M0')

  ax[0,1].plot(orig_time, orig_M1_ion)
  ax[0,1].plot(continued_time, M1_ion)
  ax[0,1].set_title('M1')
  ax[0,1].set_xlabel('Time')
  ax[0,1].set_ylabel('Integrated M1')

  ax[1,0].plot(orig_time, orig_M2par_ion)
  ax[1,0].plot(continued_time, M2par_ion)
  ax[1,0].set_title('M2 parallel')
  ax[1,0].set_xlabel('Time')
  ax[1,0].set_ylabel('Integrated M2 parallel')

  ax[1,1].plot(orig_time, orig_M2perp_ion)
  ax[1,1].plot(continued_time, M2perp_ion)
  ax[1,1].set_title('M2 perpendicular')
  ax[1,1].set_xlabel('Time')
  ax[1,1].set_ylabel('Integrated M2 perpendicular')

  plt.savefig(outDir + 'continued_integrated_moms.png')
  plt.close()

  # Get the derivatives
  dM0_ion = np.gradient(M0_ion, continued_time)
  dM1_ion = np.gradient(M1_ion, continued_time)
  dM2par_ion = np.gradient(M2par_ion, continued_time)
  dM2perp_ion = np.gradient(M2perp_ion, continued_time)

  orig_dM0_ion = np.gradient(orig_M0_ion, orig_time)
  orig_dM1_ion = np.gradient(orig_M1_ion, orig_time)
  orig_dM2par_ion = np.gradient(orig_M2par_ion, orig_time)
  orig_dM2perp_ion = np.gradient(orig_M2perp_ion, orig_time)

  fig, ax = plt.subplots(2,2, figsize=(12,12))
  fig.suptitle('Derivatives of integrated moments', fontsize=20)

  ax[0,0].plot(time, dM0_ion)
  # ax[0,0].set_xlim(0.1e-3, continued_time[-1])
  ax[0,0].set_title('dM0/dt')
  ax[0,0].set_xlabel('Time')
  ax[0,0].set_ylabel('Derivative of integrated M0')

  ax[0,1].plot(time, dM1_ion)
  # ax[0,1].set_xlim(0.1e-3, continued_time[-1])
  ax[0,1].set_title('dM1/dt')
  ax[0,1].set_xlabel('Time')
  ax[0,1].set_ylabel('Derivative of integrated M1')

  ax[1,0].plot(time, dM2par_ion)
  # ax[1,0].set_xlim(0.1e-3, continued_time[-1])
  ax[1,0].set_title('dM2par/dt')
  ax[1,0].set_xlabel('Time')
  ax[1,0].set_ylabel('Derivative of integrated M2 parallel')

  ax[1,1].plot(time, dM2perp_ion)
  # ax[1,1].set_xlim(0.1e-3, continued_time[-1])
  ax[1,1].set_title('dM2perp/dt')
  ax[1,1].set_xlabel('Time')
  ax[1,1].set_ylabel('Derivative of integrated M2 perpendicular')

  dM0_ion = np.array(dM0_ion)
  dM2_ion = np.gradient((M2par_ion + 2*M2perp_ion)/3, continued_time)
  # for i in range(len(dM0_ion)):
  #   # print(dM0_ion[i])
  #   print(dM0_ion[i])

  # print(dM0_ion[-10:])
  # print(dM2par_ion[-10:])
  # print(dM2perp_ion[-10:])



  plt.savefig(outDir + 'continued_integrated_moms_derivatives.png')
  plt.close()

  # Take the last few thousand points and fit an exponential A * exp ( - B * t) + C
  # to the data to get the decay rate of the moments

  # Get the last 1000 points
  print(len(time))
  time_fit = time[2500:] * 1e5
  M0_ion_fit = M0_ion[2500:] * 1e-19

  plt.figure()
  plt.plot(time_fit, M0_ion_fit)
  plt.xlabel('Time, $\mu s$')
  plt.ylabel('Integrated M0')
  plt.savefig(outDir + 'continued_integrated_moms_fit_data.png')
  plt.show()

  plt.figure()
  plt.plot(time_fit, np.gradient(M0_ion_fit, time_fit))
  plt.xlabel('Time, $\mu s$')
  plt.ylabel('dM0/dt')
  plt.savefig(outDir + 'continued_integrated_moms_fit_data_derivative.png')
  plt.show()


  print(time_fit[0], time_fit[-1])
  print(M0_ion_fit[0], M0_ion_fit[-1])

  # Fit the data
  def exp_func(t, A, B, C):
    return A * np.exp(-B * t) + C
  
  popt, pcov = curve_fit(exp_func, time_fit, M0_ion_fit, p0=[-1, 1e-2, 1.2])
  print("Fitted parameters: ", popt)
  print("Covariance: ", pcov)

  # plt.figure()
  # plt.plot(time_fit, M0_ion_fit, label='Data')
  # plt.plot(time_fit, exp_func(time_fit, *popt), label='Fit')
  # plt.xlabel('Time, $\mu s$')
  # plt.ylabel('Integrated M0')
  # plt.legend()
  # plt.savefig(outDir + 'continued_integrated_moms_fit.png')
  # plt.show()


  # Aiming for total M0 of _

  # Source value of 5.167229357643037 gives estimated total M0 of 5.8960673300757049e18

  # source = np.array([5.167229357643037, 8, 50, 104.99260584034226, 120]) # e20

  # total_M0 = np.array([5.8960673300757049, 6.4275217249160654, 8.7516599139106610, 11.113950849009511, 11.601065890663596]) # * 1e18

if calculate_linear_interpolation:
  #nu_frac = 2000
  # Match values of M0 at 200 microseconds
  source = np.array([130, 140, 90, 100 ])
  total_M0 = np.array([1.696, 1.7705, 1.386, 1.465])

  # Match values of M0 at 700 microseconds
  source = np.array([66.5, 75])
  total_M0 = np.array([1.41934, 1.52578])

  # source = source[-2:]
  # total_M0 = total_M0[-2:]

  # Fit the data
  slope, intercept, r_value, p_value, std_err = stats.linregress(source, total_M0)
  print("Density slope: ", slope)
  print("Density intercept: ", intercept)

  desired_M0 = 1.206625073 # * 1e18

  source_value = (desired_M0 - intercept) / slope
  print("Source value: ", source_value)

  plt.figure()
  plt.plot(source, total_M0, 'o')
  plt.xlabel('Source value')
  plt.ylabel('Total M0')
  plt.show()

if plot_source_derivative:
  # Get the derivatives
  if in_folders_orgnaized:
    filename_ion = str(dataDir+'misc/'+simName+'-ion__source_integrated_moms.gkyl')
  else:
    filename_ion = str(dataDir + simName + '-ion_source_integrated_moms.gkyl')
  pgData_ion = pg.GData(filename_ion)
  M_ion = pgData_ion.get_values()
  M0_ion = np.array(M_ion[:,0])
  M1_ion = np.array(M_ion[:,1])
  M2par_ion = np.array(M_ion[:,2])
  M2perp_ion = np.array(M_ion[:,3])
  time = np.squeeze(np.array(pgData_ion.get_grid()))

  dM0_ion = np.gradient(M0_ion, time)
  dM1_ion = np.gradient(M1_ion, time)
  dM2par_ion = np.gradient(M2par_ion, time)
  dM2perp_ion = np.gradient(M2perp_ion, time)

  fig, ax = plt.subplots(2,2, figsize=(12,12))

  ax[0,0].plot(time, dM0_ion)
  ax[0,0].set_title('dM0/dt')
  ax[0,0].set_xlabel('Time')
  ax[0,0].set_ylabel('Derivative of integrated M0')

  ax[0,1].plot(time, dM1_ion)
  ax[0,1].set_title('dM1/dt')
  ax[0,1].set_xlabel('Time')
  ax[0,1].set_ylabel('Derivative of integrated M1')

  ax[1,0].plot(time, dM2par_ion)
  ax[1,0].set_title('dM2par/dt')
  ax[1,0].set_xlabel('Time')
  ax[1,0].set_ylabel('Derivative of integrated M2 parallel')

  ax[1,1].plot(time, dM2perp_ion)
  ax[1,1].set_title('dM2perp/dt')
  ax[1,1].set_xlabel('Time')
  ax[1,1].set_ylabel('Derivative of integrated M2 perpendicular')

  plt.savefig(outDir + 'source_integrated_moms_derivatives.png')
  plt.show()

if plot_steady_state_comparison:
  # Get the latest BiMaxwellian Moments

  if in_folders_orgnaized:
    filename_ion_latest = str(dataDir+'BiMaxwellianMoments/'+simName+'-ion_BiMaxwellianMoments_'+str(frame_max_plus1-1)+'.gkyl')
  else:
    filename_ion_latest = str(dataDir+simName+'-ion_BiMaxwellianMoments_'+str(frame_max_plus1-1)+'.gkyl')
  pgData_ion = pg.GData(filename_ion_latest)
  pgInterp_ion = pg.GInterpModal(pgData_ion, polyOrder, 'ms')
  _, n = pgInterp_ion.interpolate(0)
  _, upar = pgInterp_ion.interpolate(1)
  _, Tpar = pgInterp_ion.interpolate(2)
  _, Tperp = pgInterp_ion.interpolate(3)

  if in_folders_orgnaized:
    filename_ion_1ms_ago = str(dataDir+'BiMaxwellianMoments/'+simName+'-ion_BiMaxwellianMoments_'+str(frame_max_plus1-101)+'.gkyl')
  else:
    filename_ion_1ms_ago = str(dataDir+simName+'-ion_BiMaxwellianMoments_'+str(frame_max_plus1-101)+'.gkyl')
  pgData_ion_1ms_ago = pg.GData(filename_ion_1ms_ago)
  pgInterp_ion_1ms_ago = pg.GInterpModal(pgData_ion_1ms_ago, polyOrder, 'ms')
  _, n_1ms_ago = pgInterp_ion_1ms_ago.interpolate(0)
  _, upar_1ms_ago = pgInterp_ion_1ms_ago.interpolate(1)
  _, Tpar_1ms_ago = pgInterp_ion_1ms_ago.interpolate(2)
  _, Tperp_1ms_ago = pgInterp_ion_1ms_ago.interpolate(3)

  n = n[:,0]
  upar = upar[:,0]
  Tpar = Tpar[:,0] * mi / eV
  Tperp = Tperp[:,0] * mi / eV

  n_1ms_ago = n_1ms_ago[:,0]
  upar_1ms_ago = upar_1ms_ago[:,0]
  Tpar_1ms_ago = Tpar_1ms_ago[:,0] * mi / eV
  Tperp_1ms_ago = Tperp_1ms_ago[:,0] * mi / eV

  if in_folders_orgnaized:
    filename_mc2nu_pos = pg.GData(str(dataDir+'Geometry/'+simName+'-mapc2p.gkyl'))
  else:
    filename_mc2nu_pos = pg.GData(str(dataDir+simName+'-mapc2p.gkyl'))
  interp = pg.GInterpModal(filename_mc2nu_pos, 1, 'ms')
  nodes_Z = interp.interpolate(1)[1]
  nodes_Z = np.squeeze(nodes_Z)

  fig, ax = plt.subplots(2, 2, figsize=(8,6))
  fig.suptitle('Steady state comparison')

  ax[0,0].plot(nodes_Z, n, label='Latest')
  ax[0,0].plot(nodes_Z, n_1ms_ago, '--', label='1ms ago')
  ax[0,0].set_xlabel('Field line length, radians')
  ax[0,0].set_ylabel('$n_i$, $m^{-3}$')
  # ax[0,0].legend()

  ax[0,1].plot(nodes_Z, upar, label='Latest')
  ax[0,1].plot(nodes_Z, upar_1ms_ago, '--', label='1ms ago')
  ax[0,1].set_xlabel('Field line length, radians')
  ax[0,1].set_ylabel('$U_{||}$, $m/s$')
  ax[0,1].legend()

  ax[1,0].plot(nodes_Z, Tpar, label='Latest')
  ax[1,0].plot(nodes_Z, Tpar_1ms_ago, '--', label='1ms ago')
  ax[1,0].set_xlabel('Field line length, radians')
  ax[1,0].set_ylabel('$T_{||}$, $eV$')
  # ax[1,0].legend()

  ax[1,1].plot(nodes_Z, Tperp, label='Latest')
  ax[1,1].plot(nodes_Z, Tperp_1ms_ago, '--', label='1ms ago')
  ax[1,1].set_xlabel('Field line length, radians')
  ax[1,1].set_ylabel('$T_{\perp}$, $eV$')
  # ax[1,1].legend()
  plt.tight_layout()

  plt.savefig(outDir + 'steady_state_comparison.png')
  plt.show()

  fig, ax = plt.subplots(2, 2, figsize=(8,6))
  fig.suptitle('Steady state difference')

  ax[0,0].plot(nodes_Z, (n - n_1ms_ago)/n_1ms_ago)
  ax[0,0].set_xlabel('Field line length, radians')
  ax[0,0].set_ylabel('$\Delta n / n$, $m^{-3}$')

  ax[0,1].plot(nodes_Z, (upar - upar_1ms_ago)/upar_1ms_ago)
  ax[0,1].set_xlabel('Field line length, radians')
  ax[0,1].set_ylabel('$\Delta U_{||}/ U_{||}$, $m/s$')

  ax[1,0].plot(nodes_Z, (Tpar - Tpar_1ms_ago)/Tpar_1ms_ago)
  ax[1,0].set_xlabel('Field line length, radians')
  ax[1,0].set_ylabel('$\Delta T_{||} / T_{||}$, $eV$')

  ax[1,1].plot(nodes_Z, (Tperp - Tperp_1ms_ago)/Tperp_1ms_ago)
  ax[1,1].set_xlabel('Field line length, radians')
  ax[1,1].set_ylabel('$\Delta T_{\perp}, / T_{\perp}$, $eV$')
  plt.tight_layout()
  plt.savefig(outDir + 'steady_state_difference.png')
  plt.show()







