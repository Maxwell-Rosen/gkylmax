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


# dataDir = '/home/mr1884/scratch/Link to scratch_traverse/gkylmax/traverse-wham1x-compare_unif_vs_nonunif/outputs/'
# dataDir = './data-hires-lorad/'
dataDir = './'
simName = 'gk_wham'
frame_max_plus1 = 33
time_per_frame = 10e-6

plot_potential_trace = 0
plot_bimax_moms = 0
plot_bimax_moms_2D_time_trace = 0
plot_intEnergy_trace = 1
plot_integrate_positivity = 0

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
  filename_bmag = str(dataDir+'Geometry/'+simName+'-bmag.gkyl')
  pgData_bmag = pg.GData(filename_bmag)
  pgInterp_bmag = pg.GInterpModal(pgData_bmag, polyOrder, 'ms')
  x_bmag, dataOut_bmag = pgInterp_bmag.interpolate()

  bmag_shape = dataOut_bmag.shape
  midpoint = int(bmag_shape[0]/2)
  upperhalf = dataOut_bmag[midpoint:]
  peak = np.argmax(upperhalf)
  peak_idx = midpoint+peak

  def loadphi(frame_number, filename):
    filename_phi = str(dataDir+'Field/'+filename+'-field_'+str(frame_number)+'.gkyl')
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

    filename_ion = str(dataDir+'BiMaxwellianMoments/'+simName+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
    pgData_ion = pg.GData(filename_ion)
    pgInterp_ion = pg.GInterpModal(pgData_ion, polyOrder, 'ms')
    coords, n_ion = pgInterp_ion.interpolate(0)
    coords, u_ion = pgInterp_ion.interpolate(1)
    coords, Tpar_ion = pgInterp_ion.interpolate(2)
    coords, Tperp_ion = pgInterp_ion.interpolate(3)

    filename_field = str(dataDir+'Field/'+simName+'-field_'+str(frame_number)+'.gkyl')
    pgData_field = pg.GData(filename_field)
    pgInterp_field = pg.GInterpModal(pgData_field, polyOrder, 'ms')
    coords, phi = pgInterp_field.interpolate()
    
    filename_mc2nu_pos = pg.GData(str(dataDir+'Geometry/'+simName+'-mc2nu_pos.gkyl'))
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
      filename_ion = str(dataDir+'BiMaxwellianMoments/'+simName+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
      pgData_ion = pg.GData(filename_ion)
      pgInterp_ion = pg.GInterpModal(pgData_ion, polyOrder, 'ms')
      _, n = pgInterp_ion.interpolate(0)
      _, upar = pgInterp_ion.interpolate(1)
      _, Tpar = pgInterp_ion.interpolate(2)
      _, Tperp = pgInterp_ion.interpolate(3)

      filename_field = str(dataDir+'Field/'+simName+'-field_'+str(frame_number)+'.gkyl')
      pgData_field = pg.GData(filename_field)
      pgInterp_field = pg.GInterpModal(pgData_field, polyOrder, 'ms')
      _, phi = pgInterp_field.interpolate()

      n = n[:,0]
      upar = upar[:,0]
      Tpar = Tpar[:,0]
      Tperp = Tperp[:,0]
      phi = phi[:,0]
      return n, upar, Tpar, Tperp, phi
    
    filename_mc2nu_pos = pg.GData(str(dataDir+'Geometry/'+simName+'-mc2nu_pos.gkyl'))
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
      # filename_ion = str(dataDir+'BiMaxwellianMoments/'+simName+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
      filename_ion = str(dataDir+simName+'-ion_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
      pgData_ion = pg.GData(filename_ion)
      pgInterp_ion = pg.GInterpModal(pgData_ion, polyOrder, 'ms')
      _, n = pgInterp_ion.interpolate(0)
      _, upar = pgInterp_ion.interpolate(1)
      _, Tpar = pgInterp_ion.interpolate(2)
      _, Tperp = pgInterp_ion.interpolate(3)

      # filename_field = str(dataDir+'Field/'+simName+'-field_'+str(frame_number)+'.gkyl')
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
    
    # filename_mc2nu_pos = pg.GData(str(dataDir+'Geometry/'+simName+'-mc2nu_pos.gkyl'))
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

    # filename_ion = str(dataDir+'misc/'+simName+'-ion_integrated_moms.gkyl')
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

    # filename_ion_positivity = str(dataDir+'misc/'+simName+'-ion_positivity_shift_integrated_moms.gkyl')
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


#   #....................................DEPRICATED CODE............................................#

#     densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_unif = pg.GData(densityFileName_unif)
#     pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'ms')
#     x_unif, dataOut_unif = pgInterp_unif.interpolate()


#     nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
#     pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
#     pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
#     x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pgInterp_nonunif_mapc2p.interpolate(2)


#   if plot_moments:
#     dataOut_unif, dataOut_unif_mapc2p, dataOut_nonunif, dataOut_nonunif_mapc2p = load_mapped_data('-ion_M0_')
    
#     # Create a subfigure that is 2 by 3
#     fig, ax = plt.subplots(2, 3, figsize=(20,10))
#     # Plot the density
#     ax[0,0].plot(dataOut_unif_mapc2p[:,0], dataOut_unif[:,0],'r', label=unifFile)
#     ax[0,0].plot(dataOut_nonunif_mapc2p[:,0], dataOut_nonunif[:,0],'b--', label=nonunifFile)
#     ax[0,0].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
#     ax[0,0].set_ylabel('$n_i$ (m$^{-3}$)', fontsize=xyLabelFontSize)
#     ax[0,0].legend(loc='upper left', fontsize=legendFontSize)
#     setTickFontSize(ax[0,0],tickFontSize)

#     # Plot phi
#     dataOut_unif, dataOut_unif_mapc2p, dataOut_nonunif, dataOut_nonunif_mapc2p = load_mapped_data('-field_')

#     dataOut_unif *= eV/Te0
#     dataOut_nonunif *= eV/Te0
#     # dataOut_reduced *= eV/Te0
#     # dataOut_coarse *= eV/Te0

#     ax[0,1].plot(dataOut_unif_mapc2p[:,0], dataOut_unif[:,0],'r')
#     ax[0,1].plot(dataOut_nonunif_mapc2p[:,0], dataOut_nonunif[:,0],'b--')
#     ax[0,1].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
#     ax[0,1].set_ylabel('$\phi $(V)$)', fontsize=xyLabelFontSize)
#     # ax[0,1].set_title('Frame '+str(frameNum))
#     ax[0,1].set_title("time = "+str(frameNum*1e-5)+" s")
#     setTickFontSize(ax[0,1],tickFontSize)

#     # Plot uPar
#     M0_unif, M0_map, M0_nonunif, M0_nonunif_map = load_mapped_data('-ion_M0_')
#     M1_unif, M1_map, M1_nonunif, M1_nonunif_map = load_mapped_data('-ion_M1_')

#     upar_unif = M1_unif[:,0]/M0_unif[:,0]
#     upar_nonunif = M1_nonunif[:,0]/M0_nonunif[:,0]

#     ax[1,0].plot(M0_map[:,0], upar_unif / c_s,'r')
#     ax[1,0].plot(M0_nonunif_map[:,0], upar_nonunif / c_s,'b--')
#     ax[1,0].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
#     ax[1,0].set_ylabel('$u_{\parallel} / c_s$ (m/s)', fontsize=xyLabelFontSize)
#     setTickFontSize(ax[1,0],tickFontSize)
    
#     # Plot tPerp
#     # M0_unif, M0_map, M0_nonunif, M0_nonunif_map, M0_reduced, M0_reduced_map = load_mapped_data('-ion_M0_')
#     M2perp_unif, M2perp_map, M2perp_nonunif, M2perp_nonunif_map = load_mapped_data('-ion_M2perp_')

#     tPerp_unif = M2perp_unif[:,0]/M0_unif[:,0] * mi / eV
#     tPerp_nonunif = M2perp_nonunif[:,0]/M0_nonunif[:,0] * mi / eV

#     ax[1,1].plot(M0_map[:,0], tPerp_unif,'r')
#     ax[1,1].plot(M0_nonunif_map[:,0], tPerp_nonunif,'b--')
#     ax[1,1].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
#     ax[1,1].set_ylabel('$T_{\perp}$ (eV)', fontsize=xyLabelFontSize)
#     setTickFontSize(ax[1,1],tickFontSize)

#     # Plot tPar
#     # M0_unif, M0_map, M0_nonunif, M0_nonunif_map, M0_reduced, M0_reduced_map = load_mapped_data('-ion_M0_')
#     # M1_unif, M1_map, M1_nonunif, M1_nonunif_map, M1_reduced, M1_reduced_map = load_mapped_data('-ion_M1_')
#     M2par_unif, M2par_map, M2par_nonunif, M2par_nonunif_map = load_mapped_data('-ion_M2par_')

#     tPar_unif = (M2par_unif[:,0] - M1_unif[:,0]**2/M0_unif[:,0]) * mi / eV / M0_unif[:,0]
#     tPar_nonunif = (M2par_nonunif[:,0] - M1_nonunif[:,0]**2/M0_nonunif[:,0]) * mi / eV / M0_nonunif[:,0]

#     ax[1,2].plot(M0_map[:,0], tPar_unif,'r')
#     ax[1,2].plot(M0_nonunif_map[:,0], tPar_nonunif,'b--')
#     ax[1,2].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
#     ax[1,2].set_ylabel('$T_{\parallel}$ (eV)', fontsize=xyLabelFontSize)
#     setTickFontSize(ax[1,2],tickFontSize)

#     # Plot the grid and mapc2p
#     # ax[0,2].plot(dataOut_unif_mapc2p[:,0], dataOut_unif_mapc2p[:,0],'r', label='Uniform 280x96x192', markersize=0.5)
#     # ax[0,2].plot(dataOut_nonunif_mapc2p[:,0], dataOut_nonunif_mapc2p[:,0],'b.', label='Nonuniform 280x96x192', markersize=0.5)
#     # ax[0,2].plot(dataOut_reduced_mapc2p[:,0], dataOut_reduced_mapc2p[:,0],'g.', label='Nonuniform 140x96x192', markersize=0.5)
#     # ax[0,2].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
#     # ax[0,2].set_ylabel('Mapped cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
#     # ax[0,2].legend(loc='upper left', fontsize=legendFontSize)
#     # setTickFontSize(ax[0,2],tickFontSize)
    
#     # print("mapc2p uniform grid: ", dataOut_unif_mapc2p[:,0])
#     # print("mapc2p nonuniform grid: ", dataOut_nonunif_mapc2p[:,0])
#     # print("last cell spacings uniform grid: ", dataOut_unif_mapc2p[-1,0] - dataOut_unif_mapc2p[-2,0])
#     # print("last cell spacings nonuniform grid: ", dataOut_nonunif_mapc2p[-1,0] - dataOut_nonunif_mapc2p[-2,0])
#     # print("ratio of last cell spacings: ", (dataOut_unif_mapc2p[-1,0] - dataOut_unif_mapc2p[-2,0])/(dataOut_nonunif_mapc2p[-1,0] - dataOut_nonunif_mapc2p[-2,0]))
    
#     figName = 'moments_'+str(frameNum)
#     if save_figure_as_file:
#       fig.savefig(outDir+figName+figureFileFormat)
#       plt.close()
#       print('Saved figure as '+figName+figureFileFormat)  
#     else:
#       plt.show()
  
#   if plot_distvpar:
#     # f_unif, f_map, f_nonunif, f_nonunif_map = load_mapped_data('-ion-')
#     dataName = '-ion_'
#     densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_unif = pg.GData(densityFileName_unif)
#     pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'gkhyb')
#     x_unif, dataOut_unif = pgInterp_unif.interpolate()
#     dataOut_unif = np.squeeze(dataOut_unif)

#     densityFileName_nonunif = str(dataDir+nonunifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_nonunif = pg.GData(densityFileName_nonunif)
#     pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'gkhyb')
#     x_nonunif, dataOut_nonunif = pgInterp_nonunif.interpolate()
#     dataOut_nonunif = np.squeeze(dataOut_nonunif)

#     nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
#     pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
#     pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
#     x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pgInterp_nonunif_mapc2p.interpolate(2)
#     dataOut_nonunif_mapc2p = np.squeeze(dataOut_nonunif_mapc2p)

#     unif_mapc2p_filename = str(dataDir+unifFile+'-mapc2p.gkyl')
#     pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
#     pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
#     x_unif_mapc2p, dataOut_unif_mapc2p = pgInterp_unif_mapc2p.interpolate(2)
#     dataOut_unif_mapc2p = np.squeeze(dataOut_unif_mapc2p)

#     # Convert from cell center to edges
#     zmin = x_unif_mapc2p[0][0]
#     zmax = x_unif_mapc2p[0][-1]
#     diffs  = dataOut_nonunif_mapc2p[0:-1] + np.diff(dataOut_nonunif_mapc2p)/2
#     edged_dataOut_nonunif_mapc2p = np.insert(diffs, 0, zmin)
#     edged_dataOut_nonunif_mapc2p = np.append(edged_dataOut_nonunif_mapc2p, zmax)

#     # Convert from cell center to edges
#     diffs  = dataOut_unif_mapc2p[0:-1] + np.diff(dataOut_unif_mapc2p)/2
#     edged_dataOut_unif_mapc2p = np.insert(diffs, 0, zmin)
#     edged_dataOut_unif_mapc2p = np.append(edged_dataOut_unif_mapc2p, zmax)

#     unif_Jgeo_filename = str(dataDir+unifFile+'-jacobgeo.gkyl')
#     pgData_unif_Jgeo = pg.GData(unif_Jgeo_filename)
#     pgInterp_unif_Jgeo = pg.GInterpModal(pgData_unif_Jgeo, polyOrder, 'ms')
#     x_unif_Jgeo, dataOut_unif_Jgeo = pgInterp_unif_Jgeo.interpolate()
#     dataOut_unif_Jgeo = np.squeeze(dataOut_unif_Jgeo)

#     nonunif_Jgeo_filename = str(dataDir+nonunifFile+'-jacobgeo.gkyl')
#     pgData_nonunif_Jgeo = pg.GData(nonunif_Jgeo_filename)
#     pgInterp_nonunif_Jgeo = pg.GInterpModal(pgData_nonunif_Jgeo, polyOrder, 'ms')
#     x_nonunif_Jgeo, dataOut_nonunif_Jgeo = pgInterp_nonunif_Jgeo.interpolate()
#     dataOut_nonunif_Jgeo = np.squeeze(dataOut_nonunif_Jgeo)

#     unif_distf_shape = dataOut_unif.shape
#     nonunif_distf_shape = dataOut_nonunif.shape

#     tile_unif_Jgeo = np.ones((unif_distf_shape[0], unif_distf_shape[1], unif_distf_shape[2]))
#     for i in range(unif_distf_shape[0]):
#         tile_unif_Jgeo[i,:,:] *= dataOut_unif_Jgeo[i]

#     tile_nonunif_Jgeo = np.ones((nonunif_distf_shape[0], nonunif_distf_shape[1], nonunif_distf_shape[2]))
#     for i in range(nonunif_distf_shape[0]):
#         tile_nonunif_Jgeo[i,:,:] *= dataOut_nonunif_Jgeo[i]

#     dataOut_unif = np.trapz(np.abs(dataOut_unif / tile_unif_Jgeo), axis=2)
#     dataOut_nonunif = np.trapz(np.abs(dataOut_nonunif / tile_nonunif_Jgeo), axis=2)

#     # Interpolate the non-uniform data onto a uniform grid
#     dataOut_unif_shape = dataOut_unif.shape
#     dataOut_nonunif_interp = np.zeros((dataOut_unif_shape[0], dataOut_unif_shape[1]))
#     for i in range(nonunif_distf_shape[1]):
#         dataOut_nonunif_interp[:,i] = np.interp(dataOut_unif_mapc2p, dataOut_nonunif_mapc2p, dataOut_nonunif[:,i])
#     data_difference = np.abs((dataOut_unif - dataOut_nonunif_interp))
#   #   # Need velocity space grids or to convert to edges of z for plotting

#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,12))

#     norm = LogNorm(vmin = 1e-14, vmax = 1e-4)  # Create a LogNorm instance

#     pcolormesh1 = ax1.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[1]*np.sqrt(2), dataOut_unif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh1, ax=ax1)  # Add a colorbar to the plot

#     # Label the axes
#     ax1.set_ylabel('vpar / vpar_max')
#     ax1.set_xlabel('Uniform grid Z, cylindrical coodinate (m)')
#     ax1.set_title('Frame '+str(frameNum))

#     # print(x_unif[1]/vti)
#     # print(x_nonunif[1]/vti)

#     # pcolormesh2 = ax2.pcolormesh(edged_dataOut_nonunif_mapc2p, x_nonunif[1]/vti, dataOut_nonunif.T, cmap='inferno', norm=norm)
#     pcolormesh2 = ax2.pcolormesh(edged_dataOut_nonunif_mapc2p, mapc2p_vel_vpar(x_nonunif[1]), dataOut_nonunif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh2, ax=ax2)  # Add a colorbar to the plot
#     # Label the axes
#     ax2.set_ylabel('vpar / vpar_max')
#     ax2.set_xlabel('Nonuniform grid Z, cylindrical coodinate (m)') 
#   #   pcolormesh3 = ax3.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[1] / vti, data_difference.T, cmap='inferno', norm=norm)
#   #   fig.colorbar(pcolormesh3, ax=ax3)  # Add a colorbar to the plot
#   # # Label the axes
#   #   ax3.set_ylabel('vpar / vti')
#   #   ax3.set_xlabel('Absolute difference between distribution functions Z, cylindrical coordinate (m)')
    
#     figName = 'distf_vpar_'+str(frameNum)
#     if save_figure_as_file:
#       fig.savefig(outDir+figName+figureFileFormat)
#       plt.close()
#       print('Saved figure as '+figName+figureFileFormat)  
#     else:
#       plt.show()  

#   if plot_distmu:
#     # f_unif, f_map, f_nonunif, f_nonunif_map = load_mapped_data('-ion-')
#     dataName = '-ion_'
#     densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_unif = pg.GData(densityFileName_unif)
#     pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'ms')
#     x_unif, dataOut_unif = pgInterp_unif.interpolate()
#     dataOut_unif = np.squeeze(dataOut_unif)

#     densityFileName_nonunif = str(dataDir+nonunifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_nonunif = pg.GData(densityFileName_nonunif)
#     pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'ms')
#     x_nonunif, dataOut_nonunif = pgInterp_nonunif.interpolate()
#     dataOut_nonunif = np.squeeze(dataOut_nonunif)

#     nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
#     pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
#     pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
#     x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pgInterp_nonunif_mapc2p.interpolate(2)
#     dataOut_nonunif_mapc2p = np.squeeze(dataOut_nonunif_mapc2p)

#     unif_mapc2p_filename = str(dataDir+unifFile+'-mapc2p.gkyl')
#     pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
#     pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
#     x_unif_mapc2p, dataOut_unif_mapc2p = pgInterp_unif_mapc2p.interpolate(2)
#     dataOut_unif_mapc2p = np.squeeze(dataOut_unif_mapc2p)

#     # Convert from cell center to edges
#     zmin = x_unif_mapc2p[0][0]
#     zmax = x_unif_mapc2p[0][-1]
#     diffs  = dataOut_nonunif_mapc2p[0:-1] + np.diff(dataOut_nonunif_mapc2p)/2
#     edged_dataOut_nonunif_mapc2p = np.insert(diffs, 0, zmin)
#     edged_dataOut_nonunif_mapc2p = np.append(edged_dataOut_nonunif_mapc2p, zmax)

#     # Convert from cell center to edges
#     diffs  = dataOut_unif_mapc2p[0:-1] + np.diff(dataOut_unif_mapc2p)/2
#     edged_dataOut_unif_mapc2p = np.insert(diffs, 0, zmin)
#     edged_dataOut_unif_mapc2p = np.append(edged_dataOut_unif_mapc2p, zmax)

#     unif_Jgeo_filename = str(dataDir+unifFile+'-jacobgeo.gkyl')
#     pgData_unif_Jgeo = pg.GData(unif_Jgeo_filename)
#     pgInterp_unif_Jgeo = pg.GInterpModal(pgData_unif_Jgeo, polyOrder, 'ms')
#     x_unif_Jgeo, dataOut_unif_Jgeo = pgInterp_unif_Jgeo.interpolate()
#     dataOut_unif_Jgeo = np.squeeze(dataOut_unif_Jgeo)

#     nonunif_Jgeo_filename = str(dataDir+nonunifFile+'-jacobgeo.gkyl')
#     pgData_nonunif_Jgeo = pg.GData(nonunif_Jgeo_filename)
#     pgInterp_nonunif_Jgeo = pg.GInterpModal(pgData_nonunif_Jgeo, polyOrder, 'ms')
#     x_nonunif_Jgeo, dataOut_nonunif_Jgeo = pgInterp_nonunif_Jgeo.interpolate()
#     dataOut_nonunif_Jgeo = np.squeeze(dataOut_nonunif_Jgeo)

#     unif_distf_shape = dataOut_unif.shape
#     nonunif_distf_shape = dataOut_nonunif.shape

#     tile_unif_Jgeo = np.ones((unif_distf_shape[0], unif_distf_shape[1], unif_distf_shape[2]))
#     for i in range(unif_distf_shape[0]):
#         tile_unif_Jgeo[i,:,:] *= dataOut_unif_Jgeo[i]

#     tile_nonunif_Jgeo = np.ones((nonunif_distf_shape[0], nonunif_distf_shape[1], nonunif_distf_shape[2]))
#     for i in range(nonunif_distf_shape[0]):
#         tile_nonunif_Jgeo[i,:,:] *= dataOut_nonunif_Jgeo[i]

#     dataOut_unif = np.trapz(np.abs(dataOut_unif / tile_unif_Jgeo), axis=1)
#     dataOut_nonunif = np.trapz(np.abs(dataOut_nonunif / tile_nonunif_Jgeo), axis=1)

#     # Interpolate the non-uniform data onto a uniform grid
#     dataOut_unif_shape = dataOut_unif.shape
#     dataOut_nonunif_shape = dataOut_nonunif.shape
#     dataOut_nonunif_interp = np.zeros((dataOut_unif_shape[0], dataOut_unif_shape[1]))
#     for i in range(dataOut_nonunif_shape[1]):
#         dataOut_nonunif_interp[:,i] = np.interp(dataOut_unif_mapc2p, dataOut_nonunif_mapc2p, dataOut_nonunif[:,i])
#     data_difference = np.abs((dataOut_unif - dataOut_nonunif_interp))

#     # Need velocity space grids or to convert to edges of z for plotting
#     fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,12))

#     norm = LogNorm(vmin = 1e-15, vmax = 1e-6)  # Create a LogNorm instance

#     pcolormesh1 = ax1.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[2], dataOut_unif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh1, ax=ax1)  # Add a colorbar to the plot
#     # Label the axes
#     ax1.set_ylabel('mu / mu_max')
#     ax1.set_xlabel('Uniform grid Z, cylindrical coodinate (m)')
#     ax1.set_title('Frame '+str(frameNum))

#     pcolormesh2 = ax2.pcolormesh(edged_dataOut_nonunif_mapc2p, mapc2p_vel_mu(x_nonunif[2]), dataOut_nonunif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh2, ax=ax2)  # Add a colorbar to the plot
#     # Label the axes
#     ax2.set_ylabel('mu / mu_max')
#     ax2.set_xlabel('Nonuniform grid Z, cylindrical coodinate (m)')

#   #   pcolormesh3 = ax3.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[2] / vti, data_difference.T, cmap='inferno', norm=norm)
#   #   fig.colorbar(pcolormesh3, ax=ax3)  # Add a colorbar to the plot
#   # # Label the axes
#   #   ax3.set_ylabel('vpar / vti')
#   #   ax3.set_xlabel('Absolute difference between distribution functions Z, cylindrical coordinate (m)')
    
#     figName = 'distf_mu_'+str(frameNum)
#     if save_figure_as_file:
#       fig.savefig(outDir+figName+figureFileFormat)
#       plt.close()
#       print('Saved figure as '+figName+figureFileFormat)  
#     else:
#       plt.show() 
      
#   def plot_distf_at_z_eq(z0_coordinate):
#     dataName = '-ion_'
#     densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_unif = pg.GData(densityFileName_unif)
#     pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'ms')
#     pgInterp_unif.interpolate(overwrite=True)
#     x_unif, dataOut_unif = pg.data.select(pgData_unif, z0 = z_xi(z0_coordinate))
#     dataOut_unif = np.squeeze(dataOut_unif)

#     densityFileName_nonunif = str(dataDir+nonunifFile + str(dataName) + str(frameNum) + '.gkyl')
#     pgData_nonunif = pg.GData(densityFileName_nonunif)
#     pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'ms')
#     pgInterp_nonunif.interpolate(overwrite=True)
#     x_nonunif, dataOut_nonunif = pg.data.select(pgData_nonunif, z0 = (z0_coordinate))
#     dataOut_nonunif = np.squeeze(dataOut_nonunif)

#     nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
#     pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
#     pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
#     pgInterp_nonunif_mapc2p.interpolate(overwrite=True)
#     x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pg.data.select(pgData_nonunif_mapc2p, z0 = (z0_coordinate))
#     dataOut_nonunif_mapc2p = np.squeeze(dataOut_nonunif_mapc2p)

#     unif_mapc2p_filename = str(dataDir+unifFile+'-mapc2p.gkyl')
#     pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
#     pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
#     pgInterp_unif_mapc2p.interpolate(overwrite=True)
#     x_unif_mapc2p, dataOut_unif_mapc2p = pg.data.select(pgData_unif_mapc2p, z0 = z_xi(z0_coordinate))
#     dataOut_unif_mapc2p = np.squeeze(dataOut_unif_mapc2p)

#     unif_Jgeo_filename = str(dataDir+unifFile+'-jacobgeo.gkyl')
#     pgData_unif_Jgeo = pg.GData(unif_Jgeo_filename)
#     pgInterp_unif_Jgeo = pg.GInterpModal(pgData_unif_Jgeo, polyOrder, 'ms')
#     pgInterp_unif_Jgeo.interpolate(overwrite=True)
#     x_unif_Jgeo, dataOut_unif_Jgeo = pg.data.select(pgData_unif_Jgeo, z0 = z_xi(z0_coordinate))
#     dataOut_unif_Jgeo = np.squeeze(dataOut_unif_Jgeo)

#     nonunif_Jgeo_filename = str(dataDir+nonunifFile+'-jacobgeo.gkyl')
#     pgData_nonunif_Jgeo = pg.GData(nonunif_Jgeo_filename)
#     pgInterp_nonunif_Jgeo = pg.GInterpModal(pgData_nonunif_Jgeo, polyOrder, 'ms')
#     pgInterp_nonunif_Jgeo.interpolate(overwrite=True)
#     x_nonunif_Jgeo, dataOut_nonunif_Jgeo = pg.data.select(pgData_nonunif_Jgeo, z0 = (z0_coordinate))
#     dataOut_nonunif_Jgeo = np.squeeze(dataOut_nonunif_Jgeo)

#     dataOut_unif = np.abs(dataOut_unif / dataOut_unif_Jgeo)
#     dataOut_nonunif = np.abs(dataOut_nonunif / dataOut_nonunif_Jgeo)
#     dataOut_diff = np.abs((dataOut_unif - dataOut_nonunif))

#     # Need velocity space grids or to convert to edges of z for plotting
#     fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18,6))

#     norm = LogNorm(vmin = 1e-6, vmax = np.max(dataOut_unif))  # Create a LogNorm instance

#     pcolormesh1 = ax1.pcolormesh(x_unif[1]/vti, x_unif[2]/mui0, dataOut_unif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh1, ax=ax1)  # Add a colorbar to the plot
#     # Label the axes
#     ax1.set_ylabel('mu / mui0')
#     ax1.set_xlabel('vpar / vti')
#     ax1.set_title('Frame '+str(frameNum))

#     pcolormesh2 = ax2.pcolormesh(x_nonunif[1]/vti, x_nonunif[2]/mui0, dataOut_nonunif.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh2, ax=ax2)  # Add a colorbar to the plot
#     # Label the axes
#     ax2.set_ylabel('mu / mui0')
#     ax2.set_xlabel('vpar / vti')
    
#     pcolormesh3 = ax3.pcolormesh(x_unif[1]/vti, x_unif[2] / vti, dataOut_diff.T, cmap='inferno', norm=norm)
#     fig.colorbar(pcolormesh3, ax=ax3)  # Add a colorbar to the plot
#     # Label the axes
#     ax3.set_ylabel('vpar / vti')
#     ax3.set_xlabel('Absolute difference between distribution functions Z, cylindrical coordinate (m)')

#     ax2.title.set_text('z = '+str(z_xi(z0_coordinate))+' m')

#     locName = np.round(z0_coordinate*100)
#     figName = 'distf_z'+str(locName)+'_'+str(frameNum)
#     if save_figure_as_file:
#       fig.savefig(outDir+figName+figureFileFormat)
#       plt.close()
#       print('Saved figure as '+figName+figureFileFormat)  
#     else:
#       plt.show() 

#     # On axis difference
#     reduce_set_unif = dataOut_unif[:,0]
#     reduce_set_nonunif = dataOut_nonunif[:,0]
#     x_axis = x_unif[1]/vti

#     # I need to interpolate for x 
#     plt.plot(x_axis[1:], reduce_set_unif, label='Uniform grid')
#     plt.plot(x_axis[1:], reduce_set_nonunif, label='Nonuniform grid')
#     plt.plot(x_axis[1:], np.abs(reduce_set_unif - reduce_set_nonunif), 'k--', label='Absolute difference')
#     plt.xlabel('vpar / vti')
#     plt.ylim(np.max(reduce_set_nonunif)* 1e-6, np.max(reduce_set_nonunif)*1.5)
#     plt.xlim(-2, 2)

#     plt.ylabel('f')
#     # log scale y axis
#     plt.yscale('log')
#     plt.legend()
#     plt.title('On axis difference, mu=0, frame '+str(frameNum)+ ', z = '+str(z_xi(z0_coordinate))+' m')

#     figName = 'distf_z'+str(locName)+'mu0'+str(frameNum)
#     if save_figure_as_file:
#       plt.savefig(outDir+figName+figureFileFormat)
#       plt.close()
#       print('Saved figure as '+figName+figureFileFormat)
#     else:
#       plt.show()

#   if plot_distf_at_z:
#     for z0 in z_loctions:
#       plot_distf_at_z_eq(z0)


# # # process_frame(0)
# # # Number of processes to run in parallel
# # num_processes = multiprocessing.cpu_count()
# # print('Number of processes: ', num_processes)

# # # Create a pool of processes
# # pool = multiprocessing.Pool(processes=num_processes)

# # # Map the frame_arr to the pool of processes
# # pool.map(process_frame, frame_arr)

# # # Close the pool to prevent any more tasks from being submitted
# # pool.close()

# # # Wait for all processes to finish
# # pool.join()
  