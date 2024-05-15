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

# dataDir = '/home/mr1884/scratch/Link to scratch_traverse/gkylmax/traverse-wham1x-compare_unif_vs_nonunif/outputs/'
dataDir = '/home/mr1884/scratch/Link to scratch_traverse/gkylmax/traverse-wham1x-compare_unif_vs_nonunif/'
unifFile = 'gk_wham_unif'
nonunifFile = 'gk_wham_nonunif'
reducedFile = 'gk_wham_reduced'
frame_arr = np.array([0])#np.arange(0,1)
# frame_arr = np.array([1:4])
save_figure_as_file= 0     #[ If True, save figure to file. If False, display figure on screen.

plot_moments       = 1  # Plot density, potential, upar, tperp, tpar.
plot_distvpar      = 0  # plot distribution function in vpar.
plot_distmu        = 0  # plot distribution function in mu.
plot_distf_at_z    = 0
z_loctions = np.array([0, 0.3, 0.98, 2.4])

print(frame_arr)
def process_frame(frameNum):
  outDir = './python_plots/'

  figureFileFormat = '.png'    #[ Can be .png, .pdf, .ps, .eps, .svg.

  #[ ............... End of user inputs (MAYBE) ..................... ]#

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

  #[ Thermal speeds.
  vti = np.sqrt(Ti0/mi)
  vte = np.sqrt(Te0/me)
  mui0 = 0.5*mi*(vti**2)/B_p
  c_s = np.sqrt(Te0/mi)

  #[ Gyrofrequencies and gyroradii.
  omega_ci = eV*B_p/mi
  rho_s    = c_s/omega_ci

  z_m = 0.983244   #[ Location of maximum B-field.

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

  def load_mapped_data(dataName):
    densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
    pgData_unif = pg.GData(densityFileName_unif)
    pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'ms')
    x_unif, dataOut_unif = pgInterp_unif.interpolate()

    densityFileName_nonunif = str(dataDir+nonunifFile + str(dataName) + str(frameNum) + '.gkyl')
    pgData_nonunif = pg.GData(densityFileName_nonunif)
    pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'ms')
    x_nonunif, dataOut_nonunif = pgInterp_nonunif.interpolate()
    
    densityFileName_reduced = str(dataDir+reducedFile + str(dataName) + str(frameNum) + '.gkyl')
    pgData_reduced = pg.GData(densityFileName_reduced)
    pgInterp_reduced = pg.GInterpModal(pgData_reduced, polyOrder, 'ms')
    x_reduced, dataOut_reduced = pgInterp_reduced.interpolate()

    nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
    pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
    pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
    x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pgInterp_nonunif_mapc2p.interpolate(2)

    unif_mapc2p_filename = str(dataDir+unifFile+'-mapc2p.gkyl')
    pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
    pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
    x_unif_mapc2p, dataOut_unif_mapc2p = pgInterp_unif_mapc2p.interpolate(2)

    reduced_mapc2p_filename = str(dataDir+reducedFile+'-mapc2p.gkyl')
    pgData_reduced_mapc2p = pg.GData(reduced_mapc2p_filename)
    pgInterp_reduced_mapc2p = pg.GInterpModal(pgData_reduced_mapc2p, polyOrder, 'ms')
    x_reduced_mapc2p, dataOut_reduced_mapc2p = pgInterp_reduced_mapc2p.interpolate(2)

    return dataOut_unif, dataOut_unif_mapc2p, dataOut_nonunif, dataOut_nonunif_mapc2p, dataOut_reduced, dataOut_reduced_mapc2p
  
  def z_xi(xi):
    z_m = 0.983244
    expander_order = 7
    mapping_fraction = 0.7
    center_order = 3
    z_min = -2.5
    z_max = 2.5
    if xi >= z_min and xi <= z_max:
      if xi <= -z_m:
        left = -z_m
        right = z_min
        n = expander_order
      elif xi <= 0.0:
        left = -z_m
        right = 0.0
        n = center_order
      elif xi <= z_m:
        left = z_m
        right = 0.0
        n = center_order
      else:
        left = z_m
        right = z_max
        n = expander_order
      z = (pow(right - left, 1 - n) * pow(xi - left, n) + left) * mapping_fraction + xi * (1 - mapping_fraction)
    else:
      z = xi
    return z


  #................................................................................#

  if plot_moments:
    dataOut_unif, dataOut_unif_mapc2p, dataOut_nonunif, dataOut_nonunif_mapc2p, dataOut_reduced, dataOut_reduced_mapc2p = load_mapped_data('-ion_M0_')
    
    # Create a subfigure that is 2 by 3
    fig, ax = plt.subplots(2, 3, figsize=(20,10))
    # Plot the density
    ax[0,0].plot(dataOut_unif_mapc2p[:,0], dataOut_unif[:,0],'r', label='Uniform grid')
    ax[0,0].plot(dataOut_nonunif_mapc2p[:,0], dataOut_nonunif[:,0],'b--', label='Nonuniform grid')
    ax[0,0].plot(dataOut_reduced_mapc2p[:,0], dataOut_reduced[:,0],'g:', label='Reduced grid')
    ax[0,0].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
    ax[0,0].set_ylabel('$n_i$ (m$^{-3}$)', fontsize=xyLabelFontSize)
    ax[0,0].legend(loc='upper left', fontsize=legendFontSize)
    setTickFontSize(ax[0,0],tickFontSize)

    # Plot phi
    dataOut_unif, dataOut_unif_mapc2p, dataOut_nonunif, dataOut_nonunif_mapc2p, dataOut_reduced, dataOut_reduced_mapc2p = load_mapped_data('-field_')

    dataOut_unif *= eV/Te0
    dataOut_nonunif *= eV/Te0
    dataOut_reduced *= eV/Te0

    ax[0,1].plot(dataOut_unif_mapc2p[:,0], dataOut_unif[:,0],'r', label='Uniform grid')
    ax[0,1].plot(dataOut_nonunif_mapc2p[:,0], dataOut_nonunif[:,0],'b--', label='Nonuniform grid')
    ax[0,1].plot(dataOut_reduced_mapc2p[:,0], dataOut_reduced[:,0],'g:', label='Reduced grid')
    ax[0,1].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
    ax[0,1].set_ylabel('$\phi$ (m$^{-3}$)', fontsize=xyLabelFontSize)
    ax[0,1].legend(loc='upper left', fontsize=legendFontSize)
    ax[0,1].set_title('Frame '+str(frameNum))
    setTickFontSize(ax[0,1],tickFontSize)

    # Plot uPar
    M0_unif, M0_map, M0_nonunif, M0_nonunif_map, M0_reduced, M0_reduced_map = load_mapped_data('-ion_M0_')
    M1_unif, M1_map, M1_nonunif, M1_nonunif_map, M1_reduced, M1_reduced_map = load_mapped_data('-ion_M1_')

    upar_unif = M1_unif[:,0]/M0_unif[:,0]
    upar_nonunif = M1_nonunif[:,0]/M0_nonunif[:,0]
    upar_reduced = M1_reduced[:,0]/M0_reduced[:,0]

    ax[1,0].plot(M0_map[:,0], upar_unif / c_s,'r', label='Uniform grid')
    ax[1,0].plot(M0_nonunif_map[:,0], upar_nonunif / c_s,'b--', label='Nonuniform grid')
    ax[1,0].plot(M0_reduced_map[:,0], upar_reduced / c_s,'g:', label='Reduced grid')
    ax[1,0].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
    ax[1,0].set_ylabel('$u_{\parallel} / c_s$ (m/s)', fontsize=xyLabelFontSize)
    ax[1,0].legend(loc='upper left', fontsize=legendFontSize)
    setTickFontSize(ax[1,0],tickFontSize)
    
    # Plot tPerp
    M0_unif, M0_map, M0_nonunif, M0_nonunif_map, M0_reduced, M0_reduced_map = load_mapped_data('-ion_M0_')
    M2perp_unif, M2perp_map, M2perp_nonunif, M2perp_nonunif_map, M2perp_reduced, M2perp_reduced_map = load_mapped_data('-ion_M2perp_')

    tPerp_unif = M2perp_unif[:,0]/M0_unif[:,0] * mi / eV
    tPerp_nonunif = M2perp_nonunif[:,0]/M0_nonunif[:,0] * mi / eV
    tPerp_reduced = M2perp_reduced[:,0]/M0_reduced[:,0] * mi / eV

    ax[1,1].plot(M0_map[:,0], tPerp_unif,'r', label='Uniform grid')
    ax[1,1].plot(M0_nonunif_map[:,0], tPerp_nonunif,'b--', label='Nonuniform grid')
    ax[1,1].plot(M0_reduced_map[:,0], tPerp_reduced,'g:', label='Reduced grid')
    ax[1,1].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
    ax[1,1].set_ylabel('$T_{\perp}$ (eV)', fontsize=xyLabelFontSize)
    ax[1,1].legend(loc='upper left', fontsize=legendFontSize)
    setTickFontSize(ax[1,1],tickFontSize)

    # Plot tPar
    M0_unif, M0_map, M0_nonunif, M0_nonunif_map, M0_reduced, M0_reduced_map = load_mapped_data('-ion_M0_')
    M1_unif, M1_map, M1_nonunif, M1_nonunif_map, M1_reduced, M1_reduced_map = load_mapped_data('-ion_M1_')
    M2par_unif, M2par_map, M2par_nonunif, M2par_nonunif_map, M2par_reduced, M2par_reduced_map = load_mapped_data('-ion_M2par_')

    tPar_unif = (M2par_unif[:,0] - M1_unif[:,0]**2/M0_unif[:,0]) * mi / eV / M0_unif[:,0]
    tPar_nonunif = (M2par_nonunif[:,0] - M1_nonunif[:,0]**2/M0_nonunif[:,0]) * mi / eV / M0_nonunif[:,0]
    tPar_reduced = (M2par_reduced[:,0] - M1_reduced[:,0]**2/M0_reduced[:,0]) * mi / eV / M0_reduced[:,0]

    ax[1,2].plot(M0_map[:,0], tPar_unif,'r', label='Uniform grid')
    ax[1,2].plot(M0_nonunif_map[:,0], tPar_nonunif,'b--', label='Nonuniform grid')
    ax[1,2].plot(M0_reduced_map[:,0], tPar_reduced,'g:', label='Reduced grid')
    ax[1,2].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
    ax[1,2].set_ylabel('$T_{\parallel}$ (eV)', fontsize=xyLabelFontSize)
    ax[1,2].legend(loc='upper left', fontsize=legendFontSize)
    setTickFontSize(ax[1,2],tickFontSize)

    # Plot the grid and mapc2p
    ax[0,2].plot(dataOut_unif_mapc2p[:,0], dataOut_unif_mapc2p[:,0],'r', label='Uniform grid', markersize=0.5)
    ax[0,2].plot(dataOut_nonunif_mapc2p[:,0], dataOut_nonunif_mapc2p[:,0],'b.', label='Nonuniform grid', markersize=0.5)
    ax[0,2].plot(dataOut_reduced_mapc2p[:,0], dataOut_reduced_mapc2p[:,0],'g.', label='Reduced grid', markersize=0.5)
    ax[0,2].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
    ax[0,2].set_ylabel('Mapped cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
    ax[0,2].legend(loc='upper left', fontsize=legendFontSize)
    setTickFontSize(ax[0,2],tickFontSize)
    # print("mapc2p uniform grid: ", dataOut_unif_mapc2p[:,0])
    # print("mapc2p nonuniform grid: ", dataOut_nonunif_mapc2p[:,0])
    # print("last cell spacings uniform grid: ", dataOut_unif_mapc2p[-1,0] - dataOut_unif_mapc2p[-2,0])
    # print("last cell spacings nonuniform grid: ", dataOut_nonunif_mapc2p[-1,0] - dataOut_nonunif_mapc2p[-2,0])
    # print("ratio of last cell spacings: ", (dataOut_unif_mapc2p[-1,0] - dataOut_unif_mapc2p[-2,0])/(dataOut_nonunif_mapc2p[-1,0] - dataOut_nonunif_mapc2p[-2,0]))
    
    figName = 'moments_'+str(frameNum)
    if save_figure_as_file:
      fig.savefig(outDir+figName+figureFileFormat)
      plt.close()
      print('Saved figure as '+figName+figureFileFormat)  
    else:
      plt.show()
  
  if plot_distvpar:
    # f_unif, f_map, f_nonunif, f_nonunif_map = load_mapped_data('-ion-')
    dataName = '-ion_'
    densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
    pgData_unif = pg.GData(densityFileName_unif)
    pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'ms')
    x_unif, dataOut_unif = pgInterp_unif.interpolate()
    dataOut_unif = np.squeeze(dataOut_unif)

    densityFileName_nonunif = str(dataDir+nonunifFile + str(dataName) + str(frameNum) + '.gkyl')
    pgData_nonunif = pg.GData(densityFileName_nonunif)
    pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'ms')
    x_nonunif, dataOut_nonunif = pgInterp_nonunif.interpolate()
    dataOut_nonunif = np.squeeze(dataOut_nonunif)

    nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
    pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
    pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
    x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pgInterp_nonunif_mapc2p.interpolate(2)
    dataOut_nonunif_mapc2p = np.squeeze(dataOut_nonunif_mapc2p)

    unif_mapc2p_filename = str(dataDir+unifFile+'-mapc2p.gkyl')
    pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
    pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
    x_unif_mapc2p, dataOut_unif_mapc2p = pgInterp_unif_mapc2p.interpolate(2)
    dataOut_unif_mapc2p = np.squeeze(dataOut_unif_mapc2p)

    # Convert from cell center to edges
    zmin = x_unif_mapc2p[0][0]
    zmax = x_unif_mapc2p[0][-1]
    diffs  = dataOut_nonunif_mapc2p[0:-1] + np.diff(dataOut_nonunif_mapc2p)/2
    edged_dataOut_nonunif_mapc2p = np.insert(diffs, 0, zmin)
    edged_dataOut_nonunif_mapc2p = np.append(edged_dataOut_nonunif_mapc2p, zmax)

    # Convert from cell center to edges
    diffs  = dataOut_unif_mapc2p[0:-1] + np.diff(dataOut_unif_mapc2p)/2
    edged_dataOut_unif_mapc2p = np.insert(diffs, 0, zmin)
    edged_dataOut_unif_mapc2p = np.append(edged_dataOut_unif_mapc2p, zmax)

    unif_Jgeo_filename = str(dataDir+unifFile+'-jacobgeo.gkyl')
    pgData_unif_Jgeo = pg.GData(unif_Jgeo_filename)
    pgInterp_unif_Jgeo = pg.GInterpModal(pgData_unif_Jgeo, polyOrder, 'ms')
    x_unif_Jgeo, dataOut_unif_Jgeo = pgInterp_unif_Jgeo.interpolate()
    dataOut_unif_Jgeo = np.squeeze(dataOut_unif_Jgeo)

    nonunif_Jgeo_filename = str(dataDir+nonunifFile+'-jacobgeo.gkyl')
    pgData_nonunif_Jgeo = pg.GData(nonunif_Jgeo_filename)
    pgInterp_nonunif_Jgeo = pg.GInterpModal(pgData_nonunif_Jgeo, polyOrder, 'ms')
    x_nonunif_Jgeo, dataOut_nonunif_Jgeo = pgInterp_nonunif_Jgeo.interpolate()
    dataOut_nonunif_Jgeo = np.squeeze(dataOut_nonunif_Jgeo)

    unif_distf_shape = dataOut_unif.shape
    nonunif_distf_shape = dataOut_nonunif.shape

    tile_unif_Jgeo = np.ones((unif_distf_shape[0], unif_distf_shape[1], unif_distf_shape[2]))
    for i in range(unif_distf_shape[0]):
        tile_unif_Jgeo[i,:,:] *= dataOut_unif_Jgeo[i]

    tile_nonunif_Jgeo = np.ones((nonunif_distf_shape[0], nonunif_distf_shape[1], nonunif_distf_shape[2]))
    for i in range(nonunif_distf_shape[0]):
        tile_nonunif_Jgeo[i,:,:] *= dataOut_nonunif_Jgeo[i]

    dataOut_unif = np.trapz(np.abs(dataOut_unif / tile_unif_Jgeo), axis=2)
    dataOut_nonunif = np.trapz(np.abs(dataOut_nonunif / tile_nonunif_Jgeo), axis=2)

    # Interpolate the non-uniform data onto a uniform grid
    dataOut_unif_shape = dataOut_unif.shape
    dataOut_nonunif_interp = np.zeros((dataOut_unif_shape[0], dataOut_unif_shape[1]))
    for i in range(unif_distf_shape[1]):
        dataOut_nonunif_interp[:,i] = np.interp(dataOut_unif_mapc2p, dataOut_nonunif_mapc2p, dataOut_nonunif[:,i])
    data_difference = np.abs((dataOut_unif - dataOut_nonunif_interp))
  #   # Need velocity space grids or to convert to edges of z for plotting

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10,12))

    norm = LogNorm(vmin = 1e-10, vmax = 1000)  # Create a LogNorm instance

    pcolormesh1 = ax1.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[1]/vti, dataOut_unif.T, cmap='inferno', norm=norm)
    fig.colorbar(pcolormesh1, ax=ax1)  # Add a colorbar to the plot
    # Label the axes
    ax1.set_ylabel('vpar / vti')
    ax1.set_xlabel('Uniform grid Z, cylindrical coodinate (m)')
    ax1.set_title('Frame '+str(frameNum))

    # pcolormesh2 = ax2.pcolormesh(edged_dataOut_nonunif_mapc2p, x_nonunif[1]/vti, dataOut_nonunif.T, cmap='inferno', norm=norm)
    pcolormesh2 = ax2.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[1]/vti, dataOut_nonunif_interp.T, cmap='inferno', norm=norm)
    fig.colorbar(pcolormesh2, ax=ax2)  # Add a colorbar to the plot
    # Label the axes
    ax2.set_ylabel('vpar / vti')
    ax2.set_xlabel('Nonuniform grid Z, cylindrical coodinate (m)') 
    pcolormesh3 = ax3.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[1] / vti, data_difference.T, cmap='inferno', norm=norm)
    fig.colorbar(pcolormesh3, ax=ax3)  # Add a colorbar to the plot
  # Label the axes
    ax3.set_ylabel('vpar / vti')
    ax3.set_xlabel('Absolute difference between distribution functions Z, cylindrical coordinate (m)')
    
    figName = 'distf_vpar_'+str(frameNum)
    if save_figure_as_file:
      fig.savefig(outDir+figName+figureFileFormat)
      plt.close()
      print('Saved figure as '+figName+figureFileFormat)  
    else:
      plt.show()  

  if plot_distmu:
    # f_unif, f_map, f_nonunif, f_nonunif_map = load_mapped_data('-ion-')
    dataName = '-ion_'
    densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
    pgData_unif = pg.GData(densityFileName_unif)
    pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'ms')
    x_unif, dataOut_unif = pgInterp_unif.interpolate()
    dataOut_unif = np.squeeze(dataOut_unif)

    densityFileName_nonunif = str(dataDir+nonunifFile + str(dataName) + str(frameNum) + '.gkyl')
    pgData_nonunif = pg.GData(densityFileName_nonunif)
    pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'ms')
    x_nonunif, dataOut_nonunif = pgInterp_nonunif.interpolate()
    dataOut_nonunif = np.squeeze(dataOut_nonunif)

    nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
    pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
    pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
    x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pgInterp_nonunif_mapc2p.interpolate(2)
    dataOut_nonunif_mapc2p = np.squeeze(dataOut_nonunif_mapc2p)

    unif_mapc2p_filename = str(dataDir+unifFile+'-mapc2p.gkyl')
    pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
    pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
    x_unif_mapc2p, dataOut_unif_mapc2p = pgInterp_unif_mapc2p.interpolate(2)
    dataOut_unif_mapc2p = np.squeeze(dataOut_unif_mapc2p)

    # Convert from cell center to edges
    zmin = x_unif_mapc2p[0][0]
    zmax = x_unif_mapc2p[0][-1]
    diffs  = dataOut_nonunif_mapc2p[0:-1] + np.diff(dataOut_nonunif_mapc2p)/2
    edged_dataOut_nonunif_mapc2p = np.insert(diffs, 0, zmin)
    edged_dataOut_nonunif_mapc2p = np.append(edged_dataOut_nonunif_mapc2p, zmax)

    # Convert from cell center to edges
    diffs  = dataOut_unif_mapc2p[0:-1] + np.diff(dataOut_unif_mapc2p)/2
    edged_dataOut_unif_mapc2p = np.insert(diffs, 0, zmin)
    edged_dataOut_unif_mapc2p = np.append(edged_dataOut_unif_mapc2p, zmax)

    unif_Jgeo_filename = str(dataDir+unifFile+'-jacobgeo.gkyl')
    pgData_unif_Jgeo = pg.GData(unif_Jgeo_filename)
    pgInterp_unif_Jgeo = pg.GInterpModal(pgData_unif_Jgeo, polyOrder, 'ms')
    x_unif_Jgeo, dataOut_unif_Jgeo = pgInterp_unif_Jgeo.interpolate()
    dataOut_unif_Jgeo = np.squeeze(dataOut_unif_Jgeo)

    nonunif_Jgeo_filename = str(dataDir+nonunifFile+'-jacobgeo.gkyl')
    pgData_nonunif_Jgeo = pg.GData(nonunif_Jgeo_filename)
    pgInterp_nonunif_Jgeo = pg.GInterpModal(pgData_nonunif_Jgeo, polyOrder, 'ms')
    x_nonunif_Jgeo, dataOut_nonunif_Jgeo = pgInterp_nonunif_Jgeo.interpolate()
    dataOut_nonunif_Jgeo = np.squeeze(dataOut_nonunif_Jgeo)

    unif_distf_shape = dataOut_unif.shape
    nonunif_distf_shape = dataOut_nonunif.shape

    tile_unif_Jgeo = np.ones((unif_distf_shape[0], unif_distf_shape[1], unif_distf_shape[2]))
    for i in range(unif_distf_shape[0]):
        tile_unif_Jgeo[i,:,:] *= dataOut_unif_Jgeo[i]

    tile_nonunif_Jgeo = np.ones((nonunif_distf_shape[0], nonunif_distf_shape[1], nonunif_distf_shape[2]))
    for i in range(nonunif_distf_shape[0]):
        tile_nonunif_Jgeo[i,:,:] *= dataOut_nonunif_Jgeo[i]

    dataOut_unif = np.trapz(np.abs(dataOut_unif / tile_unif_Jgeo), axis=1)
    dataOut_nonunif = np.trapz(np.abs(dataOut_nonunif / tile_nonunif_Jgeo), axis=1)

    # Interpolate the non-uniform data onto a uniform grid
    dataOut_unif_shape = dataOut_unif.shape
    dataOut_nonunif_interp = np.zeros((dataOut_unif_shape[0], dataOut_unif_shape[1]))
    for i in range(dataOut_unif_shape[1]):
        dataOut_nonunif_interp[:,i] = np.interp(dataOut_unif_mapc2p, dataOut_nonunif_mapc2p, dataOut_nonunif[:,i])
    data_difference = np.abs((dataOut_unif - dataOut_nonunif_interp))

    # Need velocity space grids or to convert to edges of z for plotting
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10,12))

    norm = LogNorm(vmin = 1e-4, vmax = 1000)  # Create a LogNorm instance

    pcolormesh1 = ax1.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[2]/mui0, dataOut_unif.T, cmap='inferno', norm=norm)
    fig.colorbar(pcolormesh1, ax=ax1)  # Add a colorbar to the plot
    # Label the axes
    ax1.set_ylabel('mu / mui0')
    ax1.set_xlabel('Uniform grid Z, cylindrical coodinate (m)')
    ax1.set_title('Frame '+str(frameNum))

    pcolormesh2 = ax2.pcolormesh(edged_dataOut_nonunif_mapc2p, x_nonunif[2]/mui0, dataOut_nonunif.T, cmap='inferno', norm=norm)
    fig.colorbar(pcolormesh2, ax=ax2)  # Add a colorbar to the plot
    # Label the axes
    ax2.set_ylabel('mu / mui0')
    ax2.set_xlabel('Nonuniform grid Z, cylindrical coodinate (m)')

    pcolormesh3 = ax3.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[2] / vti, data_difference.T, cmap='inferno', norm=norm)
    fig.colorbar(pcolormesh3, ax=ax3)  # Add a colorbar to the plot
  # Label the axes
    ax3.set_ylabel('vpar / vti')
    ax3.set_xlabel('Absolute difference between distribution functions Z, cylindrical coordinate (m)')
    
    figName = 'distf_mu_'+str(frameNum)
    if save_figure_as_file:
      fig.savefig(outDir+figName+figureFileFormat)
      plt.close()
      print('Saved figure as '+figName+figureFileFormat)  
    else:
      plt.show() 
      
  def plot_distf_at_z_eq(z0_coordinate):
    dataName = '-ion_'
    densityFileName_unif = str(dataDir+unifFile + str(dataName) + str(frameNum) + '.gkyl')
    pgData_unif = pg.GData(densityFileName_unif)
    pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'ms')
    pgInterp_unif.interpolate(overwrite=True)
    x_unif, dataOut_unif = pg.data.select(pgData_unif, z0 = z_xi(z0_coordinate))
    dataOut_unif = np.squeeze(dataOut_unif)

    densityFileName_nonunif = str(dataDir+nonunifFile + str(dataName) + str(frameNum) + '.gkyl')
    pgData_nonunif = pg.GData(densityFileName_nonunif)
    pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'ms')
    pgInterp_nonunif.interpolate(overwrite=True)
    x_nonunif, dataOut_nonunif = pg.data.select(pgData_nonunif, z0 = (z0_coordinate))
    dataOut_nonunif = np.squeeze(dataOut_nonunif)

    nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
    pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
    pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
    pgInterp_nonunif_mapc2p.interpolate(overwrite=True)
    x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pg.data.select(pgData_nonunif_mapc2p, z0 = (z0_coordinate))
    dataOut_nonunif_mapc2p = np.squeeze(dataOut_nonunif_mapc2p)

    unif_mapc2p_filename = str(dataDir+unifFile+'-mapc2p.gkyl')
    pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
    pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
    pgInterp_unif_mapc2p.interpolate(overwrite=True)
    x_unif_mapc2p, dataOut_unif_mapc2p = pg.data.select(pgData_unif_mapc2p, z0 = z_xi(z0_coordinate))
    dataOut_unif_mapc2p = np.squeeze(dataOut_unif_mapc2p)

    unif_Jgeo_filename = str(dataDir+unifFile+'-jacobgeo.gkyl')
    pgData_unif_Jgeo = pg.GData(unif_Jgeo_filename)
    pgInterp_unif_Jgeo = pg.GInterpModal(pgData_unif_Jgeo, polyOrder, 'ms')
    pgInterp_unif_Jgeo.interpolate(overwrite=True)
    x_unif_Jgeo, dataOut_unif_Jgeo = pg.data.select(pgData_unif_Jgeo, z0 = z_xi(z0_coordinate))
    dataOut_unif_Jgeo = np.squeeze(dataOut_unif_Jgeo)

    nonunif_Jgeo_filename = str(dataDir+nonunifFile+'-jacobgeo.gkyl')
    pgData_nonunif_Jgeo = pg.GData(nonunif_Jgeo_filename)
    pgInterp_nonunif_Jgeo = pg.GInterpModal(pgData_nonunif_Jgeo, polyOrder, 'ms')
    pgInterp_nonunif_Jgeo.interpolate(overwrite=True)
    x_nonunif_Jgeo, dataOut_nonunif_Jgeo = pg.data.select(pgData_nonunif_Jgeo, z0 = (z0_coordinate))
    dataOut_nonunif_Jgeo = np.squeeze(dataOut_nonunif_Jgeo)

    dataOut_unif = np.abs(dataOut_unif / dataOut_unif_Jgeo)
    dataOut_nonunif = np.abs(dataOut_nonunif / dataOut_nonunif_Jgeo)
    dataOut_diff = np.abs((dataOut_unif - dataOut_nonunif))

    # Need velocity space grids or to convert to edges of z for plotting
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18,6))

    norm = LogNorm(vmin = 1e-6, vmax = np.max(dataOut_unif))  # Create a LogNorm instance

    pcolormesh1 = ax1.pcolormesh(x_unif[1]/vti, x_unif[2]/mui0, dataOut_unif.T, cmap='inferno', norm=norm)
    fig.colorbar(pcolormesh1, ax=ax1)  # Add a colorbar to the plot
    # Label the axes
    ax1.set_ylabel('mu / mui0')
    ax1.set_xlabel('vpar / vti')
    ax1.set_title('Frame '+str(frameNum))

    pcolormesh2 = ax2.pcolormesh(x_nonunif[1]/vti, x_nonunif[2]/mui0, dataOut_nonunif.T, cmap='inferno', norm=norm)
    fig.colorbar(pcolormesh2, ax=ax2)  # Add a colorbar to the plot
    # Label the axes
    ax2.set_ylabel('mu / mui0')
    ax2.set_xlabel('vpar / vti')
    
    pcolormesh3 = ax3.pcolormesh(x_unif[1]/vti, x_unif[2] / vti, dataOut_diff.T, cmap='inferno', norm=norm)
    fig.colorbar(pcolormesh3, ax=ax3)  # Add a colorbar to the plot
    # Label the axes
    ax3.set_ylabel('vpar / vti')
    ax3.set_xlabel('Absolute difference between distribution functions Z, cylindrical coordinate (m)')

    ax2.title.set_text('z = '+str(z_xi(z0_coordinate))+' m')

    locName = np.round(z0_coordinate*100)
    figName = 'distf_z'+str(locName)+'_'+str(frameNum)
    if save_figure_as_file:
      fig.savefig(outDir+figName+figureFileFormat)
      plt.close()
      print('Saved figure as '+figName+figureFileFormat)  
    else:
      plt.show() 

    # On axis difference
    reduce_set_unif = dataOut_unif[:,0]
    reduce_set_nonunif = dataOut_nonunif[:,0]
    x_axis = x_unif[1]/vti

    # I need to interpolate for x 
    plt.plot(x_axis[1:], reduce_set_unif, label='Uniform grid')
    plt.plot(x_axis[1:], reduce_set_nonunif, label='Nonuniform grid')
    plt.plot(x_axis[1:], np.abs(reduce_set_unif - reduce_set_nonunif), 'k--', label='Absolute difference')
    plt.xlabel('vpar / vti')
    plt.ylim(np.max(reduce_set_nonunif)* 1e-6, np.max(reduce_set_nonunif)*1.5)
    plt.xlim(-2, 2)

    plt.ylabel('f')
    # log scale y axis
    plt.yscale('log')
    plt.legend()
    plt.title('On axis difference, mu=0, frame '+str(frameNum)+ ', z = '+str(z_xi(z0_coordinate))+' m')

    figName = 'distf_z'+str(locName)+'mu0'+str(frameNum)
    if save_figure_as_file:
      plt.savefig(outDir+figName+figureFileFormat)
      plt.close()
      print('Saved figure as '+figName+figureFileFormat)
    else:
      plt.show()

  if plot_distf_at_z:
    for z0 in z_loctions:
      plot_distf_at_z_eq(z0)


process_frame(0)
# # Number of processes to run in parallel
# num_processes = multiprocessing.cpu_count()
# print('Number of processes: ', num_processes)

# # Create a pool of processes
# pool = multiprocessing.Pool(processes=num_processes)

# # Map the frame_arr to the pool of processes
# pool.map(process_frame, frame_arr)

# # Close the pool to prevent any more tasks from being submitted
# pool.close()

# # Wait for all processes to finish
# pool.join()
  