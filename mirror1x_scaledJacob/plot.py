#[ ........................................................... ]#
#[
#[ Plot ion density and potential from a Boltzmann electron
#[ 1x2v mirror simulation.
#[
#[
#[ Manaure Francisquez.
#[ January 2024
#[
#[ ........................................................... ]#

import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from matplotlib.colors import LogNorm

dataDir = '/scratch/gpfs/mr1884/scratch/gkylmax/'
unifFile = 'mirror1x_compare_unif_vs_nonunif_grids/outputs/gk_mirror_adiabatic_elc_1x2v_p1_nosource_uniform'
nonunifFile = 'mirror1x_scaledJacob/outputs/gk_mirror_adiabatic_elc_1x2v_p1_nosource_nonuniform'
frameNum = 14

plot_moments       = 1  # Plot density, potential, upar, tperp, tpar.
plot_distvpar      = 1  # plot distribution function in vpar.
plot_distmu        = 1  # plot distribution function in mu.

outDir = './python_plots/'

outFigureFile    = 1     #[ If True, save figure to file. If False, display figure on screen.
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

  nonunif_mapc2p_filename = str(dataDir+nonunifFile+'-mapc2p.gkyl')
  pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
  pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
  x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pgInterp_nonunif_mapc2p.interpolate(2)

  unif_mapc2p_filename = str(dataDir+unifFile+'-mapc2p.gkyl')
  pgData_unif_mapc2p = pg.GData(unif_mapc2p_filename)
  pgInterp_unif_mapc2p = pg.GInterpModal(pgData_unif_mapc2p, polyOrder, 'ms')
  x_unif_mapc2p, dataOut_unif_mapc2p = pgInterp_unif_mapc2p.interpolate(2)

  return dataOut_unif, dataOut_unif_mapc2p, dataOut_nonunif, dataOut_nonunif_mapc2p

#................................................................................#

if plot_moments:
  dataOut_unif, dataOut_unif_mapc2p, dataOut_nonunif, dataOut_nonunif_mapc2p = load_mapped_data('-ion_M0_')
  
  # Create a subfigure that is 2 by 3
  fig, ax = plt.subplots(2, 3, figsize=(20,10))
  # Plot the density
  ax[0,0].plot(dataOut_unif_mapc2p[:,0], dataOut_unif[:,0],'r', label='Uniform grid')
  ax[0,0].plot(dataOut_nonunif_mapc2p[:,0], dataOut_nonunif[:,0],'b--', label='Nonuniform grid')
  ax[0,0].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
  ax[0,0].set_ylabel('$n_i$ (m$^{-3}$)', fontsize=xyLabelFontSize)
  ax[0,0].legend(loc='upper left', fontsize=legendFontSize)
  setTickFontSize(ax[0,0],tickFontSize)

  # Plot phi
  dataOut_unif, dataOut_unif_mapc2p, dataOut_nonunif, dataOut_nonunif_mapc2p = load_mapped_data('-field_')

  dataOut_unif *= eV/Te0
  dataOut_nonunif *= eV/Te0

  ax[0,1].plot(dataOut_unif_mapc2p[:,0], dataOut_unif[:,0],'r', label='Uniform grid')
  ax[0,1].plot(dataOut_nonunif_mapc2p[:,0], dataOut_nonunif[:,0],'b--', label='Nonuniform grid')
  ax[0,1].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
  ax[0,1].set_ylabel('$\phi$ (m$^{-3}$)', fontsize=xyLabelFontSize)
  ax[0,1].legend(loc='upper left', fontsize=legendFontSize)
  setTickFontSize(ax[0,1],tickFontSize)

  # Plot uPar
  M0_unif, M0_map, M0_nonunif, M0_nonunif_map = load_mapped_data('-ion_M0_')
  M1_unif, M1_map, M1_nonunif, M1_nonunif_map = load_mapped_data('-ion_M1_')

  upar_unif = M1_unif[:,0]/M0_unif[:,0]
  upar_nonunif = M1_nonunif[:,0]/M0_nonunif[:,0]

  ax[1,0].plot(M0_map[:,0], upar_unif / c_s,'r', label='Uniform grid')
  ax[1,0].plot(M0_nonunif_map[:,0], upar_nonunif / c_s,'b--', label='Nonuniform grid')
  ax[1,0].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
  ax[1,0].set_ylabel('$u_{\parallel} / c_s$ (m/s)', fontsize=xyLabelFontSize)
  ax[1,0].legend(loc='upper left', fontsize=legendFontSize)
  setTickFontSize(ax[1,0],tickFontSize)
  
  # Plot tPerp
  M0_unif, M0_map, M0_nonunif, M0_nonunif_map = load_mapped_data('-ion_M0_')
  M2perp_unif, M2perp_map, M2perp_nonunif, M2perp_nonunif_map = load_mapped_data('-ion_M2perp_')

  tPerp_unif = M2perp_unif[:,0]/M0_unif[:,0] * mi / eV
  tPerp_nonunif = M2perp_nonunif[:,0]/M0_nonunif[:,0] * mi / eV

  ax[1,1].plot(M0_map[:,0], tPerp_unif,'r', label='Uniform grid')
  ax[1,1].plot(M0_nonunif_map[:,0], tPerp_nonunif,'b--', label='Nonuniform grid')
  ax[1,1].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
  ax[1,1].set_ylabel('$T_{\perp}$ (eV)', fontsize=xyLabelFontSize)
  ax[1,1].legend(loc='upper left', fontsize=legendFontSize)
  setTickFontSize(ax[1,1],tickFontSize)

  # Plot tPar
  M0_unif, M0_map, M0_nonunif, M0_nonunif_map = load_mapped_data('-ion_M0_')
  M1_unif, M1_map, M1_nonunif, M1_nonunif_map = load_mapped_data('-ion_M1_')
  M2par_unif, M2par_map, M2par_nonunif, M2par_nonunif_map = load_mapped_data('-ion_M2par_')

  tPar_unif = (M2par_unif[:,0] - M1_unif[:,0]**2/M0_unif[:,0]) * mi / eV / M0_unif[:,0]
  tPar_nonunif = (M2par_nonunif[:,0] - M1_nonunif[:,0]**2/M0_nonunif[:,0]) * mi / eV / M0_nonunif[:,0]

  ax[1,2].plot(M0_map[:,0], tPar_unif,'r', label='Uniform grid')
  ax[1,2].plot(M0_nonunif_map[:,0], tPar_nonunif,'b--', label='Nonuniform grid')
  ax[1,2].set_xlabel('Cylindrical length coordinate, $Z$ (m)', fontsize=xyLabelFontSize)
  ax[1,2].set_ylabel('$T_{\parallel}$ (eV)', fontsize=xyLabelFontSize)
  ax[1,2].legend(loc='upper left', fontsize=legendFontSize)
  setTickFontSize(ax[1,2],tickFontSize)
  
  figName = 'moments_'+str(frameNum)
  if outFigureFile:
    fig.savefig(outDir+figName+figureFileFormat)
    plt.close()
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

  # Need velocity space grids or to convert to edges of z for plotting

  fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,8))

  norm = LogNorm(vmin = 1e-10, vmax = 1000)  # Create a LogNorm instance

  pcolormesh1 = ax1.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[1]/vti, dataOut_unif.T, cmap='inferno', norm=norm)
  fig.colorbar(pcolormesh1, ax=ax1)  # Add a colorbar to the plot
  # Label the axes
  ax1.set_ylabel('vpar / vti')
  ax1.set_xlabel('Uniform grid Z, cylindrical coodinate (m)')

  pcolormesh2 = ax2.pcolormesh(edged_dataOut_nonunif_mapc2p, x_nonunif[1]/vti, dataOut_nonunif.T, cmap='inferno', norm=norm)
  fig.colorbar(pcolormesh2, ax=ax2)  # Add a colorbar to the plot
  # Label the axes
  ax2.set_ylabel('vpar / vti')
  ax2.set_xlabel('Nonuniform grid Z, cylindrical coodinate (m)')
  
  figName = 'distf_vpar_'+str(frameNum)
  if outFigureFile:
    fig.savefig(outDir+figName+figureFileFormat)
    plt.close()
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

  # Need velocity space grids or to convert to edges of z for plotting

  fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,8))

  norm = LogNorm(vmin = 1e-4, vmax = 1000)  # Create a LogNorm instance

  pcolormesh1 = ax1.pcolormesh(edged_dataOut_unif_mapc2p, x_unif[2]/mui0, dataOut_unif.T, cmap='inferno', norm=norm)
  fig.colorbar(pcolormesh1, ax=ax1)  # Add a colorbar to the plot
  # Label the axes
  ax1.set_ylabel('mu / mui0')
  ax1.set_xlabel('Uniform grid Z, cylindrical coodinate (m)')

  pcolormesh2 = ax2.pcolormesh(edged_dataOut_nonunif_mapc2p, x_nonunif[2]/mui0, dataOut_nonunif.T, cmap='inferno', norm=norm)
  fig.colorbar(pcolormesh2, ax=ax2)  # Add a colorbar to the plot
  # Label the axes
  ax2.set_ylabel('mu / mui0')
  ax2.set_xlabel('Nonuniform grid Z, cylindrical coodinate (m)')
  
  figName = 'distf_mu_'+str(frameNum)
  if outFigureFile:
    fig.savefig(outDir+figName+figureFileFormat)
    plt.close()
  else:
    plt.show() 