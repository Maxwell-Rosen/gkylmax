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

dataDir = '/scratch/gpfs/mr1884/scratch/gkylmax/mirror1x_compare_unif_vs_nonunif_grids/outputs/'
unifFile = 'gk_mirror_adiabatic_elc_1x2v_p1_nosource_uniform'
nonunifFile = 'gk_mirror_adiabatic_elc_1x2v_p1_nosource_nonuniform'
frameNum = 0                 #.Frame number to process.

plot_density       = True  #[ Plot density.
plot_phiAdiabatic  = False  #[ Plot e*phi(z,t)/Te0 and e*phi(z=0,t)/Te0 for the adiabatic elc sim.

outDir = './python_plots'

outFigureFile    = False     #[ Output a figure file?.
figureFileFormat = '.png'    #[ Can be .png, .pdf, .ps, .eps, .svg.
saveData         = False    #.Indicate whether to save data in plot to HDF5 file.



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
c_s = np.sqrt(Te0/mi)
print("vti0 = ",vti)
print("mui0 = ",0.5*mi*(vti**2)/B_p)

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

#................................................................................#

if plot_density:
  densityFileName_unif = str(dataDir+unifFile+'-ion_M0_' + str(frameNum) + '.gkyl')
  pgData_unif = pg.GData(densityFileName_unif)
  pgInterp_unif = pg.GInterpModal(pgData_unif, polyOrder, 'ms')
  x_unif, dataOut_unif = pgInterp_unif.interpolate()

  densityFileName_nonunif = str(dataDir+nonunifFile+'-ion_M0_' + str(frameNum) + '.gkyl')
  pgData_nonunif = pg.GData(densityFileName_nonunif)
  pgInterp_nonunif = pg.GInterpModal(pgData_nonunif, polyOrder, 'ms')
  x_nonunif, dataOut_nonunif = pgInterp_nonunif.interpolate()

  nonunif_mapc2p_filename = str(dataDir+'mapc2p_nonunif.gkyl')
  pgData_nonunif_mapc2p = pg.GData(nonunif_mapc2p_filename)
  pgInterp_nonunif_mapc2p = pg.GInterpModal(pgData_nonunif_mapc2p, polyOrder, 'ms')
  x_nonunif_mapc2p, dataOut_nonunif_mapc2p = pgInterp_nonunif_mapc2p.interpolate(2)
  avg_x = x_nonunif_mapc2p[0][:-1] + 0.5*np.diff(x_nonunif_mapc2p[0])
  plt.plot(avg_x, dataOut_nonunif_mapc2p)
  plt.show()
        #.Complete file name.
  # #[ Plot phi(z) and e*Phi/Te0 for a single frame, or phi(z,t) and phi(z=0,t).

  # phiFile = dataDir+fileName+'-field_'    #.Complete file name.
  # nFrames = 1+pgu.findLastFrame(phiFile, '.gkyl')
  # times = 1e-6*np.arange(0, nFrames, 1)  #1.e6*pgu.getTimeStamps(phiFile,0,nFrames-1, '.gkyl')

  # denFile = dataDir+fileName+'_ion_gridDiagnostics_%d.gkyl'    #.Complete file name.
  # #[ Load the grid.

  # # xInt, _, nxInt, lxInt, dxInt, _ = pgu.getGrid(denFile % 0,polyOrder,basisType,varName='M0')
  # # xIntC, _, nxIntC, lxIntC, dxIntC, _ = pgu.getGrid(denFile % 0,polyOrder,basisType,varName='M0',location='center')
  # firstFile = str(fileName + '-field_0.gkyl')
  # print(firstFile)
  # pgData = pg.GData(firstFile)
  # pgInterp = pg.GInterpNodal(pgData, polyOrder, 'ns')
  # print(pgInterp)
  # xInt, dataOut = pgInterp.interpolate()

  # den_i = np.squeeze(pgu.getInterpData(denFile % (nFrames-1), polyOrder, basisType, varName='M0'))

  # #[ Prepare figure.
  # figProp1 = (6.4,3.5)
  # ax1Pos   = [[0.16, 0.18, 0.82, 0.80],]
  # fig1     = plt.figure(figsize=figProp1)
  # ax1      = [fig1.add_axes(pos) for pos in ax1Pos]

  # ax1[0].semilogy(xInt[0], den_i)

  # #[ Plot central value over time:
  # ax1[0].set_xlim( np.amin(xInt[0]),np.amax(xInt[0]) )
  # ax1[0].set_xlabel('Length along field line, $z$ (m)', fontsize=xyLabelFontSize)
  # ax1[0].set_ylabel('$n_i$ (m$^{-3}$)', fontsize=xyLabelFontSize)
  # for i in range(len(ax1)):
  #   setTickFontSize(ax1[i],tickFontSize)

  # figName = fileName+'_ion_M0_'+str(nFrames-1)
  # if outFigureFile:
  #   fig1.savefig(outDir+figName+figureFileFormat)
  #   plt.close()
  # else:
  #   plt.show()

#................................................................................#

if plot_phiAdiabatic:
  #[ Plot phi(z) and e*Phi/Te0 for a single frame, or phi(z,t) and phi(z=0,t).

  # dataDir = '/scratch/gpfs/manaurer/gkeyll/mirror/gk57/'
  # fileName = 'gk57-wham1x2v'    #.Root name of files to process.

  #[ Prepare figure.
  figProp7 = (6.4,4.5)
  ax7Pos   = [[0.12, 0.57, 0.725, 0.42],
              [0.12, 0.13, 0.725, 0.42]]
  cax7Pos  = [[0.855, 0.57, 0.02, 0.42]]
  fig7     = plt.figure(figsize=figProp7)
  ax7      = [fig7.add_axes(pos) for pos in ax7Pos]
  cbar_ax7 = [fig7.add_axes(pos) for pos in cax7Pos]

  phiFile = dataDir+fileName+'-field_'    #.Complete file name.
  nFrames = 1+pgu.findLastFrame(phiFile, '.gkyl')
  times = 1.e6*pgu.getTimeStamps(phiFile,0,nFrames-1, '.gkyl')

  phiFile = dataDir+fileName+'-field_%d.gkyl'    #.Complete file name.

  #[ Load the grid.
  xInt, _, nxInt, lxInt, dxInt, _ = pgu.getGrid(phiFile % 0,polyOrder,basisType)
  xIntC, _, nxIntC, lxIntC, dxIntC, _ = pgu.getGrid(phiFile % 0,polyOrder,basisType,location='center')

  phi = np.zeros((nFrames,nxIntC[0]))

  for fI in range(nFrames):
    phi[fI,:] = np.squeeze(pgu.getInterpData(phiFile % fI, polyOrder, basisType))

  #[ Create colorplot grid (coordinates have to be nodal).
  timesC  = 0.5*(times[1:]+times[0:-1])
  Xnodal = [np.outer(np.concatenate([[timesC[0]-(timesC[1]-timesC[0])],timesC,[timesC[-1]+(timesC[-1]-timesC[-2])]]), \
                     np.ones(np.shape(xInt[0]))), \
            np.outer(np.ones(np.size(timesC)+2),xInt[0])]

  hpl7a = ax7[0].pcolormesh(Xnodal[0], Xnodal[1], eV*phi/Te0, cmap='inferno')
  hcb7a = plt.colorbar(hpl7a, ax=ax7[0], cax=cbar_ax7[0], ticks=[0., 2.5, 5., 7.5, 10.])
  #[ Plots lines at mirror throat:
  ax7[0].plot([Xnodal[0][0][0], Xnodal[0][-1][0]],   [z_m, z_m], color='white', linestyle='--', alpha=0.5)
  ax7[0].plot([Xnodal[0][0][0], Xnodal[0][-1][0]], [-z_m, -z_m], color='white', linestyle='--', alpha=0.5)

  #[ Plot central value over time:
  hpl7b = ax7[1].plot(times, (eV/Te0)*0.5*(phi[:,nxIntC[0]//2-1]+phi[:,nxIntC[0]//2]), color=defaultColors[0])
  ax7[1].set_xlim( np.amin(times),np.amax(times) )
  ax7[1].set_ylim( 0., 12. ) #np.amin(y),np.amax(y) )
  hcb7a.set_label('$e\phi/T_{e0}$', rotation=270, labelpad=18, fontsize=colorBarLabelFontSize)
  hcb7a.ax.tick_params(labelsize=tickFontSize)
  plt.setp( ax7[0].get_xticklabels(), visible=False)
  ax7[0].set_ylabel(r'$z$ (m)', fontsize=xyLabelFontSize)
  ax7[1].set_xlabel(r'Time ($\mu$s)', fontsize=xyLabelFontSize)
  ax7[1].set_ylabel(r'$e\phi(z=0)/T_{e0}$', fontsize=xyLabelFontSize)
  plt.text(0.025, 0.85, r'(a)', fontsize=textFontSize, color='white', transform=ax7[0].transAxes)
  plt.text(0.025, 0.85, r'(b)', fontsize=textFontSize, color='black', transform=ax7[1].transAxes)
  for i in range(len(ax7)):
    setTickFontSize(ax7[i],tickFontSize)

  figName = fileName+'_phiVtime_0-'+str(nFrames-1)
  if outFigureFile:
    fig7.savefig(outDir+figName+figureFileFormat)
    plt.close()
  else:
    plt.show()

#................................................................................#
