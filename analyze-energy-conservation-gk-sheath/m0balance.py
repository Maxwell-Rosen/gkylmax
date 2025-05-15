#[ ........................................................... ]#
#[
#[ Check particle balance
#[
#[ Manaure Francisquez.
#[ Dec 2024
#[
#[ ........................................................... ]#


import numpy as np
import postgkyl as pg
import matplotlib.pyplot as plt
import sys
# sys.path.insert(0, '/Users/mfrancis/Documents/codebits/pets/gkeyll/postProcessingScripts/')
import pgkylUtil as pgu

dataDir = './'
outDir = './conservation_plots/'
outFileRoot = outDir + 'particle_balance_'

simName = 'gk_sheath_1x2v_p1'

cfl_dirs = ['on_main/', 'thresh-1e-10/', 'thresh-1e-8/', 'thresh-1e-5/', 'thresh-1e-3/']
legendStrings_error = [r'Main',r'Threshold = Maximum*1e-10',r'Threshold = Maximum*1e-8',r'Threshold = Maximum*1e-5',r'Threshold = Maximum*1e-3']
legendStrings_diff = [r'Threshold = Maximum*1e-10',r'Threshold = Maximum*1e-8',r'Threshold = Maximum*1e-5',r'Threshold = Maximum*1e-3']

plot_balance = True  #[ Balance of various terms.
plot_error   = True  #[ Conservation error.
plot_m0_wrt_main = False  #[ M0 moment of the distribution function with respect to the main one.
plot_m2_wrt_main = False  #[ M2 moment of the distribution function with respect to the main one.  

outFigureFile    = True   #[ Output a figure file?.
figureFileFormat = '.png'  #[ Can be .png, .pdf, .ps, .eps, .svg.


#[ ............... End of user inputs (MAYBE) ..................... ]#

#[ Some RGB colors. These are MATLAB-like.
defaultBlue    = [0, 0.4470, 0.7410]
defaultOrange  = [0.8500, 0.3250, 0.0980]
defaultGreen   = [0.4660, 0.6740, 0.1880]
defaultPurple  = [0.4940, 0.1840, 0.5560]
# defaultRed     = [0.6350, 0.0780, 0.1840]
defaultSkyBlue = [0.3010, 0.7450, 0.9330]
grey           = [0.5, 0.5, 0.5]
#[ Colors in a single array.
defaultColors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
#[ LineStyles in a single array.
lineStyles = ['-','--',':','-.','-', '--', ':', '-.']
markers    = ['None','None','None','None','None','None','None','None']

#[ Some fontsizes used in plots.
xyLabelFontSize       = 17
titleFontSize         = 17
colorBarLabelFontSize = 17
tickFontSize          = 14
legendFontSize        = 14

#.Set the font size of the ticks to a given size.
def setTickFontSize(axIn,fontSizeIn):
  axIn.tick_params(axis='both',labelsize=fontSizeIn)
  offset_txt = axIn.yaxis.get_offset_text() # Get the text object
  offset_txt.set_size(fontSizeIn) # # Set the size.
  offset_txt = axIn.xaxis.get_offset_text() # Get the text object
  offset_txt.set_size(fontSizeIn) # # Set the size.

#[ ............... End common utilities ..................... ]#

if plot_balance:
  #[ Read the M0 moment of df/dt, the source and the particle fluxes.
  
  for dI in range(len(cfl_dirs)):
    dataPath = dataDir + cfl_dirs[dI] + simName + '-'
    run_name = cfl_dirs[dI][:-1]
  
    species = 'elc'
    
    time_fdot, fdot = pgu.readDynVector(dataPath + species + '_fdot_integrated_moms.gkyl')
    time_bflux_lo, bflux_lo = pgu.readDynVector(dataPath + species + '_bflux_xlower_integrated_HamiltonianMoments.gkyl')
    time_bflux_up, bflux_up = pgu.readDynVector(dataPath + species + '_bflux_xupper_integrated_HamiltonianMoments.gkyl')
    time_src, src = pgu.readDynVector(dataPath + species + '_source_integrated_moms.gkyl')
    
    #[ Select the M0 moment.
    fdot = fdot[:,0]
    bflux_lo = bflux_lo[:,0]
    bflux_up = bflux_up[:,0]
    src = src[:,0]
    
    src[0] = 0.0 #[ Set source=0 at t=0 since we don't have fdot and bflux then.
    
    bflux_tot = bflux_lo+bflux_up #[ Total boundary flux loss.
    mom_err = src - bflux_tot - fdot #[ Error.
    
    #[ Plot each contribution.
    figProp1a = (10, 6)
    ax1aPos   = [0.09, 0.15, 0.87, 0.8]
    fig1a     = plt.figure(figsize=figProp1a)
    ax1a      = fig1a.add_axes(ax1aPos)
    
    hpl1a = list()
    hpl1a.append(ax1a.plot(time_fdot, fdot, color=defaultColors[0], linestyle=lineStyles[0], linewidth=2))
    hpl1a.append(ax1a.plot(time_bflux_lo, -(bflux_lo+bflux_up), color=defaultColors[1], linestyle=lineStyles[1], linewidth=2))
    hpl1a.append(ax1a.plot(time_src, src, color=defaultColors[2], linestyle=lineStyles[2], linewidth=2))
    hpl1a.append(ax1a.plot(time_fdot, mom_err, color=defaultColors[3], linestyle=lineStyles[3], linewidth=2))
    ax1a.set_xlabel(r'Time ($s$)',fontsize=xyLabelFontSize, labelpad=+4)
    ax1a.set_ylabel(r'0-th  moment',fontsize=xyLabelFontSize, labelpad=0)
    ax1a.set_xlim( time_fdot[0], time_fdot[-1] )
    ax1a.set_title(run_name, fontsize=titleFontSize)
    legendStrings = [r'$\dot{f}$',r'$\int_{\partial \Omega}\mathrm{d}\mathbf{S}\cdot\mathbf{\dot{R}}f$',r'$\mathcal{S}$',
    r'$E_{\dot{N}}=\mathcal{S}-\int_{\partial \Omega}\mathrm{d}\mathbf{S}\cdot\mathbf{\dot{R}}f-\dot{f}$',
    ]
    ax1a.legend([hpl1a[i][0] for i in range(len(hpl1a))], legendStrings, fontsize=legendFontSize, frameon=False)
    setTickFontSize(ax1a,tickFontSize)
    plt.tight_layout()
    
    if outFigureFile:
      plt.savefig(outFileRoot+run_name+'_particle'+figureFileFormat)
    else:
      plt.show()

#[ .......................................................... ]#

if plot_error:
  #[ Plot the error normalized for different time steps.
  eV = 1.602e-19
  me = 0.91e-30
  Lz = 4.0
  n0 = 7.0e18
  Te = 40.0*eV
  vte = np.sqrt(Te/me)
  tau = 0.5*Lz/vte #[ Transit time.
  
  figProp2a = (10, 6)
  ax2aPos   = [0.09, 0.15, 0.87, 0.8]
  fig2a     = plt.figure(figsize=figProp2a)
  ax2a      = fig2a.add_axes(ax2aPos)
  
  hpl2a = list()
  
  for dI in range(len(cfl_dirs)):
    dataPath = dataDir + cfl_dirs[dI] + simName + '-'
    
    time_fdot, fdot = pgu.readDynVector(dataPath + species + '_fdot_integrated_moms.gkyl')
    time_bflux_lo, bflux_lo = pgu.readDynVector(dataPath + species + '_bflux_xlower_integrated_HamiltonianMoments.gkyl')
    time_bflux_up, bflux_up = pgu.readDynVector(dataPath + species + '_bflux_xupper_integrated_HamiltonianMoments.gkyl')
    time_src, src = pgu.readDynVector(dataPath + species + '_source_integrated_moms.gkyl')
    time_distf, distf = pgu.readDynVector(dataPath + species + '_integrated_moms.gkyl')
    time_dt, dt = pgu.readDynVector(dataPath + 'dt.gkyl')
    
    #[ Select the M0 moment.
    fdot = fdot[:,0]
    bflux_lo = bflux_lo[:,0]
    bflux_up = bflux_up[:,0]
    src = src[:,0]
    distf = distf[:,0]
    
    bflux_tot = bflux_lo+bflux_up #[ Total boundary flux loss.
    mom_err = src - bflux_tot - fdot #[ Error.
  
    mom_err_norm = np.abs(mom_err[1:]*tau/distf[1:])
    
    hpl2a.append(ax2a.plot(time_dt, mom_err_norm, color=defaultColors[dI], linestyle=lineStyles[dI], linewidth=2))
  
  ax2a.set_xlabel(r'Time ($s$)',fontsize=xyLabelFontSize, labelpad=+4)
  ax2a.set_ylabel(r'$E_{\dot{N}}~\Delta t/N$',fontsize=xyLabelFontSize, labelpad=0)
  ax2a.set_xlim( time_fdot[0], time_fdot[-1] )
  ax2a.legend([hpl2a[i][0] for i in range(len(hpl2a))], legendStrings_error, fontsize=legendFontSize, frameon=False)
  ax2a.set_yscale('log')
  setTickFontSize(ax2a,tickFontSize)
  plt.tight_layout()
  
  if outFigureFile:
    plt.savefig(outFileRoot+'particle_change'+figureFileFormat)
  else:
    plt.show()


if plot_m0_wrt_main:
  #[ Plot the error normalized for different time steps.
  species_list = ['elc','ion']
  for sI in range(len(species_list)):
    species = species_list[sI]
    figProp2a = (12, 6)
    ax2aPos   = [0.09, 0.15, 0.87, 0.8]
    fig2a     = plt.figure(figsize=figProp2a)
    ax2a      = fig2a.add_axes(ax2aPos)
    
    hpl2a = list()
    
    for dI in range(len(cfl_dirs)):
      dataPath = dataDir + cfl_dirs[dI] + simName + '-'
      
      time_distf, distf = pgu.readDynVector(dataPath + species + '_integrated_moms.gkyl')
      
      #[ Select the M0 moment.
      distf = distf[:,0]
      if dI == 0:
        distf_main = distf
        time_main = time_distf
        continue

      main_interp = np.interp(time_distf, time_main, distf_main)
      
      hpl2a.append(ax2a.plot(time_distf, np.abs(distf - main_interp)/main_interp, color=defaultColors[dI], linestyle=lineStyles[dI], linewidth=2))
    
    ax2a.set_xlabel(r'Time ($s$)',fontsize=xyLabelFontSize, labelpad=+4)
    ax2a.set_ylabel(r'$abs(M_0 - M_{0,main}) / M_{0,main}$',fontsize=xyLabelFontSize, labelpad=0)
    ax2a.set_xlim( time_distf[0], time_distf[-1] )
    ax2a.legend([hpl2a[i][0] for i in range(len(hpl2a))], legendStrings_diff, fontsize=legendFontSize, frameon=False)
    ax2a.set_title(species, fontsize=titleFontSize)
    ax2a.set_yscale('log')
    setTickFontSize(ax2a,tickFontSize)
    plt.tight_layout()
    
    if outFigureFile:
      plt.savefig(outFileRoot+'m0_wrt_main_'+species+figureFileFormat)
    else:
      plt.show()
    plt.close(fig2a)

if plot_m2_wrt_main:
  #[ Plot the error normalized for different time steps.
  species_list = ['elc','ion']
  for sI in range(len(species_list)):
    species = species_list[sI]
    figProp2a = (12, 6)
    ax2aPos   = [0.09, 0.15, 0.87, 0.8]
    fig2a     = plt.figure(figsize=figProp2a)
    ax2a      = fig2a.add_axes(ax2aPos)
    
    hpl2a = list()
    
    for dI in range(len(cfl_dirs)):
      dataPath = dataDir + cfl_dirs[dI] + simName + '-'
      
      time_distf, distf = pgu.readDynVector(dataPath + species + '_integrated_moms.gkyl')
      
      #[ Select the M0 moment.
      distf = distf[:,2]
      if dI == 0:
        distf_main = distf
        time_main = time_distf
        continue

      main_interp = np.interp(time_distf, time_main, distf_main)
      
      hpl2a.append(ax2a.plot(time_distf, np.abs(distf - main_interp)/main_interp, color=defaultColors[dI], linestyle=lineStyles[dI], linewidth=2))
    
    ax2a.set_xlabel(r'Time ($s$)',fontsize=xyLabelFontSize, labelpad=+4)
    ax2a.set_ylabel(r'$abs(M_2 - M_{2,main}) / M_{2,main}$',fontsize=xyLabelFontSize, labelpad=0)
    ax2a.set_xlim( time_distf[0], time_distf[-1] )
    ax2a.legend([hpl2a[i][0] for i in range(len(hpl2a))], legendStrings_diff, fontsize=legendFontSize, frameon=False)
    ax2a.set_yscale('log')
    ax2a.set_title(species, fontsize=titleFontSize)
    setTickFontSize(ax2a,tickFontSize)
    plt.tight_layout()
    
    if outFigureFile:
      plt.savefig(outFileRoot+'m2_wrt_main_'+species+figureFileFormat)
    else:
      plt.show()
    plt.close(fig2a)

