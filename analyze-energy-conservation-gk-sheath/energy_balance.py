#[ ........................................................... ]#
#[
#[ Check energy balance
#[
#[ Manaure Francisquez.
#[ Jan 2025
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
outFileRoot = outDir + 'energy_balance_'

simName = 'gk_sheath_1x2v_p1'

cfl_dirs = ['on_main/', 'thresh-1e-10/', 'thresh-1e-8/', 'thresh-1e-5/', 'thresh-1e-3/']
legendStrings_error = [r'Main',r'Threshold = Maximum*1e-10',r'Threshold = Maximum*1e-8',r'Threshold = Maximum*1e-5',r'Threshold = Maximum*1e-3']
legendStrings_wrt_main = [r'Threshold = Maximum*1e-10',r'Threshold = Maximum*1e-8',r'Threshold = Maximum*1e-5',r'Threshold = Maximum*1e-3']
plot_balance = True  #[ Balance of various terms.
plot_error   = True  #[ Conservation error.
plot_error_wrt_main = False # Take difference with respect to main. Ensure main is cfl_dirs[0]

outFigureFile    = True   #[ Output a figure file?.
figureFileFormat = '.png'  #[ Can be .png, .pdf, .ps, .eps, .svg.


#[ ............... End of user inputs (MAYBE) ..................... ]#

#[ Some RGB colors. These are MATLAB-like.
defaultBlue    = [0, 0.4470, 0.7410]
defaultOrange  = [0.8500, 0.3250, 0.0980]
defaultGreen   = [0.4660, 0.6740, 0.1880]
defaultPurple  = [0.4940, 0.1840, 0.5560]
defaultRed     = [0.6350, 0.0780, 0.1840]
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
  #[ Read the Hamiltonian moment of df/dt, the source and the particle fluxes.
  
  for dI in range(len(cfl_dirs)):
    dataPath = dataDir + cfl_dirs[dI] + simName + '-'
    run_name = cfl_dirs[dI][:-1]
    
    #[ Plot each contribution.
    figProp1a = (10, 6)
    ax1aPos   = [0.09, 0.15, 0.87, 0.8]
    fig1a     = plt.figure(figsize=figProp1a)
    ax1a      = fig1a.add_axes(ax1aPos)
    
    hpl1a = list()
    hpl1a.append(ax1a.plot([-1.0,1.0], [0.0,0.0], color='grey', linestyle=':', linewidth=1))
    
    time_field_dot, field_dot = pgu.readDynVector(dataPath + 'field_energy_dot.gkyl')
    
    for species in ['elc','ion']:
    
      time_fdot, fdot_s = pgu.readDynVector(dataPath + species + '_fdot_integrated_moms.gkyl')
      time_bflux_lo, bflux_lo_s = pgu.readDynVector(dataPath + species + '_bflux_xlower_integrated_HamiltonianMoments.gkyl')
      time_bflux_up, bflux_up_s = pgu.readDynVector(dataPath + species + '_bflux_xupper_integrated_HamiltonianMoments.gkyl')
      time_src, src_s = pgu.readDynVector(dataPath + species + '_source_integrated_moms.gkyl')
      
      #[ Select the Hamiltonian moment.
      fdot_s = fdot_s[:,2]
      bflux_lo_s = bflux_lo_s[:,2]
      bflux_up_s = bflux_up_s[:,2]
      src_s = src_s[:,2]
      
      src_s[0] = 0.0 #[ Set source=0 at t=0 since we don't have fdot and bflux then.
      
      bflux_tot_s = bflux_lo_s+bflux_up_s #[ Total boundary flux loss.
    
      if species == 'elc':
        fdot      = fdot_s     
        src       = src_s      
        bflux_tot = bflux_tot_s
      else:
        fdot      = fdot      + fdot_s     
        src       = src       + src_s      
        bflux_tot = bflux_tot + bflux_tot_s
    
    mom_err = src - bflux_tot - (fdot - field_dot) #[ Error.
      
    hpl1a.append(ax1a.plot(time_fdot, fdot, color=defaultColors[0], linestyle=lineStyles[0], linewidth=2, marker=markers[0]))
    hpl1a.append(ax1a.plot(time_bflux_lo, -bflux_tot, color=defaultColors[1], linestyle=lineStyles[1], linewidth=2, marker=markers[1]))
    hpl1a.append(ax1a.plot(time_src, src, color=defaultColors[2], linestyle=lineStyles[2], linewidth=2, marker=markers[2]))
    hpl1a.append(ax1a.plot(time_field_dot, field_dot, color=defaultColors[4], linestyle=':', linewidth=2, marker='+',markevery=8))
    hpl1a.append(ax1a.plot(time_fdot, mom_err, color=defaultColors[3], linestyle=lineStyles[3], linewidth=2, marker=markers[3]))
    
    ax1a.set_xlabel(r'Time ($s$)',fontsize=xyLabelFontSize, labelpad=+4)
    ax1a.set_ylabel(r'Hamiltonian  moment',fontsize=xyLabelFontSize, labelpad=0)
    ax1a.set_xlim( time_fdot[0], time_fdot[-1] )
    ax1a.set_title(run_name, fontsize=titleFontSize)
    legendStrings = [r'$\dot{f}$',
    r'$-\int_{\partial \Omega}\mathrm{d}\mathbf{S}\cdot\mathbf{\dot{R}}f$',
    r'$\mathcal{S}$',
    r'$\dot{\phi}$',
    r'$E_{\dot{\mathcal{E}}}=\mathcal{S}-\int_{\partial \Omega}\mathrm{d}\mathbf{S}\cdot\mathbf{\dot{R}}f-(\dot{f}-\dot{\phi})$',
    ]
    ax1a.legend([hpl1a[i][0] for i in range(1,len(hpl1a))], legendStrings, fontsize=legendFontSize, frameon=False, loc='lower right')
    setTickFontSize(ax1a,tickFontSize)
    plt.tight_layout()

    
    if outFigureFile:
      plt.savefig(outFileRoot+run_name+'_energy'+figureFileFormat)
    else:
      plt.show()

if plot_error:
  #[ Plot the error normalized for different time steps.
  eV = 1.602e-19
  me = 0.91e-30
  Lz = 4.0
  n0 = 7.0e18
  N0 = 9e17 #[ Reference number of particles
  Te = 40.0*eV
  vte = np.sqrt(Te/me)
  tau = 0.5*Lz/vte #[ Transit time.
  
  figProp2a = (10, 6)
  ax2aPos   = [0.09, 0.15, 0.87, 0.8]
  fig2a     = plt.figure(figsize=figProp2a)
  ax2a      = fig2a.add_axes(ax2aPos)
  
  hpl2a = list()
  hpl2a.append(ax2a.plot([-1.0,1.0], [0.0,0.0], color='grey', linestyle=':', linewidth=1))
  ylabelString = ""
  
  for dI in range(len(cfl_dirs)):
    dataPath = dataDir + cfl_dirs[dI] + simName + '-'
    
    time_field, field = pgu.readDynVector(dataPath + 'field_energy.gkyl')
    time_field_dot, field_dot = pgu.readDynVector(dataPath + 'field_energy_dot.gkyl')
    field = field[1:]
    field_dot = field_dot[1:]
  
    for species in ['elc','ion']:
      time_fdot, fdot_s = pgu.readDynVector(dataPath + species + '_fdot_integrated_moms.gkyl')
      time_bflux_lo, bflux_lo_s = pgu.readDynVector(dataPath + species + '_bflux_xlower_integrated_HamiltonianMoments.gkyl')
      time_bflux_up, bflux_up_s = pgu.readDynVector(dataPath + species + '_bflux_xupper_integrated_HamiltonianMoments.gkyl')
      time_src, src_s = pgu.readDynVector(dataPath + species + '_source_integrated_moms.gkyl')
  
      time_distf, distf_s = pgu.readDynVector(dataPath + species + '_integrated_moms.gkyl')
    
      time_dt, dt = pgu.readDynVector(dataPath + 'dt.gkyl')
      
      #[ Select the Hamiltonian moment and remove the t=0 data point.
      fdot_s = fdot_s[1:,2]
      bflux_lo_s = 0.*bflux_lo_s[1:,2]
      bflux_up_s = 0.*bflux_up_s[1:,2]
      src_s = src_s[1:,2]
      distf_s = distf_s[1:,2]
      
      bflux_tot_s = bflux_lo_s+bflux_up_s #[ Total boundary flux loss.
  
      if species == 'elc':
        fdot      = fdot_s     
        src       = src_s      
        bflux_tot = bflux_tot_s
        distf     = distf_s      
      else:
        fdot      = fdot      + fdot_s     
        src       = src       + src_s      
        bflux_tot = bflux_tot + bflux_tot_s
        distf     = distf     + distf_s      
  
    mom_err = src - bflux_tot - (fdot - field_dot) #[ Error.
  
  #  mom_err_norm = mom_err
  #  ylabelString = r'$E_{\dot{\mathcal{E}}}$'
    mom_err_norm = mom_err*dt/(distf-field)
    ylabelString = r'$E_{\dot{\mathcal{E}}}~\Delta t/\mathcal{E}$'
  #  mom_err_norm = mom_err*tau/N0
  #  ylabelString = r'$E_{\dot{\mathcal{E}}}~\tau/\mathcal{E}_0$'
  
    mom_err_norm_mean = np.mean(mom_err_norm[np.size(mom_err_norm)//2:])
    # if (dI == 0):
    #   print(mom_err_norm_mean)
    # else:
    #   print(mom_err_norm_mean, mom_err_norm_mean_prev/mom_err_norm_mean)
    mom_err_norm_mean_prev = mom_err_norm_mean
    
    hpl2a.append(ax2a.plot(time_dt, mom_err_norm, color=defaultColors[dI], linestyle=lineStyles[dI], linewidth=2))
  

  ax2a.set_xlabel(r'Time ($s$)',fontsize=xyLabelFontSize, labelpad=+4)
  ax2a.set_ylabel(ylabelString,fontsize=xyLabelFontSize, labelpad=0)
  ax2a.set_xlim( time_fdot[0], time_fdot[-1] )
  ax2a.legend([hpl2a[i][0] for i in range(1,len(hpl2a))], legendStrings_error, fontsize=legendFontSize, frameon=False)
  setTickFontSize(ax2a,tickFontSize)
  
  if outFigureFile:
    plt.savefig(outFileRoot+'Energy_change'+figureFileFormat)
  else:
    plt.show()
  
if plot_error_wrt_main:
  #[ Plot the error normalized for different time steps.
  
  eV = 1.602e-19
  me = 0.91e-30
  Lz = 4.0
  n0 = 7.0e18
  N0 = 9e17 #[ Reference number of particles
  Te = 40.0*eV
  vte = np.sqrt(Te/me)
  tau = 0.5*Lz/vte #[ Transit time.
  
  figProp2a = (10, 6)
  ax2aPos   = [0.09, 0.15, 0.87, 0.8]
  fig2a     = plt.figure(figsize=figProp2a)
  ax2a      = fig2a.add_axes(ax2aPos)
  
  hpl2a = list()
  hpl2a.append(ax2a.plot([-1.0,1.0], [0.0,0.0], color='grey', linestyle=':', linewidth=1))
  ylabelString = ""
  
  for dI in range(len(cfl_dirs)):
    dataPath = dataDir + cfl_dirs[dI] + simName + '-'
    
    time_field, field = pgu.readDynVector(dataPath + 'field_energy.gkyl')
    time_field_dot, field_dot = pgu.readDynVector(dataPath + 'field_energy_dot.gkyl')
    field = field[1:]
    field_dot = field_dot[1:]
  
    for species in ['elc','ion']:
      time_fdot, fdot_s = pgu.readDynVector(dataPath + species + '_fdot_integrated_moms.gkyl')
      time_bflux_lo, bflux_lo_s = pgu.readDynVector(dataPath + species + '_bflux_xlower_integrated_HamiltonianMoments.gkyl')
      time_bflux_up, bflux_up_s = pgu.readDynVector(dataPath + species + '_bflux_xupper_integrated_HamiltonianMoments.gkyl')
      time_src, src_s = pgu.readDynVector(dataPath + species + '_source_integrated_moms.gkyl')
  
      time_distf, distf_s = pgu.readDynVector(dataPath + species + '_integrated_moms.gkyl')
    
      time_dt, dt = pgu.readDynVector(dataPath + 'dt.gkyl')
      
      #[ Select the Hamiltonian moment and remove the t=0 data point.
      fdot_s = fdot_s[1:,2]
      bflux_lo_s = 0.*bflux_lo_s[1:,2]
      bflux_up_s = 0.*bflux_up_s[1:,2]
      src_s = src_s[1:,2]
      distf_s = distf_s[1:,2]
      
      bflux_tot_s = bflux_lo_s+bflux_up_s #[ Total boundary flux loss.
  
      if species == 'elc':
        fdot      = fdot_s     
        src       = src_s      
        bflux_tot = bflux_tot_s
        distf     = distf_s      
      else:
        fdot      = fdot      + fdot_s     
        src       = src       + src_s      
        bflux_tot = bflux_tot + bflux_tot_s
        distf     = distf     + distf_s      
  
    mom_err = src - bflux_tot - (fdot - field_dot) #[ Error.
  
  #  mom_err_norm = mom_err
  #  ylabelString = r'$E_{\dot{\mathcal{E}}}$'
    mom_err_norm = mom_err*dt/(distf-field)
    ylabelString = r'$\Delta E_{\dot{\mathcal{E}}}~\Delta t/\mathcal{E}$ versus main'
  #  mom_err_norm = mom_err*tau/N0
  #  ylabelString = r'$E_{\dot{\mathcal{E}}}~\tau/\mathcal{E}_0$'
  
    mom_err_norm_mean = np.mean(mom_err_norm[np.size(mom_err_norm)//2:])
    # if (dI == 0):
    #   print(mom_err_norm_mean)
    # else:
    #   print(mom_err_norm_mean, mom_err_norm_mean_prev/mom_err_norm_mean)
    mom_err_norm_mean_prev = mom_err_norm_mean
    if (dI==0):
      mom_err_norm_main = mom_err_norm
      time_main = time_dt
      continue

    err_main_interp = np.interp(time_dt, time_main, mom_err_norm_main)
    
    hpl2a.append(ax2a.plot(time_dt, np.abs(mom_err_norm - err_main_interp), color=defaultColors[dI], linestyle=lineStyles[dI], linewidth=2))
  
  ax2a.set_xlabel(r'Time ($s$)',fontsize=xyLabelFontSize, labelpad=+4)
  ax2a.set_ylabel(ylabelString,fontsize=xyLabelFontSize, labelpad=0)
  ax2a.set_xlim( time_fdot[0], time_fdot[-1] )
  ax2a.legend([hpl2a[i][0] for i in range(1,len(hpl2a))], legendStrings_wrt_main, fontsize=legendFontSize, frameon=False)
  setTickFontSize(ax2a,tickFontSize)
  ax2a.set_yscale('log')
  
  if outFigureFile:
    plt.savefig(outFileRoot+'Energy_vs_main'+figureFileFormat)
  else:
    plt.show()
  
