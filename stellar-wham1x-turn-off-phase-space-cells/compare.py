
import numpy as np
import postgkyl as pg
import matplotlib.pyplot as plt
import sys
# sys.path.insert(0, '/Users/mfrancis/Documents/codebits/pets/gkeyll/postProcessingScripts/')
import pgkylUtil as pgu

dataDir = './'
outDir = './plots/'
outFileRoot = outDir + 'particle_balance_'

simName = 'gk_wham'

cfl_dirs = ['orig/', '1e-10/', '1e-8/', '1e-5/', '1e-3/']
legendStrings_error = [r'Main',r'Threshold = Maximum*1e-10',r'Threshold = Maximum*1e-8',r'Threshold = Maximum*1e-5',r'Threshold = Maximum*1e-3']
legendStrings_diff = [r'Threshold = Maximum*1e-10',r'Threshold = Maximum*1e-8',r'Threshold = Maximum*1e-5',r'Threshold = Maximum*1e-3']
species_list = ['ion']

plot_m0_wrt_main = True  #[ M0 moment of the distribution function with respect to the main one.
plot_m2_wrt_main = True  #[ M2 moment of the distribution function with respect to the main one.  

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

if plot_m0_wrt_main:
  #[ Plot the error normalized for different time steps.
  for sI in range(len(species_list)):
    species = species_list[sI]
    figProp2a = (12, 6)
    ax2aPos   = [0.09, 0.15, 0.87, 0.8]
    fig2a     = plt.figure(figsize=figProp2a)
    ax2a      = fig2a.add_axes(ax2aPos)
    
    hpl2a = list()
    
    for dI in range(len(cfl_dirs)):
      dataPath = dataDir + cfl_dirs[dI] + 'misc/' + simName + '-'
      
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
  for sI in range(len(species_list)):
    species = species_list[sI]
    figProp2a = (12, 6)
    ax2aPos   = [0.09, 0.15, 0.87, 0.8]
    fig2a     = plt.figure(figsize=figProp2a)
    ax2a      = fig2a.add_axes(ax2aPos)
    
    hpl2a = list()
    
    for dI in range(len(cfl_dirs)):
      dataPath = dataDir + cfl_dirs[dI] + 'misc/' + simName + '-'
      
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
