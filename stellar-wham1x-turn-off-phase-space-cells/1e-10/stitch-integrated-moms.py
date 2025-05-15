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
origFile = '../stellar-wham1x-288z-run/misc/gk_wham-ion_integrated_moms.gkyl'

in_folders_orgnaized = 0 # 1 if it's in the folders, 0 if data is in present directory

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

fig, ax = plt.subplots(2,2, figsize=(12,20))
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
plt.show()
