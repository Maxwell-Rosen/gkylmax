
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
dataDir = './288-cells'

frame_max_plus1 = 107
time_per_frame = 1e-6

plot_potential_trace = 0
plot_bimax_moms = 0
plot_integrate_positivity = 1

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
