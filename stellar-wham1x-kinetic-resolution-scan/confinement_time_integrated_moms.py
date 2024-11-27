2#[ ........................................................... ]#
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
filename = 'gk_wham'
frame_max_plus1 = 31
time_per_frame = 1e-6

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

resolution_list = np.array([32, 64, 96, 128, 192, 288])
# resolution_list = np.array([32])
confinement_time = np.zeros(len(resolution_list))
for ir in range(len(resolution_list)):
    resolution = resolution_list[ir]
    print('Resolution = ', resolution)
    frame_nums = np.arange(10, frame_max_plus1)
    densities = np.zeros(len(frame_nums))
    filename_read = dataDir + str(resolution) + '/misc/' + \
      filename + '-elc_integrated_moms.gkyl'
    print(filename_read)
    frame_data = pg.GData(filename_read)
    integrated_moms = frame_data.get_values()
    time = frame_data.get_grid()[0]
    M0_elc = np.array(integrated_moms[:,0])
    plt.plot(time, M0_elc, '-', label=str(resolution))


    def fit_function(t, a, b, c):
       return a * np.exp(b * t) + c
        # return m*t + b
    # find the index of time which is closest to 10e-6
    time_index = np.abs(time - 10e-6).argmin()
    fit_time = time[time_index:] - time[time_index]
    fit_densities = M0_elc[time_index:]
    print(np.shape(fit_time), np.shape(fit_densities))
    # Initial guess for the fit parameters.
    slope = (fit_densities[-1] - fit_densities[0])/(fit_time[-1] - fit_time[0])
    initial_guess = [fit_densities[0], slope / fit_densities[0], 0]
    print('Initial guess: ', initial_guess)
    popt, pcov = curve_fit(fit_function, fit_time, fit_densities, initial_guess)
    confinement_time[ir] = 1 / popt[1]
    # plt.plot(fit_time, fit_densities, 'o', label=str(resolution))
    # plt.plot(fit_time, fit_function(fit_time, *popt), '-', label=str(resolution))
plt.xlabel('Time')
plt.ylabel('Density')
plt.yscale('log')
plt.legend()


# Fit to a power law
def power_law(x, a, b):
    return a * x**b
popt, pcov = curve_fit(power_law, resolution_list, abs(confinement_time), p0=[1e-11, 3])
print('Power law fit: ', popt)
# Calculate R^2 from this fit
y_hat = power_law(resolution_list, *popt)
y_bar = np.mean(abs(confinement_time))
ss_tot = np.sum((abs(confinement_time) - y_bar)**2)
ss_res = np.sum((abs(confinement_time) - y_hat)**2)
r_squared = 1 - ss_res / ss_tot
print('R^2: ', r_squared)

plt.figure()
plt.plot(resolution_list, abs(confinement_time), 'bo')
# plt.plot(resolution_list, power_law(resolution_list, *popt), 'rx--')
plt.xlabel('Resolution')
plt.ylabel('Confinement time')
plt.yscale('log')
plt.xscale('log')
# plt.legend(['Data', 'Fit: $(%e) x^{%f}$: $R^2 = %f$' % (popt[0], popt[1], r_squared)])
plt.savefig(outDir + 'confinement time_vs_resolution.png')
plt.show()