
import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from matplotlib.colors import LogNorm
import multiprocessing
from scipy.integrate import cumulative_trapezoid as cumtrapz
import imageio.v2 as imageio
import multiprocessing

# dataDir = '/home/mr1884/scratch/Link to scratch_traverse/gkylmax/traverse-wham1x-compare_unif_vs_nonunif/outputs/'
# dataDir = './data-hires-lorad/'
dataDir = './'
mc2pFolder = './mc2p_fix/'
mc2pFilename = 'gk_mirror'
mc2pUniformFolder = './mc2p_uniform_fix/'
mc2pUniformFilename = 'gk_mirror'
numericFolder = './mc2p/'
numericFilename = 'gk_mirror'
frame_max_plus1 = 301
time_per_frame = 1e-6 / 3.

Nz_arr = ['32', '64', '128', '192', '288']
nonuniform_frac = ['0.000', '0.333', '0.666', '0.999']

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

#[ Thermal speeds.
vti = np.sqrt(Ti0/mi)
vte = np.sqrt(Te0/me)
mui0 = 0.5*mi*(vti**2)/B_p
c_s = np.sqrt(Te0/mi)

linecolors = [
    'b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown',
    'pink', 'lime', 'teal', 'navy', 'maroon', 'gold', 'orchid', 'gray', 'olive', 'cyan',
    'deepskyblue', 'indigo', 'coral', 'darkgreen', 'crimson', 'sienna', 'plum', 'tomato', 'turquoise', 'slateblue',
    'darkorange', 'mediumseagreen', 'hotpink', 'peru', 'dodgerblue', 'darkred', 'mediumvioletred', 'springgreen', 'steelblue', 'khaki',
    'mediumorchid', 'firebrick', 'seagreen', 'darkviolet', 'lightcoral', 'royalblue', 'chartreuse', 'wheat', 'darkcyan', 'salmon'
]
linestyles = np.tile(['-', '--', '-.', ':'], 15)

#   #................................................................................#

filename_M0 = str('gk_wham-ion_BiMaxwellianMoments_300.gkyl')
pgData_M0 = pg.GData(filename_M0)
pgInterp_M0 = pg.GInterpModal(pgData_M0, polyOrder, 'ms')
x, dens_nunif = pgInterp_M0.interpolate(0)

# filename_M0 = str('gk_wham_1x2v_p1-ion_BiMaxwellianMoments_0.gkyl')
# pgData_M0 = pg.GData(filename_M0)
# pgInterp_M0 = pg.GInterpModal(pgData_M0, polyOrder, 'ms')
# x, dens_unif = pgInterp_M0.interpolate(0)

data = pg.GData(str('gk_wham-mc2nu_pos.gkyl'))
interp = pg.GInterpModal(data, 1, 'ms')
nodes_Z_nunif = interp.interpolate(2)[1]
nodes_Z_nunif = np.squeeze(nodes_Z_nunif)

# data = pg.GData(str('gk_wham_1x2v_p1-mc2nu_pos.gkyl'))
# interp = pg.GInterpModal(data, 1, 'ms')
# nodes_Z_unif = interp.interpolate(2)[1]
# nodes_Z_unif = np.squeeze(nodes_Z_unif)

plt.figure()

print(np.shape(dens_nunif))
# plt.plot(nodes_Z_unif, dens_unif, label='Uniform', color='blue')
plt.plot(nodes_Z_nunif, dens_nunif, label='Nonuniform', color='red')
plt.legend()
plt.title('Density')
plt.xlabel('Z, m')
plt.ylabel('Density, $m^{-3}$')
plt.tight_layout()
plt.savefig('Density comparison.png')