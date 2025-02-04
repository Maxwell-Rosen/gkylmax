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

Nz_arr = ['128', '192', '288']
nonuniform_frac = ['0.000', '0.333', '0.666', '0.999']

plot_potential_trace_combined = 1
plot_potential_trace_seperate = 1
plot_dndt_trace = 1
plot_intM0 = 1

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

linecolors_10 = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown']
linestyles_10 = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--']

linecolors_12 = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown', 'pink', 'gray']
linestyles_12 = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.', ':']


#   #................................................................................#

def plot_potential_trace_combined():
  print("Plotting potential trace")
  plt.figure()
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]

      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'
      fieldname = filebase+'Field/gk_wham-field'
      bmagname = filebase+'Geometry/gk_wham-bmag.gkyl'

      pgData_bmag = pg.GData(bmagname)
      pgInterp_bmag = pg.GInterpModal(pgData_bmag, polyOrder, 'ms')
      x_bmag, dataOut_bmag = pgInterp_bmag.interpolate()

      bmag_shape = dataOut_bmag.shape
      midpoint = int(bmag_shape[0]/2)
      upperhalf = dataOut_bmag[midpoint:]
      peak = np.argmax(upperhalf)
      peak_idx = midpoint+peak

      def loadphi(frame_number, filename):
        filename_phi = str(filebase+'Field/gk_wham-field_'+str(frame_number)+'.gkyl')
        pgData_phi = pg.GData(filename_phi)
        pgInterp_phi = pg.GInterpModal(pgData_phi, polyOrder, 'ms')
        x_phi, dataOut_phi = pgInterp_phi.interpolate()
        return dataOut_phi
      
      def get_temp():
        return Te0

      potential = np.zeros(frame_max_plus1)
      Temp = np.zeros(frame_max_plus1)
      for i in range(frame_max_plus1):
        dataOut_phi = loadphi(i, mc2pFilename)
        Temp[i] = get_temp()
        midphi = dataOut_phi[midpoint]
        phi_peak = dataOut_phi[peak_idx]
        potential[i] = (midphi[0] - phi_peak[0]) / Temp[i]

      plt.plot(np.arange(frame_max_plus1)*time_per_frame, potential, color = linecolors_12[linenumber], linestyle = linestyles_12[linenumber], \
          label = 'Nz = '+Nz+', nonuniform frac = '+nonuniform_frac_str)
      linenumber += 1
  
  
  plt.xlabel('Time, seconds')
  plt.ylabel('Potential difference, $e \phi / T_e(\psi_{min},z=0)$')
  plt.title('Potential difference between midplane and peak magnetic field')
  plt.legend()
  plt.savefig(outDir+'potential_trace_combined'+figureFileFormat)
  plt.close()

def plot_potential_trace_seperate():
  print("Plotting potential trace seperate")
  plt.figure(figsize=(20, 14))
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]

      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'
      fieldname = filebase+'Field/gk_wham-field'
      bmagname = filebase+'Geometry/gk_wham-bmag.gkyl'

      pgData_bmag = pg.GData(bmagname)
      pgInterp_bmag = pg.GInterpModal(pgData_bmag, polyOrder, 'ms')
      x_bmag, dataOut_bmag = pgInterp_bmag.interpolate()

      bmag_shape = dataOut_bmag.shape
      midpoint = int(bmag_shape[0]/2)
      upperhalf = dataOut_bmag[midpoint:]
      peak = np.argmax(upperhalf)
      peak_idx = midpoint+peak

      def loadphi(frame_number, filename):
        filename_phi = str(filebase+'Field/gk_wham-field_'+str(frame_number)+'.gkyl')
        pgData_phi = pg.GData(filename_phi)
        pgInterp_phi = pg.GInterpModal(pgData_phi, polyOrder, 'ms')
        x_phi, dataOut_phi = pgInterp_phi.interpolate()
        return dataOut_phi
      
      def get_temp():
        return Te0

      potential = np.zeros(frame_max_plus1)
      Temp = np.zeros(frame_max_plus1)
      for i in range(frame_max_plus1):
        dataOut_phi = loadphi(i, mc2pFilename)
        Temp[i] = get_temp()
        midphi = dataOut_phi[midpoint]
        phi_peak = dataOut_phi[peak_idx]
        potential[i] = (midphi[0] - phi_peak[0]) / Temp[i]

      plt.subplot(3, 4, linenumber+1)
      plt.plot(np.arange(frame_max_plus1)*time_per_frame, potential, color = linecolors_12[linenumber], linestyle = linestyles_12[linenumber], \
          label = 'Nz = '+Nz+', nonuniform frac = '+nonuniform_frac_str)
      linenumber += 1
      plt.legend()
      # plt.yscale('log')
      plt.title("N_z = "+Nz+", nonuniform frac = "+nonuniform_frac_str)
      plt.xlabel('Time, seconds')
      plt.ylabel('Potential difference, $e \phi / T_e(\psi_{min},z=0)$')
  
  plt.legend()
  plt.tight_layout()
  plt.savefig(outDir+'potential_trace_seperate'+figureFileFormat)
  plt.close()

def plot_dndt_trace():
  print("Plotting dndt trace")
  plt.figure( figsize=(20, 14) )
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]

      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'
      def loadM0(frame_number):
        filename_M0 = str(filebase+'M/gk_wham-ion_M0_'+str(frame_number)+'.gkyl')
        pgData_M0 = pg.GData(filename_M0)
        pgInterp_M0 = pg.GInterpModal(pgData_M0, polyOrder, 'ms')
        x_M0, dataOut_M0 = pgInterp_M0.interpolate()
        return x_M0[0], dataOut_M0[:,0]
      
      intM0 = np.zeros(frame_max_plus1)
      for i in range(frame_max_plus1):
        x_M0, dataOut_M0 = loadM0(i)
        intM0[i] = np.trapz(dataOut_M0, x_M0[0:-1])

      start_ix = int(len(intM0)/2)
      dndt = -np.gradient(intM0[start_ix:], time_per_frame)
      
      plt.subplot(3, 4, linenumber+1)
      plt.plot(np.arange(frame_max_plus1-start_ix)*time_per_frame, dndt, color = linecolors_12[linenumber], linestyle = linestyles_12[linenumber], \
          label = 'Nz = '+Nz+', nonuniform frac = '+nonuniform_frac_str)
      linenumber += 1
      plt.legend()
      # plt.yscale('log')
      plt.title("N_z = "+Nz+", nonuniform frac = "+nonuniform_frac_str)
      plt.xlabel('Time, seconds')
      plt.ylabel('dN/dt')

  plt.tight_layout()
  plt.savefig(outDir+'dndt_trace'+figureFileFormat)

def plot_intM0():
  print("Plotting integral of M0 trace")
  plt.figure( figsize=(20, 14) )
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]

      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'
      def loadM0(frame_number):
        filename_M0 = str(filebase+'M/gk_wham-ion_M0_'+str(frame_number)+'.gkyl')
        pgData_M0 = pg.GData(filename_M0)
        pgInterp_M0 = pg.GInterpModal(pgData_M0, polyOrder, 'ms')
        x_M0, dataOut_M0 = pgInterp_M0.interpolate()
        return x_M0[0], dataOut_M0[:,0]
      
      intM0 = np.zeros(frame_max_plus1)
      for i in range(frame_max_plus1):
        x_M0, dataOut_M0 = loadM0(i)
        intM0[i] = np.trapz(dataOut_M0, x_M0[0:-1])
      
      plt.subplot(3, 4, linenumber+1)
      plt.plot(np.arange(frame_max_plus1)*time_per_frame, intM0, color = linecolors_12[linenumber], linestyle = linestyles_12[linenumber], \
          label = 'Nz = '+Nz+', nonuniform frac = '+nonuniform_frac_str)
      linenumber += 1
      plt.legend()
      # plt.yscale('log')
      plt.title("N_z = "+Nz+", nonuniform frac = "+nonuniform_frac_str)
      plt.xlabel('Time, seconds')
      plt.ylabel('$\int M_0 dx$')

  plt.tight_layout()
  plt.savefig(outDir+'intM0_trace'+figureFileFormat)


## Main function

def main():
  processes = []
  if plot_potential_trace_combined:
    p = multiprocessing.Process(target=plot_potential_trace_combined)
    processes.append(p)
    p.start()
  if plot_potential_trace_seperate:
    p = multiprocessing.Process(target=plot_potential_trace_seperate)
    processes.append(p)
    p.start()
  if plot_dndt_trace:
    p = multiprocessing.Process(target=plot_dndt_trace)
    processes.append(p)
    p.start()
  if plot_intM0:
    p = multiprocessing.Process(target=plot_intM0)
    processes.append(p)
    p.start()

  for p in processes:
    p.join()

if __name__ == "__main__":
  main()
