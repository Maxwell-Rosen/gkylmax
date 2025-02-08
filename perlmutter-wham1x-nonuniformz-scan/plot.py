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
frame_max_plus1 = 301#301
time_per_frame = 1e-6 / 3.

# Nz_arr = ['32', '64', '96', '128', '192', '288']
# nonuniform_frac = ['0.000', '0.333', '0.666', '0.999']
frame_arr = [100, 300]
Nz_arr = ['288', '48', '64', '96']
nonuniform_frac = ['0.000', '0.100', '0.200', '0.300', '0.400', '0.500', '0.600', '0.666', '0.700', '0.800', '0.900', '0.999']

outDir = './python-plots/'

figureFileFormat = '.png'    #[ Can be .png, .pdf, .ps, .eps, .svg.

figsize_dim = (20, 14)
# figsize_dim=(10, 7)
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
linecolors = np.tile(linecolors, [len(frame_arr),1])
linecolors = np.reshape(linecolors, [len(frame_arr)*len(linecolors[0])], order='F')

linestyles = np.tile(['-', '--', '-.', ':'], 15)

#   #................................................................................#

def skip_sims(Nz, nonuniform_frac):
  # if Nz == '32' and nonuniform_frac == '0.000':
  #   return True
  # if Nz == '288' and nonuniform_frac != '0.000':
  #   return True
  # if Nz == '48' and nonuniform_frac =='0.000':
  #   return True
  # if Nz == '64' and (nonuniform_frac =='0.000' or \
  #             nonuniform_frac =='0.100' or \
  #             nonuniform_frac =='0.200' or \
  #             nonuniform_frac =='0.300' or \
  #             nonuniform_frac =='0.400'): \
  #   return True
  return False

def plot_only_sims(Nz, nonuniform_frac):
  if Nz == '288' and (nonuniform_frac == '0.000' or \
                      nonuniform_frac == '0.666'):
    return True
  if Nz == '64' and nonuniform_frac == '0.666':
    return True
  if Nz == '96' and nonuniform_frac == '0.666':
  # if Nz == '48' and (nonuniform_frac == '0.100' or \
  #                     nonuniform_frac == '0.200' or \
  #                     nonuniform_frac == '0.300' or \
  #                     nonuniform_frac == '0.400' or \
  #                     nonuniform_frac == '0.500' or \
  #                     nonuniform_frac == '0.600' or \
  #                     nonuniform_frac == '0.700' or \
  #                     nonuniform_frac == '0.800' or \
  #                     nonuniform_frac == '0.900' or \
  #                     nonuniform_frac == '0.999'):       
  #   return True
  # if Nz == '64' and (nonuniform_frac == '0.500' or \
  #                     nonuniform_frac == '0.600' or \
  #                     nonuniform_frac == '0.666'):
    return True  
  # if Nz == '64' and (nonuniform_frac == '0.666'):
  #   return True
  return False

def load_mom(frame_number, filebase):
  filename_M0 = str(filebase+'/BiMaxwellianMoments/gk_wham-ion_BiMaxwellianMoments_'\
                    +str(frame_number)+'.gkyl')
  pgData_M0 = pg.GData(filename_M0)
  pgInterp_M0 = pg.GInterpModal(pgData_M0, polyOrder, 'ms')
  x, dens = pgInterp_M0.interpolate(0)
  x, u = pgInterp_M0.interpolate(1)
  x, Tpar = pgInterp_M0.interpolate(2)
  x, Tperp = pgInterp_M0.interpolate(3)
  return x[0], dens[:,0], u[:,0], Tpar[:,0]*mi/eV/1000, Tperp[:,0]*mi/eV/1000

def loadphi(frame_number, filebase):
  filename_phi = str(filebase+'Field/gk_wham-field_'+str(frame_number)+'.gkyl')
  pgData_phi = pg.GData(filename_phi)
  pgInterp_phi = pg.GInterpModal(pgData_phi, polyOrder, 'ms')
  x_phi, dataOut_phi = pgInterp_phi.interpolate()
  return dataOut_phi

def loadbmag(filebase):
  bmagname = filebase+'Geometry/gk_wham-bmag.gkyl'
  pgData_bmag = pg.GData(bmagname)
  pgInterp_bmag = pg.GInterpModal(pgData_bmag, polyOrder, 'ms')
  x_bmag, dataOut_bmag = pgInterp_bmag.interpolate()
  return x_bmag, dataOut_bmag

def get_temp():
  return Te0

def loadM0(frame_number, filebase):
  filename_M0 = str(filebase+'M/gk_wham-ion_M0_'+str(frame_number)+'.gkyl')
  pgData_M0 = pg.GData(filename_M0)
  pgInterp_M0 = pg.GInterpModal(pgData_M0, polyOrder, 'ms')
  x_M0, dataOut_M0 = pgInterp_M0.interpolate()
  return x_M0[0], dataOut_M0[:,0]

def loadM0_integrated(filebase):
  filename_M0 = str(filebase+'misc/gk_wham-ion_integrated_moms.gkyl')
  pgData_M0 = pg.GData(filename_M0)
  M = pgData_M0.get_values()
  M0 = np.array(M[:,0])
  time = np.squeeze(np.array(pgData_M0.get_grid()))
  return time, M0

def expand_1D_array(original_array):
  new_length = 2 * len(original_array) - 1
  new_array = np.zeros(new_length)
  new_array[0] = original_array[0]
  for i in range(1, len(original_array)):
      new_array[2*i - 1] = (original_array[i - 1] + original_array[i]) / 2
      new_array[2*i] = original_array[i]
  return new_array

def loadZ(filebase):
    data = pg.GData(str(filebase+'/Geometry/gk_wham-nodes.gkyl'))
    vals = data.get_values()
    nodes_Z = vals[:,1]
    nodes_Z = expand_1D_array(nodes_Z)
    nodes_Z = nodes_Z[1:]
    return nodes_Z

def plot_potential_trace_combined():
  print("Plotting potential trace combined")
  plt.figure( figsize = figsize_dim)
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]
      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'

      skip = skip_sims(Nz, nonuniform_frac_str)
      noskip = plot_only_sims(Nz, nonuniform_frac_str)
      if skip or not noskip:
        continue

      x_bmag, dataOut_bmag = loadbmag(filebase)
      bmag_shape = dataOut_bmag.shape
      midpoint = int(bmag_shape[0]/2)
      upperhalf = dataOut_bmag[midpoint:]
      peak = np.argmax(upperhalf)
      peak_idx = midpoint+peak

      potential = np.zeros(frame_max_plus1)
      Temp = np.zeros(frame_max_plus1)
      for i in range(frame_max_plus1):
        dataOut_phi = loadphi(i, filebase)
        Temp[i] = get_temp()
        midphi = dataOut_phi[midpoint]
        phi_peak = dataOut_phi[peak_idx]
        potential[i] = (midphi[0] - phi_peak[0]) / Temp[i]*eV

      plt.plot(np.arange(frame_max_plus1)*time_per_frame, potential, \
               color = linecolors[linenumber], linestyle = linestyles[linenumber], \
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

      skip = skip_sims(Nz, nonuniform_frac_str)
      if skip:
        linenumber += 1
        continue

      x_bmag, dataOut_bmag = loadbmag(filebase)
      bmag_shape = dataOut_bmag.shape
      midpoint = int(bmag_shape[0]/2)
      upperhalf = dataOut_bmag[midpoint:]
      peak = np.argmax(upperhalf)
      peak_idx = midpoint+peak

      potential = np.zeros(frame_max_plus1)
      Temp = np.zeros(frame_max_plus1)
      for i in range(frame_max_plus1):
        dataOut_phi = loadphi(i, filebase)
        Temp[i] = get_temp()
        midphi = dataOut_phi[midpoint]
        phi_peak = dataOut_phi[peak_idx]
        potential[i] = (midphi[0] - phi_peak[0]) / Temp[i]*eV

      plt.subplot(len(Nz_arr), len(nonuniform_frac), linenumber+1)
      plt.plot(np.arange(frame_max_plus1)*time_per_frame, potential, \
               color = linecolors[linenumber], linestyle = linestyles[linenumber])
      linenumber += 1
      plt.title("N_z = "+Nz+", nonuniform frac = "+nonuniform_frac_str)
      plt.xlabel('Time, seconds')
      plt.ylabel('Potential difference, $e \phi / T_e(\psi_{min},z=0)$')
  
  plt.tight_layout()
  plt.savefig(outDir+'potential_trace_seperate'+figureFileFormat)
  plt.close()

def plot_dndt_trace():
  print("Plotting dndt trace")
  plt.figure( figsize = figsize_dim)
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]
      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'
      
      skip = skip_sims(Nz, nonuniform_frac_str)
      if skip:
        linenumber += 1
        continue

      intM0 = np.zeros(frame_max_plus1)
      nodes_Z = loadZ(filebase)
      for i in range(frame_max_plus1):
        x_M0, dataOut_M0 = loadM0(i, filebase)
        intM0[i] = np.trapz(dataOut_M0, nodes_Z)
      start_ix = int(len(intM0)/2)
      t = np.arange(frame_max_plus1)*time_per_frame
      dndt = np.zeros(frame_max_plus1-10)
      for i in range(len(intM0)-10):
        dndt[i] = -(intM0[i+10] - intM0[i]) / (t[i+10] - t[i])
      

      # t, M0 = loadM0_integrated(filebase)
      # start_ix = int(len(M0)/2)
      # dt = t[1] - t[0]
      # #  Take a gradient using finite difference
      # dndt = np.zeros(len(M0))
      # for i in range(len(M0)-10):
      #   dndt[i] = -(M0[i+10] - M0[i]) / (t[i+10] - t[i])

      plt.subplot(len(Nz_arr), len(nonuniform_frac), linenumber+1)
      # plt.plot(t[start_ix:-10], dndt[start_ix:-10],\
      #           color = linecolors[linenumber], linestyle = linestyles[linenumber])
      plt.plot(t[start_ix:-10], dndt[start_ix:],\
                color = linecolors[linenumber], linestyle = linestyles[linenumber])
      linenumber += 1
      plt.title("N_z = "+Nz+", nonuniform frac = "+nonuniform_frac_str)
      plt.xlabel('Time, seconds')
      plt.ylabel('dN/dt')

  plt.tight_layout()
  plt.savefig(outDir+'dndt_trace'+figureFileFormat)

def plot_dndt_trace_combined():
  print("Plotting dndt trace combined")
  plt.figure( figsize = figsize_dim)
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]
      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'
      
      skip = skip_sims(Nz, nonuniform_frac_str)
      noskip = plot_only_sims(Nz, nonuniform_frac_str)
      if skip or not noskip:
        continue

      intM0 = np.zeros(frame_max_plus1)
      nodes_Z = loadZ(filebase)
      for i in range(frame_max_plus1):
        x_M0, dataOut_M0 = loadM0(i, filebase)
        intM0[i] = np.trapz(dataOut_M0, nodes_Z)
      start_ix = int(len(intM0)/2)
      t = np.arange(frame_max_plus1)*time_per_frame
      dndt = np.zeros(frame_max_plus1-10)
      for i in range(len(intM0)-10):
        dndt[i] = -(intM0[i+10] - intM0[i]) / (t[i+10] - t[i])
      

      # t, M0 = loadM0_integrated(filebase)
      # start_ix = int(len(M0)/2)
      # dt = t[1] - t[0]
      # #  Take a gradient using finite difference
      # dndt = np.zeros(len(M0))
      # for i in range(len(M0)-10):
      #   dndt[i] = -(M0[i+10] - M0[i]) / (t[i+10] - t[i])

      # plt.plot(t[start_ix:-10], dndt[start_ix:-10],\
      #           color = linecolors[linenumber], linestyle = linestyles[linenumber])
      plt.plot(t[start_ix:-10], dndt[start_ix:],\
                color = linecolors[linenumber], linestyle = linestyles[linenumber],
                label = 'Nz = '+Nz+', nonuniform frac = '+nonuniform_frac_str)
      linenumber += 1

  plt.xlabel('Time, seconds')
  plt.ylabel('dN/dt')
  plt.title('dN/dt')
  plt.legend()
  plt.tight_layout()
  plt.savefig(outDir+'dndt_trace_combined'+figureFileFormat)

def plot_intM0():
  print("Plotting integral of M0 trace")
  plt.figure( figsize = figsize_dim)
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]
      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'
      skip = skip_sims(Nz, nonuniform_frac_str)
      if skip:
        linenumber += 1
        continue

      intM0 = np.zeros(frame_max_plus1)
      nodes_Z = loadZ(filebase)
      for i in range(frame_max_plus1):
        x_M0, dataOut_M0 = loadM0(i, filebase)
        intM0[i] = np.trapz(dataOut_M0, nodes_Z)

      # t, M0 = loadM0_integrated(filebase)
      
      plt.subplot(len(Nz_arr), len(nonuniform_frac), linenumber+1)
      # plt.plot(t, M0, \
      #          color = linecolors[linenumber], linestyle = linestyles[linenumber])
            
      plt.plot(np.arange(frame_max_plus1)*time_per_frame, intM0, \
               color = linecolors[linenumber], linestyle = linestyles[linenumber])

      linenumber += 1
      plt.title("N_z = "+Nz+", nonuniform frac = "+nonuniform_frac_str)
      plt.xlabel('Time, seconds')
      plt.ylabel('$\int M_0 dx$')

  plt.tight_layout()
  plt.savefig(outDir+'intM0_trace'+figureFileFormat)

def plot_intM0_combined():
  print("Plotting integral of M0 trace combined")
  plt.figure( figsize = figsize_dim)
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]
      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'

      skip = skip_sims(Nz, nonuniform_frac_str)
      noskip = plot_only_sims(Nz, nonuniform_frac_str)
      if skip or not noskip:
        continue

      intM0 = np.zeros(frame_max_plus1)
      nodes_Z = loadZ(filebase)
      for i in range(frame_max_plus1):
        x_M0, dataOut_M0 = loadM0(i, filebase)
        intM0[i] = np.trapz(dataOut_M0, nodes_Z)

      # t, M0 = loadM0_integrated(filebase)
      
      # plt.plot(t, M0, \
      #          color = linecolors[linenumber], linestyle = linestyles[linenumber])
      plt.plot(np.arange(frame_max_plus1)*time_per_frame, intM0, \
               color = linecolors[linenumber], linestyle = linestyles[linenumber],\
                label = 'Nz = '+Nz+', nonuniform frac = '+nonuniform_frac_str)
      linenumber += 1

  plt.xlabel('Time, seconds')
  plt.ylabel('$\int M_0 dx$')
  plt.title('Integral of M0')
  plt.legend()
  plt.tight_layout()
  plt.savefig(outDir+'intM0_trace_combined'+figureFileFormat)

def plot_final_density():
  print("Plotting final density")
  plt.figure( figsize = figsize_dim)
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]
      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'

      skip = skip_sims(Nz, nonuniform_frac_str)
      if skip:
        linenumber += 1
        continue
      
      nodes_Z = loadZ(filebase)
      xn, dens, U, Tpar, Tperp = load_mom(frame_max_plus1-1, filebase)  
      # xn, dens, U, Tpar, Tperp = load_mom(0, filebase)


      plt.subplot(len(Nz_arr), len(nonuniform_frac), linenumber+1)
      plt.plot(nodes_Z, dens, \
          color = linecolors[linenumber], linestyle = linestyles[linenumber])
      linenumber += 1
      plt.yscale('log')
      plt.title("N_z = "+Nz+", nonuniform frac = "+nonuniform_frac_str)
      plt.xlabel('Z, m')
      plt.ylabel('Density, $m^{-3}$')
  plt.tight_layout()
  plt.savefig(outDir+'density_final'+figureFileFormat)

def plot_final_upar():
  print("Plotting final upar")
  plt.figure( figsize = figsize_dim)
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]
      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'

      skip = skip_sims(Nz, nonuniform_frac_str)
      if skip:
        linenumber += 1
        continue
      
      nodes_Z = loadZ(filebase)
      xn, dens, U, Tpar, Tperp = load_mom(frame_max_plus1-1, filebase)

      plt.subplot(len(Nz_arr), len(nonuniform_frac), linenumber+1)
      plt.plot(nodes_Z, U, \
          color = linecolors[linenumber], linestyle = linestyles[linenumber])
      linenumber += 1
      plt.title("N_z = "+Nz+", nonuniform frac = "+nonuniform_frac_str)
      plt.xlabel('Z, m')
      plt.ylabel('$U_{||}, m/s$')
  plt.tight_layout()
  plt.savefig(outDir+'upar_final'+figureFileFormat)

def plot_final_Tpar():
  print("Plotting final Tpar")
  plt.figure( figsize = figsize_dim)
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]
      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'

      skip = skip_sims(Nz, nonuniform_frac_str)
      if skip:
        linenumber += 1
        continue
      
      nodes_Z = loadZ(filebase)
      xn, dens, U, Tpar, Tperp = load_mom(frame_max_plus1-1, filebase)

      plt.subplot(len(Nz_arr), len(nonuniform_frac), linenumber+1)
      plt.plot(nodes_Z, Tpar, \
          color = linecolors[linenumber], linestyle = linestyles[linenumber])
      linenumber += 1
      plt.title("N_z = "+Nz+", nonuniform frac = "+nonuniform_frac_str)
      plt.xlabel('Z, m')
      plt.ylabel('$T_{||}, keV$')
  plt.tight_layout()
  plt.savefig(outDir+'Tpar_final'+figureFileFormat)

def plot_final_Tperp():
  print("Plotting final Tperp")
  plt.figure( figsize = figsize_dim)
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      Nz = Nz_arr[iZ]
      nonuniform_frac_str = nonuniform_frac[iS]
      filebase = './'+Nz+'/'+nonuniform_frac_str+'/'

      skip = skip_sims(Nz, nonuniform_frac_str)
      if skip:
        linenumber += 1
        continue
      
      nodes_Z = loadZ(filebase)
      xn, dens, U, Tpar, Tperp = load_mom(frame_max_plus1-1, filebase)

      plt.subplot(len(Nz_arr), len(nonuniform_frac), linenumber+1)
      plt.plot(nodes_Z, Tperp, \
          color = linecolors[linenumber], linestyle = linestyles[linenumber])
      linenumber += 1
      plt.title("N_z = "+Nz+", nonuniform frac = "+nonuniform_frac_str)
      plt.xlabel('Z, m')
      plt.ylabel('$T_{\perp}, keV$')
  plt.tight_layout()
  plt.savefig(outDir+'Tperp_final'+figureFileFormat)

def plot_final_BiMax_moms_combined():
  print("Plotting final BiMaxwellian moments")
  plt.figure( figsize = figsize_dim)
  linenumber = 0
  for iZ in range(len(Nz_arr)):
    for iS in range(len(nonuniform_frac)):
      for iF in range(len(frame_arr)):
        Nz = Nz_arr[iZ]
        nonuniform_frac_str = nonuniform_frac[iS]
        frame = frame_arr[iF]
        time = frame*time_per_frame
        filebase = './'+Nz+'/'+nonuniform_frac_str+'/'

        skip = skip_sims(Nz, nonuniform_frac_str)
        noskip = plot_only_sims(Nz, nonuniform_frac_str)
        if skip or not noskip:
          continue

        nodes_Z = loadZ(filebase)
        xn, dens, U, Tpar, Tperp = load_mom(frame, filebase)

        plt.subplot(2,2,1)
        plt.plot(nodes_Z, dens, \
            color = linecolors[linenumber], linestyle = linestyles[linenumber], \
            label = 'Nz = '+Nz+', nonuniform frac = '+nonuniform_frac_str\
              +', t='+str(time)+' s')
        plt.subplot(2,2,2)
        plt.plot(nodes_Z, U, \
            color = linecolors[linenumber], linestyle = linestyles[linenumber], \
            label = 'Nz = '+Nz+', nonuniform frac = '+nonuniform_frac_str\
              +', t='+str(time)+' s')
        plt.subplot(2,2,3)
        plt.plot(nodes_Z, Tpar, \
            color = linecolors[linenumber], linestyle = linestyles[linenumber], \
            label = 'Nz = '+Nz+', nonuniform frac = '+nonuniform_frac_str\
              +', t='+str(time)+' s')
        plt.subplot(2,2,4)
        plt.plot(nodes_Z, Tperp, \
            color = linecolors[linenumber], linestyle = linestyles[linenumber], \
            label = 'Nz = '+Nz+', nonuniform frac = '+nonuniform_frac_str\
              +', t='+str(time)+' s')
        linenumber += 1

  plt.subplot(2,2,1)
  plt.xlabel('Z, m')
  plt.ylabel('Density, $m^{-3}$')
  plt.title('Final density')
  plt.subplot(2,2,2)
  plt.xlabel('Z, m')
  plt.ylabel('$U_{||}, m/s$')
  plt.title('Final upar')
  plt.legend()
  plt.subplot(2,2,3)
  plt.xlabel('Z, m')
  plt.ylabel('$T_{||}, keV$')
  plt.title('Final Tpar')
  plt.subplot(2,2,4)
  plt.xlabel('Z, m')
  plt.ylabel('$T_{\perp}, keV$')
  plt.title('Final Tperp')
  plt.tight_layout()
  plt.savefig(outDir+'BiMaxMoms_combined'+figureFileFormat)


## Main function
def main():
  processes = []
  plot_functions = [
    # plot_potential_trace_seperate,
    # plot_dndt_trace,
    # plot_intM0,
    # plot_final_density,
    # plot_final_upar,
    # plot_final_Tpar,
    # plot_final_Tperp,
    plot_potential_trace_combined,
    plot_dndt_trace_combined,
    plot_intM0_combined,
    plot_final_BiMax_moms_combined

  ]
  for plot_func in plot_functions:
    p = multiprocessing.Process(target=plot_func)
    processes.append(p)
    p.start()
  for p in processes:
    p.join()
  print("Done")

if __name__ == "__main__":
  main()
