
import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from matplotlib.colors import LogNorm

frame_number = 4600
polyOrder = 1

# Constants
eV = 1.602176634e-19  # electron volt in Joules
me = 9.10938356e-31  # electron mass in kg
mi = 2.014 * 1.6726219e-27  # ion mass in kg


def loadphi(frame_number):
  filename_phi = str('Field/gk_wham-field_'+str(frame_number)+'.gkyl')
  pgData_phi = pg.GData(filename_phi)
  pgInterp_phi = pg.GInterpModal(pgData_phi, polyOrder, 'ms')
  x_phi, dataOut_phi = pgInterp_phi.interpolate()
  return dataOut_phi

def get_temp(frame_number):
  filename_elc = str('BiMaxwellianMoments/gk_wham-elc_BiMaxwellianMoments_'+str(frame_number)+'.gkyl')
  pgData_elc = pg.GData(filename_elc)
  pgInterp_elc = pg.GInterpModal(pgData_elc, polyOrder, 'ms')
  coords, Tpar_elc = pgInterp_elc.interpolate(2)
  coords, Tperp_elc = pgInterp_elc.interpolate(3)
  midpoint = int(len(coords[0])/2)
  Temp = (Tpar_elc[0,midpoint,0] + 2*Tperp_elc[0,midpoint,0])/3 * me / eV
  return Temp

def plot_ephiTe(frame_number):
  phi = loadphi(frame_number)
  Temp = get_temp(frame_number)

  ephiTe = phi / Temp  # Convert to electric potential in eV

  # Create a figure and axis
  fig, ax = plt.subplots(figsize=(10, 6))

  # Plot the potential
  im = ax.imshow(ephiTe, extent=[-2, 2, 5e-5, 3e-3], origin='lower', aspect='auto', cmap='inferno')
  
  # Add a colorbar
  cbar = fig.colorbar(im, ax=ax)
  cbar.set_label('Electric Potential (Thermal V)', rotation=270, labelpad=15)

  # Set labels and title
  ax.set_xlabel('Z (m)')
  ax.set_ylabel('psi')
  ax.set_title(f'Electric Potential at Frame {frame_number}, T_e = {Temp:.2f} eV')

  # Save the plot
  plt.savefig(f'python-plots/ephiTe_frame_{frame_number}.png')
  plt.close(fig)


plot_ephiTe(frame_number)