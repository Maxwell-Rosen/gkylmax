import numpy as np
from scipy.integrate import quad
from scipy.special import erf
#[ Append postgkyl wrappers.
from scipy import special
from scipy import optimize
import postgkyl as pg


# Calculate the Pastukhov potential confinement time from an adiabatic electron simulation
# directory = '/scratch/gpfs/mr1884/scratch/gkylmax/stellar-lorentzian1x-exploration/stellar-lorentzian1x-orbit-average-nu-ie-maxwellian-odist-time-dilation-pos-extra-tdil20'
directory = '/scratch/gpfs/mr1884/scratch/gkylmax/stellar-lorentzian1x-kinetic-exploration/initial-test'
sim_name = 'gk_lorentzian_mirror'
species = 'elc'
frame_num = 140
source_frame = 0
z_mirror_throat = 0.98
simulation_edge = 2.5


# Extract details about the magnetic field
bmag_data = pg.GData(	directory + f'/{sim_name}-bmag.gkyl',	mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl')
pg.GInterpModal(bmag_data).interpolate(0, overwrite=True)
_, bmag_sel = pg.data.select(bmag_data, z0='0.0')
Bmag0 = float(np.squeeze(bmag_sel))
_, bmag_sel = pg.data.select(bmag_data, z0=f'{z_mirror_throat}')
BmagThroat = float(np.squeeze(bmag_sel))
_, bmag_sel = pg.data.select(bmag_data, z0=f'{simulation_edge}')
BmagExpander = float(np.squeeze(bmag_sel))
R = BmagThroat / Bmag0  # Read from postgkyl

# Measure values of the bimaxwellian moments (midplane density, upar at the throat, upar at the wall)

# Find the integrated density inside the trap
m0_data = pg.GData(	directory + f'/{sim_name}-{species}_M0_{frame_num}.gkyl',	mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl')
pg.GInterpModal(m0_data).interpolate(comp=0, overwrite=True)
pg.data.select(m0_data, z0=f'{-z_mirror_throat}:{z_mirror_throat}', overwrite=True)
_, intM0dx = pg.tools.integrate(m0_data, 0, overwrite=True)
intM0dx = float(np.squeeze(intM0dx))

m0_data_ion = pg.GData(	directory + f'/{sim_name}-ion_M0_{frame_num}.gkyl',	mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl')
pg.GInterpModal(m0_data_ion).interpolate(comp=0, overwrite=True)
pg.data.select(m0_data_ion, z0=f'{-z_mirror_throat}:{z_mirror_throat}', overwrite=True)
_, intM0dx_ion_raw = pg.tools.integrate(m0_data_ion, 0, overwrite=True)
intM0dx_ion = float(np.squeeze(intM0dx_ion_raw))

# Find the integrated source inside the trap
m0src_data = pg.GData(	directory + f'/{sim_name}-{species}_source_M0_{source_frame}.gkyl',	mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl')
pg.GInterpModal(m0src_data).interpolate(comp=0, overwrite=True)
pg.data.select(m0src_data, z0=f'{-z_mirror_throat}:{z_mirror_throat}', overwrite=True)
_, intM0Sdx = pg.tools.integrate(m0src_data, 0, overwrite=True)
intM0Sdx = float(np.squeeze(intM0Sdx))

m0src_data_ion = pg.GData(	directory + f'/{sim_name}-ion_source_M0_{source_frame}.gkyl',	mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl')
pg.GInterpModal(m0src_data_ion).interpolate(comp=0, overwrite=True)
pg.data.select(m0src_data_ion, z0=f'{-z_mirror_throat}:{z_mirror_throat}', overwrite=True)
_, intM0Sdx_ion_raw = pg.tools.integrate(m0src_data_ion, 0, overwrite=True)
intM0Sdx_ion = float(np.squeeze(intM0Sdx_ion_raw))

nusum_data = pg.GData(directory + f'/{sim_name}-{species}_lbo_nu_sum_{frame_num}.gkyl', mapc2p_name=directory + f'/{sim_name}-mc2nu_pos_deflated.gkyl')
pg.GInterpModal(nusum_data).interpolate(comp=0, overwrite=True)
_, nu_sum_sel = pg.data.select(nusum_data, z0=0.0)
nuElc = float(np.squeeze(nu_sum_sel))

eps0, mu0 = 8.8541878176204e-12, 1.2566370614359e-06
eV        = 1.602176487e-19
qe, qi    = -1.602176487e-19, 1.602176487e-19
me, mp    = 9.10938215e-31, 1.672621637e-27

Zpfl = 1.0 # Z in Najmabadi (1984). ee collisions would be 1.0, ee+ei collisions would be 1.5 (lambda_ei = lambda_ee, n_i = n_e)

def Rosen_Dougherty_confinement_time(P, R, ZpFl, coeff = 0):
    w_term = np.sqrt(1 + 2*P/(R*ZpFl)) 
    a_term  = np.sqrt(P + np.log(w_term)) 

    Loss_Rosen = 1/nuElc * \
        1/(2*ZpFl / (np.log((w_term+1)/(w_term-1)) - coeff) * (1-erf(a_term)))
    
    return Loss_Rosen

def tau_pe(R):
  tau_pe = intM0dx / intM0Sdx
  return tau_pe

def tau_pi(R):
  # These don't agree. I think we need to calculate distribution function density, but we are using guiding center, so we need to calculate the polarization density
  tau_pi = intM0dx_ion / intM0Sdx_ion
  return tau_pi

def rootEq_dough(x,R):
  return tau_pe(R)-Rosen_Dougherty_confinement_time(x,R,Zpfl,1.117)

ephi_over_Te = optimize.ridder(rootEq_dough, Zpfl, 20, args=(R))

print(f"Dougherty e*phi/Te = {ephi_over_Te:.3f}")

print(f"tau_pe = {tau_pe(R):.3f} s")
print(f"tau_pi = {tau_pi(R):.3f} s")