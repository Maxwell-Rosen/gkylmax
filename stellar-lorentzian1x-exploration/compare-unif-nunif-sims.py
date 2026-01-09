import numpy as np
import matplotlib.pyplot as plt
import postgkyl as pg
from matplotlib.colors import LogNorm
import multiprocessing
from scipy.integrate import cumulative_trapezoid as cumtrapz
import imageio.v2 as imageio
from scipy.optimize import curve_fit

frame = 20

# sim_uniform = 'stellar-lorentzian1x-orbit-average-time-dilation-unif-Nz800-fmin'
# sim_nonuniform = 'stellar-lorentzian1x-orbit-average-time-dilation-nunif-Nz800-fmin-compare'
# sim_nonuniform = 'stellar-lorentzian1x-orbit-average-time-dilation-nunif-Nz200-fmin-compare'
sim_uniform = 'stellar-lorentzian1x-orbit-average-beams'
# sim_nonuniform = 'stellar-lorentzian1x-orbit-average-beams-nunif-hires'
sim_nonuniform = 'stellar-lorentzian1x-orbit-average-beams-nunif-hires-less-intense'

uniformBiMax_name = sim_uniform+'/gk_lorentzian_mirror-ion_BiMaxwellianMoments_'+str(frame)+'.gkyl'
nonunifBiMax_name = sim_nonuniform+'/gk_lorentzian_mirror-ion_BiMaxwellianMoments_'+str(frame)+'.gkyl'

uniformField_name = sim_uniform+'/gk_lorentzian_mirror-field_'+str(frame)+'.gkyl'
nonunifField_name = sim_nonuniform+'/gk_lorentzian_mirror-field_'+str(frame)+'.gkyl'
uniformPosMap_name = sim_uniform+'/gk_lorentzian_mirror-mc2nu_pos.gkyl'
nonunifPosMap_name = sim_nonuniform+'/gk_lorentzian_mirror-mc2nu_pos.gkyl'

uniformBiMax_pgdata = pg.GData(uniformBiMax_name)
nonunifBiMax_pgdata = pg.GData(nonunifBiMax_name)
uniformField_pgdata = pg.GData(uniformField_name)
nonunifField_pgdata = pg.GData(nonunifField_name)
uniformPosMap_pgdata = pg.GData(uniformPosMap_name)
nonunifPosMap_pgdata = pg.GData(nonunifPosMap_name)

uniformBiMax_pginterp = pg.GInterpModal(uniformBiMax_pgdata, 1, 'ms')
nonunifBiMax_pginterp = pg.GInterpModal(nonunifBiMax_pgdata, 1, 'ms')
uniformField_pginterp = pg.GInterpModal(uniformField_pgdata, 1, 'ms')
nonunifField_pginterp = pg.GInterpModal(nonunifField_pgdata, 1, 'ms')
uniformPosMap_pginterp = pg.GInterpModal(uniformPosMap_pgdata, 1, 'ms')
nonunifPosMap_pginterp = pg.GInterpModal(nonunifPosMap_pgdata, 1, 'ms')

coords, uniformDens = uniformBiMax_pginterp.interpolate(0)
coords, uniformUpar = uniformBiMax_pginterp.interpolate(1)
coords, uniformTpar = uniformBiMax_pginterp.interpolate(2)
coords, uniformTperp = uniformBiMax_pginterp.interpolate(3)

coords, nonunifDens = nonunifBiMax_pginterp.interpolate(0)
coords, nonunifUpar = nonunifBiMax_pginterp.interpolate(1)
coords, nonunifTpar = nonunifBiMax_pginterp.interpolate(2)
coords, nonunifTperp = nonunifBiMax_pginterp.interpolate(3)

coords, uniformField = uniformField_pginterp.interpolate()
coords, nonunifField = nonunifField_pginterp.interpolate()

coords, uniformPosMap = uniformPosMap_pginterp.interpolate(2)
coords, nonunifPosMap = nonunifPosMap_pginterp.interpolate(2)

# Convert T into eV
# Universal constants (matching input_file.c)
eps0 = 8.8541878128e-12  # F/m (permittivity of free space)
mu0 = 4 * np.pi * 1e-7  # H/m (permeability of free space)
eV = 1.602176634e-19  # J (elementary charge / electron volt)
mp = 1.67262192369e-27  # kg (proton mass)
me = 9.1093837015e-31  # kg (electron mass)

# Plasma parameters (from input_file.c)
mi = 2.014 * mp  # kg (deuteron mass)

uniformTpar = uniformTpar * mi / eV
uniformTperp = uniformTperp * mi / eV
nonunifTpar = nonunifTpar * mi / eV
nonunifTperp = nonunifTperp * mi / eV

# Convert potential into e phi / Te
Te = 940 * eV  # J (electron temperature)
uniformField = uniformField * eV / Te
nonunifField = nonunifField * eV / Te

fig, axes = plt.subplots(3, 2, figsize=(12, 10))

# Density
ax = axes[0, 0]
ax.plot(uniformPosMap[:,0], uniformDens[:,0], label='Uniform', linewidth=2)
ax.plot(nonunifPosMap[:,0], nonunifDens[:,0], label='Nonuniform', linestyle='--', linewidth=2)
ax.set_ylabel(r"Density ($m^{-3}$)", fontsize=11)
# ax.set_yscale('log')
ax.legend(loc='best', fontsize=9)
ax.grid(True, alpha=0.3)

# Upar
ax = axes[0, 1]
ax.plot(uniformPosMap[:,0], uniformUpar[:,0], linewidth=2)
ax.plot(nonunifPosMap[:,0], nonunifUpar[:,0], linestyle='--', linewidth=2)
ax.set_ylabel(r"$u_\parallel$ (m/s)", fontsize=11)
ax.grid(True, alpha=0.3)

# Tpar
ax = axes[1, 0]
ax.plot(uniformPosMap[:,0], uniformTpar[:,0], linewidth=2)
ax.plot(nonunifPosMap[:,0], nonunifTpar[:,0], linestyle='--', linewidth=2)
ax.set_ylabel(r"$T_\parallel$ (eV)", fontsize=11)
ax.grid(True, alpha=0.3)

# Tperp
ax = axes[1, 1]
ax.plot(uniformPosMap[:,0], uniformTperp[:,0], linewidth=2)
ax.plot(nonunifPosMap[:,0], nonunifTperp[:,0], linestyle='--', linewidth=2)
ax.set_ylabel(r"$T_\perp$ (eV)", fontsize=11)
ax.grid(True, alpha=0.3)

# Electric potential
ax = axes[2, 0]
ax.plot(uniformPosMap[:,0], uniformField[:,0], linewidth=2)
ax.plot(nonunifPosMap[:,0], nonunifField[:,0], linestyle='--', linewidth=2)
ax.set_xlabel("Normalized field line length", fontsize=11)
ax.set_ylabel(r"$e\phi / T_e$", fontsize=11)
ax.grid(True, alpha=0.3)

# Remove the unused subplot
fig.delaxes(axes[2, 1])

plt.suptitle(f"Grid Comparison at Frame {frame}", fontsize=13, y=0.995)
plt.tight_layout()
# plt.savefig('./plots/moments_comparison.pdf')
plt.show()


# Follow this up by plotting /gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl from each of the simulations on top of each other
uniformPosDeflated_pgdata = pg.GData(sim_uniform+'/gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl')
nonunifPosDeflated_pgdata = pg.GData(sim_nonuniform+'/gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl')

uniformPosDeflated_pginterp = pg.GInterpModal(uniformPosDeflated_pgdata, 1, 'ms')
nonunifPosDeflated_pginterp = pg.GInterpModal(nonunifPosDeflated_pgdata, 1, 'ms')

coords, uniformPosDeflated = uniformPosDeflated_pginterp.interpolate()
coords, nonunifPosDeflated = nonunifPosDeflated_pginterp.interpolate()

# Plot gk_lorentzian_mirror-bmag.gkyl on top of each other using the posmap for the x-axis
uniformBmag_pgdata = pg.GData(sim_uniform+'/gk_lorentzian_mirror-bmag.gkyl')
nonunifBmag_pgdata = pg.GData(sim_nonuniform+'/gk_lorentzian_mirror-bmag.gkyl')

uniformBmag_pginterp = pg.GInterpModal(uniformBmag_pgdata, 1, 'ms')
nonunifBmag_pginterp = pg.GInterpModal(nonunifBmag_pgdata, 1, 'ms')

coords, uniformBmag = uniformBmag_pginterp.interpolate()
coords, nonunifBmag = nonunifBmag_pginterp.interpolate()

# Plot gk_lorentzian_mirror-B3.gkyl on top of each other using the posmap for the x-axis
uniformB3_pgdata = pg.GData(sim_uniform+'/gk_lorentzian_mirror-B3.gkyl')
nonunifB3_pgdata = pg.GData(sim_nonuniform+'/gk_lorentzian_mirror-B3.gkyl')

uniformB3_pginterp = pg.GInterpModal(uniformB3_pgdata, 1, 'ms')
nonunifB3_pginterp = pg.GInterpModal(nonunifB3_pgdata, 1, 'ms')

coords, uniformB3 = uniformB3_pginterp.interpolate()
coords, nonunifB3 = nonunifB3_pginterp.interpolate()

# Plot gk_lorentzian_mirror-jacobgeo.gkyl on top of each other using the posmap for the x-axis
uniformJacob_pgdata = pg.GData(sim_uniform+'/gk_lorentzian_mirror-jacobgeo.gkyl')
nonunifJacob_pgdata = pg.GData(sim_nonuniform+'/gk_lorentzian_mirror-jacobgeo.gkyl')

uniformJacob_pginterp = pg.GInterpModal(uniformJacob_pgdata, 1, 'ms')
nonunifJacob_pginterp = pg.GInterpModal(nonunifJacob_pgdata, 1, 'ms')

coords, uniformJacob = uniformJacob_pginterp.interpolate()
coords, nonunifJacob = nonunifJacob_pginterp.interpolate()

# Create figure with the four additional plots
fig2, axes2 = plt.subplots(4, 2, figsize=(14, 12))

# Calculate global x-limits from all position maps
x_min = min(uniformPosMap[:,0].min(), nonunifPosMap[:,0].min())
x_max = max(uniformPosMap[:,0].max(), nonunifPosMap[:,0].max())

# Function to split data at x=0
def split_at_zero(x, y):
    """Split data into left (x<0) and right (x>=0) portions"""
    # Ensure x is 1D for comparison
    if x.ndim > 1:
        x_1d = x.flatten()
    else:
        x_1d = x
    
    # Ensure y is 1D for indexing
    if y.ndim > 1:
        y_1d = y.flatten()
    else:
        y_1d = y
    
    left_mask = x_1d < 0
    right_mask = x_1d >= 0
    return x_1d[left_mask], y_1d[left_mask], x_1d[right_mask], y_1d[right_mask]

# Plot deflated position map
# Handle different possible shapes of the deflated position map
if uniformPosDeflated.ndim > 1 and uniformPosDeflated.shape[1] > 1:
    u_x_left, u_y_left, u_x_right, u_y_right = split_at_zero(np.linspace(uniformPosMap[0], uniformPosMap[-1], uniformPosMap.shape[0]), uniformPosDeflated[:,2])
    n_x_left, n_y_left, n_x_right, n_y_right = split_at_zero(np.linspace(uniformPosMap[0], uniformPosMap[-1], nonunifPosDeflated.shape[0]), nonunifPosDeflated[:,2])
else:
    u_x_left, u_y_left, u_x_right, u_y_right = split_at_zero(np.linspace(uniformPosMap[0], uniformPosMap[-1], uniformPosDeflated.shape[0]), uniformPosDeflated[:,0])
    n_x_left, n_y_left, n_x_right, n_y_right = split_at_zero(np.linspace(uniformPosMap[0], uniformPosMap[-1], nonunifPosDeflated.shape[0]), nonunifPosDeflated[:,0])

ax_left = axes2[0, 0]
ax_left.plot(u_x_left, u_y_left, label='Uniform', linewidth=2)
ax_left.plot(n_x_left, n_y_left, label='Nonuniform', linestyle='--', linewidth=2)
ax_left.set_ylabel(r"Deflated Position Map", fontsize=11)
ax_left.legend(loc='best', fontsize=9)
ax_left.grid(True, alpha=0.3)
ax_left.set_title('Left side (linear scale)', fontsize=10)
ax_left.set_xlim(x_min, 0)

ax_right = axes2[0, 1]
ax_right.plot(u_x_right, u_y_right, linewidth=2)
ax_right.plot(n_x_right, n_y_right, linestyle='--', linewidth=2)
ax_right.set_ylabel(r"Deflated Position Map", fontsize=11)
ax_right.yaxis.tick_right()
ax_right.yaxis.set_label_position("right")
ax_right.set_yscale('log')
ax_right.grid(True, alpha=0.3)
ax_right.set_title('Right side (log scale)', fontsize=10)
ax_right.set_xlim(0, x_max)

# Plot bmag
u_x_left, u_y_left, u_x_right, u_y_right = split_at_zero(uniformPosMap[:,0], uniformBmag[:,0])
n_x_left, n_y_left, n_x_right, n_y_right = split_at_zero(nonunifPosMap[:,0], nonunifBmag[:,0])

ax_left = axes2[1, 0]
ax_left.plot(u_x_left, u_y_left, linewidth=2)
ax_left.plot(n_x_left, n_y_left, linestyle='--', linewidth=2)
ax_left.set_ylabel(r"$|B|$ (T)", fontsize=11)
ax_left.grid(True, alpha=0.3)
ax_left.set_xlim(x_min, 0)

ax_right = axes2[1, 1]
ax_right.plot(u_x_right, u_y_right, linewidth=2)
ax_right.plot(n_x_right, n_y_right, linestyle='--', linewidth=2)
ax_right.set_ylabel(r"$|B|$ (T)", fontsize=11)
ax_right.yaxis.tick_right()
ax_right.yaxis.set_label_position("right")
ax_right.set_yscale('log')
ax_right.grid(True, alpha=0.3)
ax_right.set_xlim(0, x_max)

# Plot B3
u_x_left, u_y_left, u_x_right, u_y_right = split_at_zero(uniformPosMap[:,0], uniformB3[:,0])
n_x_left, n_y_left, n_x_right, n_y_right = split_at_zero(nonunifPosMap[:,0], nonunifB3[:,0])

ax_left = axes2[2, 0]
ax_left.plot(u_x_left, u_y_left, linewidth=2)
ax_left.plot(n_x_left, n_y_left, linestyle='--', linewidth=2)
ax_left.set_ylabel(r"$B^3$ (T$^3$)", fontsize=11)
ax_left.grid(True, alpha=0.3)
ax_left.set_xlim(x_min, 0)

ax_right = axes2[2, 1]
ax_right.plot(u_x_right, u_y_right, linewidth=2)
ax_right.plot(n_x_right, n_y_right, linestyle='--', linewidth=2)
ax_right.set_ylabel(r"$B^3$ (T$^3$)", fontsize=11)
ax_right.yaxis.tick_right()
ax_right.yaxis.set_label_position("right")
ax_right.set_yscale('log')
ax_right.grid(True, alpha=0.3)
ax_right.set_xlim(0, x_max)

# Plot Jacobian inverse
u_x_left, u_y_left, u_x_right, u_y_right = split_at_zero(uniformPosMap[:,0], uniformJacob[:,0])
n_x_left, n_y_left, n_x_right, n_y_right = split_at_zero(nonunifPosMap[:,0], nonunifJacob[:,0])

ax_left = axes2[3, 0]
ax_left.plot(u_x_left, u_y_left, linewidth=2)
ax_left.plot(n_x_left, n_y_left, linestyle='--', linewidth=2)
ax_left.set_xlabel("Normalized field line length", fontsize=11)
ax_left.set_ylabel(r"$J_{\mathrm{geo}}$", fontsize=11)
ax_left.grid(True, alpha=0.3)
ax_left.set_xlim(x_min, 0)

ax_right = axes2[3, 1]
ax_right.plot(u_x_right, u_y_right, linewidth=2)
ax_right.plot(n_x_right, n_y_right, linestyle='--', linewidth=2)
ax_right.set_xlabel("Normalized field line length", fontsize=11)
ax_right.set_ylabel(r"$J_{\mathrm{geo}}$", fontsize=11)
ax_right.yaxis.tick_right()
ax_right.yaxis.set_label_position("right")
ax_right.set_yscale('log')
ax_right.grid(True, alpha=0.3)
ax_right.set_xlim(0, x_max)

plt.suptitle("Magnetic Field and Grid Comparison", fontsize=13, y=0.995)
plt.tight_layout()
plt.subplots_adjust(wspace=0.0)  # Reduce horizontal space between subplots
# plt.savefig('./plots/geometry_comparison.pdf')
plt.show()