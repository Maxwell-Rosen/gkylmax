# Compare a handwritten calculation of the cross species moments for collisions with the reulsts in gk_lorentzian_mirror-ion_lbo_cross_elc_prim_moms_0.gkyl.
# Cross vt2 is computed as
#vt2_sr = vts^2 + delta / 2 * (1+beta)/(1 + m_s / m_r) * (v_tr^2 - m_s / m_r * v_ts^2 + (u_s - u_r)^2 / dv)
# where dv = 2, delta = 2, beta = 1
# Compute vt2_ie. We can find the relevant information in the maxwellian moments file

import numpy as np
import postgkyl as pg
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ---- Physical constants ----
GKYL_PROTON_MASS   = 1.6726219e-27    # kg
GKYL_ELECTRON_MASS = 9.10938e-31     # kg
eV = 1.602176634e-19                  # J

# ---- Simulation parameters (from sim.c) ----
mi = 2.014 * GKYL_PROTON_MASS
me = GKYL_ELECTRON_MASS

# LBO operator parameters (1x2v GK)
dv    = 2   # number of velocity dimensions (vpar + mu)
delta = 2
beta  = 1

# ---- File paths ----
dataDir = '/scratch/gpfs/mr1884/scratch/gkylmax/stellar-lorentzian1x-kinetic-exploration/debug-nuie/'
simName = 'gk_lorentzian_mirror'

poly_order = 1   # DG polynomial order used in the simulation

def interp_comp(fname, comp_slice):
    """Load a component slice from a 1-D Gkeyll file and interpolate to grid points.

    comp_slice : str like '0:2' or '2:4' selects which DG components to load.
    Returns (z_grid, values) where z_grid has Ncells*poly_order points and
    values has the same length.
    """
    d = pg.GData(dataDir + simName + '-' + fname, comp=comp_slice)
    interp = pg.GInterpModal(d, poly_order, 'ms')
    z_grid, vals = interp.interpolate()
    # z_grid[0] has Npts+1 nodes; evaluation points are at sub-cell midpoints
    z_eval = 0.5 * (z_grid[0][:-1] + z_grid[0][1:])
    return z_eval, vals[:, 0]

# ---- Load and interpolate self primary moments ----
# prim_moms layout: comp 0:2 = u_par DG coeffs, comp 2:4 = vtsq DG coeffs
z,  u_ion  = interp_comp('ion_lbo_self_prim_moms_0.gkyl', '0:2')
_,  vts2_i = interp_comp('ion_lbo_self_prim_moms_0.gkyl', '2:4')

_,  u_elc  = interp_comp('elc_lbo_self_prim_moms_0.gkyl', '0:2')
_,  vts2_e = interp_comp('elc_lbo_self_prim_moms_0.gkyl', '2:4')

# ---- Load and interpolate Gkeyll cross primary moments ----
_, u_cross_ion_gk    = interp_comp('ion_lbo_cross_elc_prim_moms_0.gkyl', '0:2')
_, vts2_cross_ion_gk = interp_comp('ion_lbo_cross_elc_prim_moms_0.gkyl', '2:4')

_, u_cross_elc_gk    = interp_comp('elc_lbo_cross_ion_prim_moms_0.gkyl', '0:2')
_, vts2_cross_elc_gk = interp_comp('elc_lbo_cross_ion_prim_moms_0.gkyl', '2:4')

# ---- Handwritten cross-vtsq formula ----
# vt2_sr = vts^2 + delta/2*(1+beta)/(1+ms/mr) * (vtr^2 - ms/mr*vts^2 + (u_s-u_r)^2/dv)
prefac = (delta / 2.0) * (1.0 + beta)   # = 2 for delta=2, beta=1

# Case 1: s = ion,  r = elc  (compare with ion_lbo_cross_elc_prim_moms)
ms_over_mr_ie = mi / me
vts2_cross_formula_ion = (
    vts2_i
    + prefac / (1.0 + ms_over_mr_ie)
    * (vts2_e - ms_over_mr_ie * vts2_i + (u_ion - u_elc)**2 / dv)
)

# Case 2: s = elc,  r = ion  (compare with elc_lbo_cross_ion_prim_moms)
ms_over_mr_ei = me / mi
vts2_cross_formula_elc = (
    vts2_e
    + prefac / (1.0 + ms_over_mr_ei)
    * (vts2_i - ms_over_mr_ei * vts2_e + (u_elc - u_ion)**2 / dv)
)

# ---- Plot ----
fig = plt.figure(figsize=(13, 10))
gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.48, wspace=0.38)

# --- Row 0: parallel drift u_cross ---
ax00 = fig.add_subplot(gs[0, 0])
ax00.set_title(r'ION cross drift  ($s$=ion, $r$=elc)')
ax00.plot(z, u_ion,          color='C0', ls='--', lw=1.2, label=r'$u_{\rm ion}$ self')
ax00.plot(z, u_elc,          color='C1', ls='--', lw=1.2, label=r'$u_{\rm elc}$ self')
ax00.plot(z, u_cross_ion_gk, color='C0', ls='-',  lw=2,   label=r'$u_{\rm cross,ion}$ Gkeyll')
ax00.set_xlabel(r'$z$ (m)')
ax00.set_ylabel(r'$u_\parallel$ (m s$^{-1}$)')
ax00.legend(fontsize=8, frameon=False)
ax00.axhline(0, color='k', lw=0.5, ls=':')

ax01 = fig.add_subplot(gs[0, 1])
ax01.set_title(r'ELC cross drift  ($s$=elc, $r$=ion)')
ax01.plot(z, u_ion,           color='C0', ls='--', lw=1.2, label=r'$u_{\rm ion}$ self')
ax01.plot(z, u_elc,           color='C1', ls='--', lw=1.2, label=r'$u_{\rm elc}$ self')
ax01.plot(z, u_cross_elc_gk,  color='C1', ls='-',  lw=2,   label=r'$u_{\rm cross,elc}$ Gkeyll')
ax01.set_xlabel(r'$z$ (m)')
ax01.set_ylabel(r'$u_\parallel$ (m s$^{-1}$)')
ax01.legend(fontsize=8, frameon=False)
ax01.axhline(0, color='k', lw=0.5, ls=':')

# --- Row 1: thermal velocity squared vtsq_cross ---
ax10 = fig.add_subplot(gs[1, 0])
ax10.set_title(r'ION cross vtsq:  formula ($s$=ion, $r$=elc) vs Gkeyll')
ax10.plot(z, vts2_i,               color='C0', ls='--', lw=1.2, label=r'$v_{t,i}^2$ self')
ax10.plot(z, vts2_e,               color='C1', ls='--', lw=1.2, label=r'$v_{t,e}^2$ self')
ax10.plot(z, vts2_cross_ion_gk,    color='C0', ls='-',  lw=2,   label=r'$v_{t,{\rm cross,ion}}^2$ Gkeyll')
ax10.plot(z, vts2_cross_formula_ion, color='k', ls=':',  lw=2,
          label=r'Formula ($s$=ion, $r$=elc)')
ax10.set_xlabel(r'$z$ (m)')
ax10.set_ylabel(r'$v_t^2$ (m$^2$ s$^{-2}$)')
ax10.legend(fontsize=8, frameon=False)
ax10.axhline(0, color='gray', lw=0.5, ls=':')
n_neg = int(np.sum(vts2_cross_formula_ion < 0))
if n_neg > 0:
    ax10.text(0.02, 0.05, f'Formula < 0 in {n_neg}/{len(z)} cells',
              transform=ax10.transAxes, fontsize=8, color='red')

ax11 = fig.add_subplot(gs[1, 1])
ax11.set_title(r'ELC cross vtsq:  formula ($s$=elc, $r$=ion) vs Gkeyll')
ax11.plot(z, vts2_e,               color='C1', ls='--', lw=1.2, label=r'$v_{t,e}^2$ self')
ax11.plot(z, vts2_i,               color='C0', ls='--', lw=1.2, label=r'$v_{t,i}^2$ self')
ax11.plot(z, vts2_cross_elc_gk,    color='C1', ls='-',  lw=2,   label=r'$v_{t,{\rm cross,elc}}^2$ Gkeyll')
ax11.plot(z, vts2_cross_formula_elc, color='k', ls=':',  lw=2,
          label=r'Formula ($s$=elc, $r$=ion)')
ax11.set_xlabel(r'$z$ (m)')
ax11.set_ylabel(r'$v_t^2$ (m$^2$ s$^{-2}$)')
ax11.legend(fontsize=8, frameon=False)
rel_err_elc = np.abs(vts2_cross_formula_elc - vts2_cross_elc_gk) / np.abs(vts2_cross_elc_gk)
ax11.text(0.02, 0.05, f'Max rel. err: {rel_err_elc.max():.2e}',
          transform=ax11.transAxes, fontsize=8, color='darkgreen')

# fig.suptitle(
#     r'Cross-species LBO primary moments  ($\delta$=2, $\beta$=1, $d_v$=2)'
#     '\n'
#     r'$v_{t,sr}^2 = v_{t,s}^2 + \frac{\delta(1+\beta)/2}{1+m_s/m_r}'
#     r'\!\left[v_{t,r}^2 - \frac{m_s}{m_r}v_{t,s}^2 + \frac{(u_s-u_r)^2}{d_v}\right]$',
#     fontsize=11, y=1.02)

outfile = dataDir + 'lbo-cross-moms-compare.png'
plt.show()
# plt.savefig(outfile, dpi=150, bbox_inches='tight')
print(f'Saved {outfile}')

# ---- Print summary at midplane ----
midz = len(z) // 2
print(f'\nAt z ~ {z[midz]:.3f} m  (cell {midz}):')
print(f'  Ion self:  u = {u_ion[midz]:.3g} m/s,   vtsq = {vts2_i[midz]:.4g} m^2/s^2')
print(f'  Elc self:  u = {u_elc[midz]:.3g} m/s,   vtsq = {vts2_e[midz]:.4g} m^2/s^2')
print()
print(f'  Gkeyll ion_cross_elc:   u = {u_cross_ion_gk[midz]:.3g},  vtsq = {vts2_cross_ion_gk[midz]:.4g}')
print(f'  Formula (s=ion,r=elc):  vtsq = {vts2_cross_formula_ion[midz]:.4g}  '
      f'[{"NEGATIVE" if vts2_cross_formula_ion[midz]<0 else "positive"}]')
print()
print(f'  Gkeyll elc_cross_ion:   u = {u_cross_elc_gk[midz]:.3g},  vtsq = {vts2_cross_elc_gk[midz]:.4g}')
print(f'  Formula (s=elc,r=ion):  vtsq = {vts2_cross_formula_elc[midz]:.4g}')
print(f'  Relative error (elc case, midplane): {rel_err_elc[midz]:.3e}')
print(f'  Max relative error (elc case):       {rel_err_elc.max():.3e}')
