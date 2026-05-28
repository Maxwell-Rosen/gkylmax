"""
source_matching_2x2v_boltzmann.py
==================================
Source-matching analysis for gk_lorentzian_mirror (2x2v, Boltzmann electrons).

Electrons are adiabatic (Boltzmann response). Only the ion species is kinetic;
all source diagnostics reflect ions only.

The simulation uses a non-uniform configuration-space grid. Physical (psi, Z)
coordinates are extracted from the mc2nu coordinate-mapping file:
  mc2nu[..., 0]  →  physical psi values  (component 0)
  mc2nu[..., 1]  →  physical Z   values  (component 1)

Source diagnostics loaded:
  ion_source_M2  →  M_2^{src}(psi,Z) = ∫ v² S_ion dv³  [m⁻³ s⁻¹ (m/s)²]
  ion_source_M0  →  M_0^{src}(psi,Z) = ∫    S_ion dv³  [m⁻³ s⁻¹]

Quantities computed:
  dP/dpsi      – ion source power per unit psi (Z-integrated with J_c)
  dΓ/dpsi      – ion source particle rate per unit psi
  psi_0        – power-weighted median flux surface  (Option 1 criterion)
  Δpsi_eff     – P_tot / (dP/dpsi)|_{psi_0}          (power effective width)
  Δpsi_eff_N   – Γ_tot / (dΓ/dpsi)|_{psi_0}          (particle effective width)
  Γ_1x/Γ_2x   – Δpsi_eff_N / Δpsi_eff               (source factorizability diagnostic)
  kappa        – (dP/dpsi)|_{psi_0,2x} / (dP/dpsi)|_{psi_0,1x,raw}

Reference: source-matching framework, Section "Source matching across
dimensionalities" of the thesis.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

from postgkyl.data import GData, GInterpModal

_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

# ── 0. USER CONFIGURATION ─────────────────────────────────────────────────────

SIM_DIR    = "/global/homes/m/mhrosen/scratch/gkylmax/perlmutter-lorentzian2x-orbit-average/R32-without-nuie/"   # directory containing .gkyl files
PREFIX     = "gk_lorentzian_mirror"
TF         = 0           # time frame for source diagnostics (0 = time independent source)

SIM_DIR_1X = "/global/homes/m/mhrosen/scratch/gkylmax/perlmutter-lorentzian2x-orbit-average/R32-without-nuie-1x/"   # 1x2v simulation directory
PREFIX_1X  = "gk_lorentzian_mirror"
TF_1X      = 0

POLY_ORDER = 1
BASIS_TYPE = "ms"

# Ion mass — only needed to convert M_2^{src} [m⁻³ s⁻¹ (m/s)²] → power density [W m⁻³]
M_ION = 2.014 * 1.673e-27   # e.g. deuterium [kg]

# ── 1. FILE PREFIXES ──────────────────────────────────────────────────────────

fp    = os.path.join(SIM_DIR,    PREFIX)
fp_1x = os.path.join(SIM_DIR_1X, PREFIX_1X)

# ── 2. HELPERS ────────────────────────────────────────────────────────────────

def load_physical_grids(mc2nu_file):
    """
    Extract physical (psi, Z) coordinates from the mc2nu grid-mapping file.

    mc2nu maps the computational uniform grid to the physical non-uniform
    configuration-space grid. Calling interpolate with all component indices
    returns a values array of shape (Npsi, NZ, cdim):
      [:, 0, 0]  →  physical psi at each psi index  (Z index irrelevant)
      [0, :, 1]  →  physical Z   at each Z   index  (psi index irrelevant)
    """
    gdata = GData(mc2nu_file)
    cdim  = gdata.get_num_dims()
    dg    = GInterpModal(gdata, POLY_ORDER, BASIS_TYPE)
    _, vals = dg.interpolate(tuple(range(cdim)))   # (Npsi, NZ, cdim)
    psi_phys = vals[:, 0, 0]   # physical psi coordinate, shape (Npsi,)
    Z_phys   = vals[0, :, 1]   # physical Z   coordinate, shape (NZ,)
    return psi_phys, Z_phys


def load_2x_values(filename, comp=0):
    """
    Load one DG component from a 2x configuration-space .gkyl file.
    Returns an array of shape (Npsi, NZ), squeezing the dummy y-dimension.
    Values are in raw (SI) units – no normalization applied.
    """
    gdata = GData(filename)
    dg    = GInterpModal(gdata, POLY_ORDER, BASIS_TYPE)
    dg.interpolate(comp, overwrite=True)
    vals  = gdata.get_values()[..., 0]   # drop component axis
    return np.squeeze(vals)              # (Npsi, NZ)


def load_1x_Z_grid(mc2nu_file):
    """
    Extract the physical Z coordinate from a 1x mc2nu grid-mapping file.

    For a 1x simulation cdim=1, so mc2nu stores one component: the physical Z
    values along the single computational axis. Shape after interpolate: (NZ, 1).
    """
    gdata = GData(mc2nu_file)
    dg    = GInterpModal(gdata, POLY_ORDER, BASIS_TYPE)
    _, vals = dg.interpolate((0,))   # (NZ, 1)
    return vals[:, 0]                # physical Z, shape (NZ_1x,)


def load_1x_values(filename, comp=0):
    """
    Load one DG component from a 1x configuration-space .gkyl file.
    Returns a 1-D array of shape (NZ,) in raw (SI) units.
    """
    gdata = GData(filename)
    dg    = GInterpModal(gdata, POLY_ORDER, BASIS_TYPE)
    dg.interpolate(comp, overwrite=True)
    vals  = gdata.get_values()[..., 0]
    return np.squeeze(vals)          # (NZ,)


# ── 3. LOAD PHYSICAL GRID FROM mc2nu ──────────────────────────────────────────

mc2nu_file = f"{fp}-mc2nu_pos_deflated.gkyl"
print(f"Loading physical grid from: {mc2nu_file}")
psi, Z = load_physical_grids(mc2nu_file)
Npsi, NZ = len(psi), len(Z)
print(f"  Grid:  Npsi={Npsi},  NZ={NZ}")
print(f"  psi ∈ [{psi[0]:.6f},  {psi[-1]:.6f}]")
print(f"  Z   ∈ [{Z[0]:.4f},  {Z[-1]:.4f}]")

# ── 4. LOAD SOURCE DIAGNOSTICS ────────────────────────────────────────────────

print(f"\nLoading source diagnostics at frame {TF}...")

# Second velocity moment of the ion source: M_2^{src} = ∫ v² S_ion dv³  [m⁻³ s⁻¹ (m/s)²]
M2_src = load_2x_values(f"{fp}-ion_source_M2_{TF}.gkyl")

# Power density:  P_ion(psi, Z) = (1/2) m_ion M_2^{src}   [W m⁻³]
P_density = 0.5 * M_ION * M2_src

# Zeroth velocity moment of the ion source: M_0^{src} = ∫ S_ion dv³   [m⁻³ s⁻¹]
M0_src = load_2x_values(f"{fp}-ion_source_M0_{TF}.gkyl")

# Configuration-space Jacobian J_c = √(det g)
Jac = load_2x_values(f"{fp}-jacobgeo.gkyl")

print(f"  P_density   shape: {P_density.shape}")
print(f"  M0_src      shape: {M0_src.shape}")
print(f"  Jacobian    shape: {Jac.shape}")

# ── 5. dP/dpsi  and  dΓ/dpsi ─────────────────────────────────────────────────
#
#   dP/dpsi(psi)  = 2π ∫_Z  P_ion(psi,Z) · J_c(psi,Z)  dZ
#   dΓ/dpsi(psi)  = 2π ∫_Z  M_0^{src}(psi,Z) · J_c(psi,Z)  dZ
#
# The 2π factor comes from integrating out the toroidal angle α under
# axisymmetry:  ∫ dα = 2π.

dP_dpsi     = 2*np.pi * _trapz(P_density * Jac, x=Z, axis=1)   # (Npsi,)
dGamma_dpsi = 2*np.pi * _trapz(M0_src    * Jac, x=Z, axis=1)   # (Npsi,)

P_tot     = _trapz(dP_dpsi,     x=psi)
Gamma_tot = _trapz(dGamma_dpsi, x=psi)

# ── 6. PSI_0: POWER-WEIGHTED MEDIAN ──────────────────────────────────────────
#
#   Choose psi_0 s.t.  ∫_{psi_min}^{psi_0} dP/dpsi dpsi = (1/2) P_tot

cumP = scipy.integrate.cumulative_trapezoid(dP_dpsi, x=psi, initial=0.0)

half = 0.5 * P_tot
i0   = np.clip(int(np.searchsorted(cumP, half)) - 1, 0, Npsi - 2)
frac = (half - cumP[i0]) / (cumP[i0 + 1] - cumP[i0])
psi_0 = psi[i0] + frac * (psi[i0 + 1] - psi[i0])

dPdpsi_psi0     = float(np.interp(psi_0, psi, dP_dpsi))
dGammadpsi_psi0 = float(np.interp(psi_0, psi, dGamma_dpsi))

# ── 7. EFFECTIVE WIDTHS ───────────────────────────────────────────────────────
#
#   Δpsi_eff   = P_tot   / (dP/dpsi)|_{psi_0}
#   Δpsi_eff_N = Γ_tot   / (dΓ/dpsi)|_{psi_0}
#
#   The ratio Δpsi_eff_N / Δpsi_eff = Γ_1x / Γ_2x measures source
#   non-factorizability: it equals 1 when S ∝ g(psi) h(Z, v).

Delta_psi_eff   = P_tot     / dPdpsi_psi0
Delta_psi_eff_N = Gamma_tot / dGammadpsi_psi0
mismatch_ratio  = Delta_psi_eff_N / Delta_psi_eff

# ── 8. LOAD 1x SOURCE DATA ───────────────────────────────────────────────────
#
#   The 1x simulation runs along Z at a fixed psi_0. Its source shape in Z
#   is inherited from the 2x source evaluated at psi_0 (with kappa=1).
#   Here we load the 1x source diagnostics and compute the Z-integrated
#   power and particle rates directly (no psi integral — the 1x domain is
#   a single field line).

print(f"\nLoading 1x source diagnostics at frame {TF_1X}...")
Z_1x      = load_1x_Z_grid(f"{fp_1x}-mc2nu_pos_deflated.gkyl")
M2_src_1x = load_1x_values(f"{fp_1x}-ion_source_M2_{TF_1X}.gkyl")
M0_src_1x = load_1x_values(f"{fp_1x}-ion_source_M0_{TF_1X}.gkyl")
Jac_1x    = load_1x_values(f"{fp_1x}-jacobgeo.gkyl")

P_density_1x = 0.5 * M_ION * M2_src_1x

dPdpsi_1x_raw     = 2*np.pi * _trapz(P_density_1x * Jac_1x, x=Z_1x)
dGammadpsi_1x_raw = 2*np.pi * _trapz(M0_src_1x    * Jac_1x, x=Z_1x)

print(f"  Z_1x        range: [{Z_1x[0]:.4f}, {Z_1x[-1]:.4f}],  NZ={len(Z_1x)}")
print(f"  dP/dpsi (1x,raw)  = {dPdpsi_1x_raw:.4e}  W/[psi]")
print(f"  dΓ/dpsi (1x,raw)  = {dGammadpsi_1x_raw:.4e}  s⁻¹/[psi]")

# ── 9. KAPPA ──────────────────────────────────────────────────────────────────
#
#   kappa   = (dP/dpsi)|_{psi_0, 2x} / (dP/dpsi)|_{psi_0, 1x, raw}
#   kappa_N = (dΓ/dpsi)|_{psi_0, 2x} / (dΓ/dpsi)|_{psi_0, 1x, raw}
#
#   "raw" means the 1x source was constructed with kappa=1 (no rescaling),
#   inheriting the 2x source shape at psi_0. kappa is then the factor by
#   which the 1x source must be multiplied to match the 2x power at psi_0.

kappa   = dPdpsi_psi0     / dPdpsi_1x_raw
kappa_N = dGammadpsi_psi0 / dGammadpsi_1x_raw

# ── 10. SUMMARY ───────────────────────────────────────────────────────────────

print()
print("=" * 62)
print("  SOURCE MATCHING  —  gk_lorentzian_mirror  (2x2v, Boltz. e⁻)")
print("=" * 62)
print(f"  P_tot (2x2v)               = {P_tot:.4e}  W")
print(f"  Γ_tot (2x2v)               = {Gamma_tot:.4e}  s⁻¹")
print()
print(f"  psi_0  (power median)      = {psi_0:.6f}")
print(f"  dP/dpsi  |_{{psi_0}}         = {dPdpsi_psi0:.4e}  W/[psi]")
print(f"  dΓ/dpsi  |_{{psi_0}}         = {dGammadpsi_psi0:.4e}  s⁻¹/[psi]")
print()
print(f"  Δpsi_eff  (power)          = {Delta_psi_eff:.6f}  [psi]")
print(f"  Δpsi_eff  (particle)       = {Delta_psi_eff_N:.6f}  [psi]")
print(f"  Γ_1x/Γ_2x mismatch ratio   = {mismatch_ratio:.4f}",
      " (1 → factorizable source)" if abs(mismatch_ratio - 1) < 0.05 else "")
print()
print(f"  kappa   (power)            = {kappa:.6f}")
print(f"  kappa_N (particle)         = {kappa_N:.6f}")
print("=" * 62)

# ── 10. PLOTS ─────────────────────────────────────────────────────────────────

cumGamma = scipy.integrate.cumulative_trapezoid(dGamma_dpsi, x=psi, initial=0.0)
PSI2D, Z2D = np.meshgrid(psi, Z, indexing="ij")

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle(
    "Source matching — gk_lorentzian_mirror  (2x2v, Boltzmann e⁻)",
    fontsize=12, y=1.01,
)

# ── Row 0: POWER ──────────────────────────────────────────────────────────────

ax = axes[0, 0]
ax.plot(psi, dP_dpsi * 1e-3, "k-", lw=1.5)
ax.axvline(psi_0, color="tab:red", ls="--", lw=1.5,
           label=rf"$\psi_0 = {psi_0:.4f}$")
ax.set_xlabel(r"$\psi$")
ax.set_ylabel(r"$dP/d\psi$  [kW / $\psi$-unit]")
ax.set_title("Ion source power profile")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
ax.plot(psi, cumP / P_tot, "k-", lw=1.5)
ax.axhline(0.5, color="gray", ls=":", lw=1, label="50 %")
ax.axvline(psi_0, color="tab:red", ls="--", lw=1.5,
           label=rf"$\psi_0 = {psi_0:.4f}$")
ax.set_xlabel(r"$\psi$")
ax.set_ylabel(r"$P(<\!\psi) \;/\; P_{\rm tot}$")
ax.set_title("Cumulative power fraction")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

ax = axes[0, 2]
pc = ax.pcolormesh(PSI2D, Z2D, P_density * 1e-6, cmap="inferno", shading="auto")
ax.axvline(psi_0, color="w", ls="--", lw=1.5, label=r"$\psi_0$")
ax.set_xlabel(r"$\psi$")
ax.set_ylabel(r"$Z$")
ax.set_title(r"Ion source power density  $\frac{1}{2}m_i M_2^{\rm src}$")
plt.colorbar(pc, ax=ax, label=r"MW m$^{-3}$")
ax.legend(fontsize=8, loc="upper right")

# ── Row 1: PARTICLES ──────────────────────────────────────────────────────────

ax = axes[1, 0]
ax.plot(psi, dGamma_dpsi, "k-", lw=1.5)
ax.axvline(psi_0, color="tab:red", ls="--", lw=1.5,
           label=rf"$\psi_0 = {psi_0:.4f}$")
ax.set_xlabel(r"$\psi$")
ax.set_ylabel(r"$d\Gamma/d\psi$  [s$^{-1}$ / $\psi$-unit]")
ax.set_title("Ion source particle rate profile")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
ax.plot(psi, cumP     / P_tot,     "k--", lw=1.2, alpha=0.7, label=r"$P$ (power)")
ax.plot(psi, cumGamma / Gamma_tot, "k-",  lw=1.5,            label=r"$\Gamma$ (particle)")
ax.axhline(0.5, color="gray", ls=":", lw=1)
ax.axvline(psi_0, color="tab:red", ls="--", lw=1.5,
           label=rf"$\psi_0 = {psi_0:.4f}$")
ax.set_xlabel(r"$\psi$")
ax.set_ylabel(r"Cumulative fraction")
ax.set_title(
    rf"Non-factorizability:  $\Gamma_{{1x}}/\Gamma_{{2x}} = {mismatch_ratio:.3f}$"
)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

ax = axes[1, 2]
pc = ax.pcolormesh(PSI2D, Z2D, M0_src, cmap="plasma", shading="auto")
ax.axvline(psi_0, color="w", ls="--", lw=1.5, label=r"$\psi_0$")
ax.set_xlabel(r"$\psi$")
ax.set_ylabel(r"$Z$")
ax.set_title(r"Ion source particle density  $M_0^{\rm src}$")
plt.colorbar(pc, ax=ax, label=r"m$^{-3}$ s$^{-1}$")
ax.legend(fontsize=8, loc="upper right")

fig.tight_layout()
out_fig = "source_matching_2x2v_boltzmann.png"
fig.savefig(out_fig, dpi=150, bbox_inches="tight")
print(f"\nFigure saved → {out_fig}")
plt.show()
