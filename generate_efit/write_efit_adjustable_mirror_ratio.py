#!/usr/bin/env python3
"""
Generate EFIT/geqdsk files for a magnetic mirror configuration.

The magnetic field is produced by circular coils on a cylinder of radius
``COIL_RADIUS``:

  * Two **high-field** (mirror-throat) coils at  z = +/- COIL_DIST/2,
    each carrying current ``I_high``.
  * ``N_MIDDLE_COILS`` **low-field** coils uniformly distributed between
    the mirror throats, each carrying current ``I_low``.  These flatten
    the field profile near the device centre (z = 0).

The script optimises the two currents so that, on the flux surface
psi = PSI_EVAL:

    |B|_max   (at the mirror throats, z ~ +/- COIL_DIST/2)  =  TARGET_BMAX
    |B|_min   (at the device centre,  z = 0)                 =  TARGET_BMIN

giving mirror ratio  R_m  =  Bmax / Bmin.

The poloidal flux psi = R * A_phi is computed from the exact vector
potential of each circular current loop (Jackson, sec 5.5):

    A_phi = -(mu0 * I / pi) * sqrt(a/R) *
            [ (k^2 - 2)/(2k) * K(k^2)  +  E(k^2)/k ]

    k^2 = 4 a R / [ (a + R)^2 + (Z - Z_coil)^2 ]

where K and E are the complete elliptic integrals of the first and
second kind.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")                       # non-interactive backend
import matplotlib.pyplot as plt
import scipy.constants
import scipy.integrate as integrate
import scipy.optimize as sco
from scipy.interpolate import interp1d
from scipy.special import ellipk, ellipe
import fortranformat as ff
from datetime import date
from typing import Any, Iterable, TextIO


# ======================================================================
#  Configuration
# ======================================================================

# --- Mirror field targets -------------------------------------------------
MIRROR_RATIO = 5.0                          # desired Bmax / Bmin
TARGET_BMIN     = 0.5                          # [T]  target |B| at z = 0
TARGET_BMAX     = MIRROR_RATIO * TARGET_BMIN      # [T]  target |B| at throats

# --- Coil geometry --------------------------------------------------------
COIL_RADIUS    = 0.45       # [m]  radius of every coil  (must be > RMAX)
N_MIDDLE_COILS = 48         # number of low-field coils between the mirrors
COIL_DIST      = 2.0        # [m]  axial separation of mirror coils (= h)

# --- Flux surface for B evaluation ----------------------------------------
PSI_EVAL = 1e-3             # value of psi on which to measure Bmax / Bmin

# --- Background / toroidal field (zero for a pure mirror) -----------------
B0 = 0.0
R0 = 0.02

# --- RZ computational grid ------------------------------------------------
NW = 257                    # number of R grid points
NH = 257                    # number of Z grid points
RMIN, RMAX = 0.001, 0.4
ZMIN, ZMAX = -2.0,  2.0

Rgrid = np.linspace(RMIN, RMAX, NW)
Zgrid = np.linspace(ZMIN, ZMAX, NH)

# Derived quantities used in the geqdsk header
RDIM    = RMAX - RMIN
ZDIM    = ZMAX - ZMIN
RLEFT   = RMIN
ZMID    = 0.5 * (ZMAX + ZMIN)
RMAXIS  = 0.0
ZMAXIS  = 0.0
RCENTR  = 0.0
BCENTR  = B0
CURRENT = 0.0

# Physical constant
MU0 = scipy.constants.mu_0


# ======================================================================
#  Flux and field functions
# ======================================================================

def _aphi_single_coil(R, Z_rel, I, a=COIL_RADIUS):
    """Toroidal component of the vector potential from one circular coil.

    Parameters
    ----------
    R     : array_like -- cylindrical radius of the field point(s)
    Z_rel : array_like -- Z_field - Z_coil  (signed axial distance)
    I     : float      -- coil current  [A]
    a     : float      -- coil radius   [m]

    Returns
    -------
    Aphi : ndarray -- same shape as broadcast(R, Z_rel)
    """
    R     = np.asarray(R, dtype=float)
    Z_rel = np.asarray(Z_rel, dtype=float)
    R_safe = np.maximum(R, 1e-17)                   # avoid 1/0 on axis

    k2 = 4.0 * a * R_safe / ((a + R_safe)**2 + Z_rel**2)
    k  = np.sqrt(k2)

    Aphi = (-MU0 * I / np.pi) * np.sqrt(a / R_safe) * (
        (k2 - 2.0) / (2.0 * k) * ellipk(k2) + ellipe(k2) / k
    )
    return Aphi


def compute_psi(R, Z, I_high, I_low):
    """Poloidal flux  psi = R * A_phi  from all coils.

    Parameters
    ----------
    R, Z   : array_like -- field-point coordinates (broadcastable)
    I_high : float      -- current in each mirror-throat coil
    I_low  : float      -- current in each central coil
    """
    R = np.asarray(R, dtype=float)
    Z = np.asarray(Z, dtype=float)

    # Mirror-throat coils at z = +/- COIL_DIST/2
    Aphi  = _aphi_single_coil(R, Z - COIL_DIST / 2, I_high)
    Aphi += _aphi_single_coil(R, Z + COIL_DIST / 2, I_high)

    # Low-field coils uniformly spaced between the throats (excl. endpoints)
    z_coils = np.linspace(-COIL_DIST / 2, COIL_DIST / 2,
                          N_MIDDLE_COILS + 2)[1:-1]
    for zc in z_coils:
        Aphi += _aphi_single_coil(R, Z - zc, I_low)

    return np.maximum(R, 1e-17) * Aphi


def compute_fields_on_grid(I_high, I_low):
    """Evaluate psi, Bz, Br, |B| on the full (Rgrid x Zgrid) mesh.

    Uses  Bz = (1/R) dpsi/dR  and  Br = -(1/R) dpsi/dZ.

    Returns
    -------
    psiRZ, BzRZ, BrRZ, BRZ : ndarray of shape (NW, NH)
    """
    RR, ZZ = np.meshgrid(Rgrid, Zgrid, indexing="ij")   # shape (NW, NH)

    psiRZ = compute_psi(RR, ZZ, I_high, I_low)

    # Central finite differences via np.gradient
    dpsi_dR = np.gradient(psiRZ, Rgrid, axis=0)
    dpsi_dZ = np.gradient(psiRZ, Zgrid, axis=1)

    R_safe = np.maximum(RR, Rgrid[1])        # avoid 1/0 at smallest R
    BzRZ =  dpsi_dR / R_safe
    BrRZ = -dpsi_dZ / R_safe
    BRZ  = np.sqrt(BrRZ**2 + BzRZ**2)

    return psiRZ, BzRZ, BrRZ, BRZ


def B_along_flux_surface(psiRZ, BRZ, psi_target=PSI_EVAL):
    """Extract |B| along the flux surface psi = psi_target.

    For each Z slice, root-find in R to locate the flux surface, then
    interpolate |B| at that radius.

    Returns
    -------
    B_fs : ndarray of shape (NH,)
        |B| on the surface (NaN where the surface does not exist).
    """
    B_fs = np.full(NH, np.nan)
    for j in range(NH):
        psi_col = psiRZ[:, j]

        # Check that psi_target is bracketed by the R profile
        if psi_col.min() > psi_target or psi_col.max() < psi_target:
            continue

        psi_interp = interp1d(Rgrid, psi_col, kind="cubic")
        try:
            R_target = sco.brentq(lambda r: psi_interp(r) - psi_target,
                                  RMIN, RMAX)
        except ValueError:
            continue

        B_interp = interp1d(Rgrid, BRZ[:, j], kind="cubic")
        B_fs[j] = B_interp(R_target)

    return B_fs


# ======================================================================
#  Optimisation
# ======================================================================

def _residual(params):
    """Return [Bmax - AIM_BMAX, Bmin - AIM_BMIN] for given (I_high, I_low)."""
    I_high, I_low = params

    psiRZ, _, _, BRZ = compute_fields_on_grid(I_high, I_low)
    B_fs = B_along_flux_surface(psiRZ, BRZ)

    Bmax = np.nanmax(B_fs)
    Bmin = B_fs[NH // 2]                     # z = 0  (device centre)

    print(f"  I_high = {I_high:12.1f}   I_low = {I_low:12.1f}   "
          f"Bmax = {Bmax:.4f}   Bmin = {Bmin:.4f}   ratio = {Bmax / Bmin:.2f}")

    return [Bmax - TARGET_BMAX, Bmin - TARGET_BMIN]


def optimise_currents(I_high0=812_998.0, I_low0=0.0):
    """Find (I_high, I_low) that realise the target Bmax and Bmin.

    Uses scipy.optimize.fsolve (Powell hybrid method) to solve the 2x2
    nonlinear system.
    """
    print(f"Target:  Bmax = {TARGET_BMAX:.2f} T,  Bmin = {TARGET_BMIN:.2f} T  "
          f"(mirror ratio = {MIRROR_RATIO})\n")

    sol, info, ier, msg = sco.fsolve(_residual, [I_high0, I_low0],
                                     full_output=True)
    if ier != 1:
        print(f"WARNING: fsolve did not converge -- {msg}")

    I_high, I_low = sol
    print(f"\nOptimised currents:  I_high = {I_high:.1f} A,  "
          f"I_low = {I_low:.1f} A\n")
    return I_high, I_low


# ======================================================================
#  geqdsk (EFIT) file writer
# ======================================================================

def _write_fortran_line(data: Iterable[Any], fh: TextIO, fmt: str):
    """Write one Fortran-formatted record followed by a newline."""
    fh.write(ff.FortranRecordWriter(fmt).write(data))
    fh.write("\n")


def _write_1d(arr, fh, per_line=5):
    """Write a 1-D float array, ``per_line`` values per line."""
    for i, v in enumerate(arr):
        fh.write(f"{v:16.9E}")
        if (i + 1) % per_line == 0:
            fh.write("\n")
    if len(arr) % per_line != 0:
        fh.write("\n")


def write_geqdsk(filename, psiRZ, I_high, I_low):
    """Write the geqdsk (EFIT) file.

    Parameters
    ----------
    filename : str
    psiRZ    : ndarray (NW, NH)  -- poloidal flux on the RZ grid
    I_high, I_low : float        -- optimised coil currents
    """
    NPSI = NW

    # Boundary flux values (evaluated analytically outside the grid)
    SIMAG = float(compute_psi(2.0, 0.0, I_high, I_low))
    SIBRY = float(compute_psi(4.0, 0.0, I_high, I_low))

    # 1-D profile arrays on the psi grid
    PSIGRID = np.linspace(SIMAG, SIBRY, NPSI)
    FPOL    = np.full(NPSI, B0 * R0)                # F = R*Bphi  (zero)
    FFPRIM  = np.zeros(NPSI)
    PPRIME  = np.full(NPSI, -1e-6)
    PRES    = integrate.cumulative_trapezoid(PPRIME, PSIGRID, initial=0)
    QPSI    = np.zeros(NPSI)

    PSIZR = psiRZ.T                                  # geqdsk: (NH, NW)

    # 48-character header comment
    label = "GKYLMAX"
    creation_date = date.today().strftime("%d/%m/%Y")
    comment = f"{label:11s}{creation_date:10s}   {'# 0':>8s}{'  0ms':16s}"

    with open(filename, "w", newline="") as f:
        # Line 1:  header comment + idum + NW + NH
        _write_fortran_line((comment, 3, NW, NH), f, "(a48,3i4)")

        # Helper: write exactly 5 floats on one line
        def w5(*vals):
            for v in vals:
                f.write(f"{float(v):16.9E}")
            f.write("\n")

        w5(RDIM,    ZDIM,   RCENTR, RLEFT, ZMID)
        w5(RMAXIS,  ZMAXIS, SIMAG,  SIBRY, BCENTR)
        w5(CURRENT, SIMAG,  0.0,    RMAXIS, 0.0)
        w5(ZMAXIS,  0.0,    SIBRY,  0.0,    0.0)

        # 1-D profiles
        for arr in (FPOL, PRES, FFPRIM, PPRIME):
            _write_1d(arr, f)

        # 2-D flux  PSIZR[j, k]  (NH rows x NW columns)
        count = 0
        for j in range(NH):
            for k in range(NW):
                f.write(f"{PSIZR[j, k]:16.9E}")
                count += 1
                if count % 5 == 0:
                    f.write("\n")
        if count % 5 != 0:
            f.write("\n")

        # Safety factor
        _write_1d(QPSI, f)

    print(f"Wrote  {filename}")


# ======================================================================
#  Diagnostic plots
# ======================================================================

def make_plots(psiRZ, BRZ, B_fs, out_prefix="mirror_field"):
    """Generate and save a 2x2 panel of diagnostic plots."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # (a) |B| on the target flux surface
    ax = axes[0, 0]
    ax.plot(Zgrid, B_fs)
    ax.axhline(TARGET_BMAX, ls="--", color="r", lw=0.8,
               label=f"Bmax target = {TARGET_BMAX:.1f} T")
    ax.axhline(TARGET_BMIN, ls="--", color="b", lw=0.8,
               label=f"Bmin target = {TARGET_BMIN:.1f} T")
    ax.set_xlabel("Z  [m]")
    ax.set_ylabel("|B|  [T]")
    ax.set_title(f"|B| on flux surface  psi = {PSI_EVAL}")
    ax.legend(fontsize=8)

    # (b) Psi contours
    ax = axes[0, 1]
    cs = ax.contour(Rgrid, Zgrid, psiRZ.T, levels=50)
    ax.set_xlabel("R  [m]")
    ax.set_ylabel("Z  [m]")
    ax.set_title("psi(R, Z)")
    plt.colorbar(cs, ax=ax)

    # (c) log10 |B| heat-map
    ax = axes[1, 0]
    with np.errstate(divide="ignore"):
        logB = np.log10(BRZ.T)
    pcm = ax.pcolormesh(Rgrid, Zgrid, logB, shading="auto")
    ax.set_xlabel("R  [m]")
    ax.set_ylabel("Z  [m]")
    ax.set_title("log10 |B|")
    plt.colorbar(pcm, ax=ax)

    # (d) |B| vs Z at two radii
    ax = axes[1, 1]
    ax.plot(Zgrid, BRZ[0, :],       label=f"R = {Rgrid[0]:.4f} m")
    ax.plot(Zgrid, BRZ[NW // 2, :], label=f"R = {Rgrid[NW // 2]:.3f} m")
    ax.set_xlabel("Z  [m]")
    ax.set_ylabel("|B|  [T]")
    ax.set_title("|B|(Z) at fixed R")
    ax.legend(fontsize=8)

    plt.tight_layout()
    figname = f"{out_prefix}_diagnostics.png"
    plt.savefig(figname, dpi=150)
    print(f"Saved  {figname}")


# ======================================================================
#  Main
# ======================================================================

def main():
    # Clean up old output files
    for fname in os.listdir("."):
        if fname.endswith(".geqdsk") and fname.startswith("conditioned_coil_R"):
            os.remove(fname)

    # Step 1 -- optimise coil currents for target mirror ratio
    I_high, I_low = optimise_currents()

    # Step 2 -- compute final fields on the RZ grid
    psiRZ, BzRZ, BrRZ, BRZ = compute_fields_on_grid(I_high, I_low)
    B_fs = B_along_flux_surface(psiRZ, BRZ)

    Bmax  = np.nanmax(B_fs)
    Bmin  = B_fs[NH // 2]
    ratio = Bmax / Bmin

    print(f"Final:  Bmax = {Bmax:.4f} T,  Bmin = {Bmin:.4f} T,  "
          f"mirror ratio = {ratio:.2f}")

    # Step 3 -- write geqdsk file
    out_name = f"conditioned_coil_R{int(round(ratio))}.geqdsk"
    write_geqdsk(out_name, psiRZ, I_high, I_low)

    # Step 4 -- save diagnostic plots
    make_plots(psiRZ, BRZ, B_fs)


if __name__ == "__main__":
    main()
