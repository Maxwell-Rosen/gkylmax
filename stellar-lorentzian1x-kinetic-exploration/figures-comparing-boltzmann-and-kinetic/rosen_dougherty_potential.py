"""Rosen-Dougherty theoretical potential calculations for electron confinement.

This module computes the theoretical e*phi/Te value for mirror/barrier systems
using the Rosen-Dougherty confinement time formulation.

Uses the same postgkyl loading framework as calc_pastukhov_potential.py.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import postgkyl as pg
from scipy import optimize, special
from scipy.integrate import quad
from scipy.special import erf


# Physical constants
ELECTRON_MASS = 9.1093837015e-31
ELEMENTARY_CHARGE = 1.602176634e-19

# Coordinates for averages/sampling
Z_CENTER = 0.0
Z_MIRROR_THROAT = 0.98
Z_SHEATH = 2.5


def find_map_file(run_path: Path, sim_prefix: str) -> str | None:
    """Find coordinate map file in order of preference."""
    candidates = [
        run_path / f"{sim_prefix}-mc2nu_pos_deflated.gkyl",
        run_path / f"{sim_prefix}-mc2nu_pos.gkyl",
        run_path / f"{sim_prefix}-mapc2p.gkyl",
    ]
    for path in candidates:
        if path.exists():
            return str(path)
    return None


def load_magnetic_field_ratio(run_path: Path, sim_prefix: str) -> float:
    """Load mirror ratio R = B_throat / B_0 from magnetic field file.
    
    Uses postgkyl pattern: GData -> interpolate(comp, overwrite=True) -> select(z0)
    """
    bmag_file = str(run_path / f"{sim_prefix}-bmag.gkyl")
    map_file = find_map_file(run_path, sim_prefix)

    # Load at z=0
    bmag_data = pg.GData(bmag_file, mapc2p_name=map_file)
    pg.GInterpModal(bmag_data).interpolate(0, overwrite=True)
    _, bmag_sel = pg.data.select(bmag_data, z0=f"{Z_CENTER}")
    bmag0 = float(np.squeeze(bmag_sel))

    # Load at throat
    bmag_data = pg.GData(bmag_file, mapc2p_name=map_file)
    pg.GInterpModal(bmag_data).interpolate(0, overwrite=True)
    _, bmag_sel = pg.data.select(bmag_data, z0=f"{Z_MIRROR_THROAT}")
    bmag_throat = float(np.squeeze(bmag_sel))

    R = bmag_throat / bmag0
    return R


def compute_electron_collision_freq(
    density_at_0: float, te0_ev: float
) -> float:
    """Compute electron-electron collision frequency.
    
    Uses Coulomb collision formula from Najmabadi (1984).
    
    Parameters
    ----------
    density_at_0 : float
        Electron density at z=0 in SI units [m^-3]
    te0_ev : float
        Electron temperature at z=0 in eV
    
    Returns
    -------
    float
        Collision frequency in SI units [s^-1]
    """
    eps0 = 8.8541878176204e-12

    n0_cm3 = density_at_0 * 1e-6
    log_lambda = 24 - np.log(n0_cm3**0.5 / te0_ev)

    te0_j = te0_ev * ELEMENTARY_CHARGE
    nu_ec = (
        log_lambda
        * density_at_0
        * (4 * np.pi / (2**1.5 * te0_j**1.5 * ELECTRON_MASS**0.5))
        * (ELEMENTARY_CHARGE**2 / (4 * np.pi * eps0)) ** 2
    )
    return nu_ec


def compute_particle_confinement_time(
    run_path: Path, sim_prefix: str, frame: int, z_trap_limit: float = 0.98
) -> float:
    """Compute particle confinement time tau_pi = int(M0 dx) / int(source_M0 dx).
    
    Parameters
    ----------
    run_path : Path
        Path to simulation data
    sim_prefix : str
        Simulation prefix (e.g., 'gk_lorentzian_mirror')
    frame : int
        Frame number for distribution
    z_trap_limit : float
        Z limit for integration region (default mirror throat at 0.98)
    
    Returns
    -------
    float
        Confinement time in seconds
    
    Uses postgkyl pattern: GData -> interpolate -> select(z0=range) -> integrate
    """
    map_file = find_map_file(run_path, sim_prefix)

    # Load M0 (density) for the frame
    m0_file = str(run_path / f"{sim_prefix}-ion_M0_{frame}.gkyl")
    m0_data = pg.GData(m0_file, mapc2p_name=map_file)
    pg.GInterpModal(m0_data).interpolate(comp=0, overwrite=True)
    pg.data.select(
        m0_data, z0=f"{-z_trap_limit}:{z_trap_limit}", overwrite=True
    )
    _, int_m0_dx = pg.tools.integrate(m0_data, 0, overwrite=True)
    int_m0 = float(np.squeeze(int_m0_dx))

    # Load source_M0 (source) for reference frame (usually 0)
    m0src_file = str(run_path / f"{sim_prefix}-ion_source_M0_0.gkyl")
    m0src_data = pg.GData(m0src_file, mapc2p_name=map_file)
    pg.GInterpModal(m0src_data).interpolate(comp=0, overwrite=True)
    pg.data.select(
        m0src_data, z0=f"{-z_trap_limit}:{z_trap_limit}", overwrite=True
    )
    _, int_m0src_dx = pg.tools.integrate(m0src_data, 0, overwrite=True)
    int_m0src = float(np.squeeze(int_m0src_dx))

    if int_m0src <= 0.0:
        raise ValueError(f"Invalid source integral: {int_m0src}")

    tau_pi = int_m0 / int_m0src
    return tau_pi


def rosen_dougherty_confinement_time(
    ephi_over_te: float, R: float, nu_e: float, Zpfl: float = 1.0, coeff: float = 0.0
) -> float:
    """Compute Rosen-Dougherty confinement time for given e*phi/Te.
    
    This is the full confinement time including the 1/nu_e factor, matching
    the reference formulation in calc_pastukhov_potential.py.
    
    Reference formula (from calc_pastukhov_potential.py):
    Loss_Rosen = (1/nu_e) * 1 / (2*Zpfl / (log(...) - coeff) * (1-erf(...)))
    
    Which simplifies to:
    Loss_Rosen = (1/nu_e) * (log(...) - coeff) / (2*Zpfl * (1-erf(...)))
    
    Parameters
    ----------
    ephi_over_te : float
        Normalized potential e*phi/Te (dimensionless)
    R : float
        Mirror ratio (B_throat / B_0)
    nu_e : float
        Electron collision frequency [s^-1]
    Zpfl : float
        Charge parameter (default 1.0 for electrons only)
    coeff : float
        Coefficient in log term (default 0.0 for pure Rosen-Dougherty)
    
    Returns
    -------
    float
        Confinement time in seconds
    """
    w_term = np.sqrt(1.0 + 2.0 * ephi_over_te / (R * Zpfl))
    a_term = np.sqrt(ephi_over_te + np.log(w_term))
    
    # CORRECTED formula: (1-erf) goes in DENOMINATOR, not numerator
    loss_rosen = (1.0 / nu_e) * (
        (np.log((w_term + 1.0) / (w_term - 1.0)) - coeff)
        / (2.0 * Zpfl * (1.0 - erf(a_term)))
    )
    return loss_rosen


def compute_theoretical_potential(
    run_path: Path,
    sim_prefix: str,
    frame: int,
    te0_ev: float,
    zpfl: float = 1.0,
    coeff: float = 1.117,
) -> tuple[float, float, float]:
    """Compute theoretical e*phi/Te using Rosen-Dougherty formulation.
    
    Solves: tau_pi = rosen_dougherty_confinement_time(P, R, nu_e, Zpfl, coeff)
    for P = e*phi/Te.
    
    Parameters
    ----------
    run_path : Path
        Path to simulation directory
    sim_prefix : str
        Simulation prefix (e.g., 'gk_lorentzian_mirror')
    frame : int
        Frame number for kinetic data
    te0_ev : float
        Electron temperature at z=0 in eV
    zpfl : float
        Charge parameter (default 1.0)
    coeff : float
        Coefficient in log term (default 1.117 from Dougherty)
    
    Returns
    -------
    tuple[float, float, float]
        (ephi_over_Te_theory, R, nu_e)
    
    Raises
    ------
    ValueError
        If root finding fails to find a solution
    """
    # Load magnetic field ratio
    R = load_magnetic_field_ratio(run_path, sim_prefix)

    # For Boltzmann electrons, use electron density from ion M0 (quasineutrality)
    map_file = find_map_file(run_path, sim_prefix)
    ion_m0_file = str(run_path / f"{sim_prefix}-ion_M0_{frame}.gkyl")
    
    m0_data = pg.GData(ion_m0_file, mapc2p_name=map_file)
    pg.GInterpModal(m0_data).interpolate(comp=0, overwrite=True)
    _, density_sel = pg.data.select(m0_data, z0=f"{Z_CENTER}")
    n0 = float(np.squeeze(density_sel))

    nu_e = compute_electron_collision_freq(n0, te0_ev)

    # Compute particle confinement time
    tau_pi = compute_particle_confinement_time(run_path, sim_prefix, frame)

    # Define root equation: tau_pi = rosen_dougherty_confinement_time(P, R, nu_e, Zpfl, coeff)
    def root_eq(ephi_over_te: float) -> float:
        rd_time = rosen_dougherty_confinement_time(ephi_over_te, R, nu_e, zpfl, coeff)
        return tau_pi - rd_time

    # Try root finding with increasing bracket sizes
    bracket_upper_limits = [20.0, 100.0, 500.0, 1000.0]
    ephi_over_te_theory = None
    
    for upper_limit in bracket_upper_limits:
        try:
            ephi_over_te_theory = optimize.ridder(
                root_eq, zpfl, upper_limit, maxiter=100
            )
            break  # Success!
        except ValueError:
            continue  # Try next bracket
    
    if ephi_over_te_theory is None:
        # Try brent method as fallback
        try:
            from scipy.optimize import brent
            ephi_over_te_theory = brent(root_eq, brack=(float(zpfl), 100.0))
        except Exception as e:
            raise ValueError(
                f"Failed to find root for e*phi/Te. Check tau_pi={tau_pi}, R={R}, nu_e={nu_e}. "
                f"Last error: {e}"
            ) from e

    return ephi_over_te_theory, R, nu_e
