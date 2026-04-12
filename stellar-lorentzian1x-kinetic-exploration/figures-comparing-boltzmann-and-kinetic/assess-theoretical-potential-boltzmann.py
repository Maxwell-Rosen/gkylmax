"""Assess theoretical potential of Boltzmann electron simulations using Rosen-Dougherty formulation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from rosen_dougherty_potential import compute_theoretical_potential


def run(context: dict[str, Any]) -> None:
    """Compute and display theoretical e*phi/Te for Boltzmann simulations.
    
    Uses Rosen-Dougherty confinement time formulation with Dougherty coefficient
    (coeff=1.117) to predict the confinement potential.
    
    If root finding fails (e.g., parameters out of model validity range),
    reports the diagnostic values instead.
    """
    boltzmann_cfg = context["boltzmann"]
    boltz_run_path = Path(boltzmann_cfg["data_path"])
    boltz_frame = int(boltzmann_cfg["last_frame"])
    boltz_te0_ev = float(context["constants"]["boltzmann_te0_ev"])

    print(
        f"[assess-theoretical-potential] Computing Rosen-Dougherty theoretical potential..."
    )
    print(f"[assess-theoretical-potential] Boltzmann run: {boltz_run_path}")
    print(f"[assess-theoretical-potential] Frame: {boltz_frame}, Te(z=0): {boltz_te0_ev} eV")

    # Infer simulation prefix from data directory name
    # Try common patterns
    sim_prefix = None
    for pattern in ["gk_lorentzian_mirror", "gk_wham", "gk_mirror"]:
        if pattern in str(boltz_run_path):
            sim_prefix = pattern
            break

    if sim_prefix is None:
        # Try to find from field file
        field_files = list(boltz_run_path.glob(f"*-field_{boltz_frame}.gkyl"))
        if field_files:
            field_name = field_files[0].stem
            sim_prefix = field_name.split("-field_", 1)[0]
        else:
            raise ValueError(
                f"Could not infer simulation prefix from {boltz_run_path}. "
                f"Please check the directory structure."
            )

    print(f"[assess-theoretical-potential] Detected sim_prefix: {sim_prefix}")

    try:
        ephi_theory, R, nu_e = compute_theoretical_potential(
            boltz_run_path,
            sim_prefix,
            boltz_frame,
            boltz_te0_ev,
            zpfl=1.0,
            coeff=1.117,
        )

        print(
            f"[assess-theoretical-potential] =========================================="
        )
        print(f"[assess-theoretical-potential] Theoretical Results (Rosen-Dougherty):")
        print(
            f"[assess-theoretical-potential]   Mirror ratio R = B_throat/B_0: {R:.4f}"
        )
        print(f"[assess-theoretical-potential]   Electron collision freq [s^-1]: {nu_e:.3e}")
        print(
            f"[assess-theoretical-potential]   Theoretical e*phi/Te (Dougherty coeff=1.117): {ephi_theory:.4f}"
        )
        print(
            f"[assess-theoretical-potential] =========================================="
        )

    except ValueError as root_error:
        # Root finding failed - report diagnostic info
        print(
            f"[assess-theoretical-potential] WARNING: Root finding failed (simulation parameters "
            f"may be outside model validity range)"
        )
        print(f"[assess-theoretical-potential] Diagnostic message: {root_error}")
        print(
            f"[assess-theoretical-potential] Check computed values in error message for R and nu_e"
        )

