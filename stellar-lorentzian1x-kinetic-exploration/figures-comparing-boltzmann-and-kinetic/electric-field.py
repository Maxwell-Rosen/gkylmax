"""Plot physical electric field overlay (V/m) for Boltzmann vs kinetic runs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from field_plot_common import (
    BOLTZMANN_COLOR,
    KINETIC_COLOR,
    Z_CENTER,
    Z_MIRROR_THROAT,
    Z_SHEATH,
    prepare_profiles,
)


def run(context: dict[str, Any]) -> None:
    out_dir = Path(context["figures_dir"]) / "field-and-potential"
    out_dir.mkdir(parents=True, exist_ok=True)

    data = prepare_profiles(context)
    z_kin = data["z_kin"]
    phi_kin = data["phi_kin"]
    z_boltz = data["z_boltz"]
    phi_boltz = data["phi_boltz"]

    e_kin = -np.gradient(phi_kin, z_kin)
    e_boltz = -np.gradient(phi_boltz, z_boltz)

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.plot(z_kin, e_kin, color=KINETIC_COLOR, lw=2.0, label="Kinetic electrons")
    ax.plot(z_boltz, e_boltz, color=BOLTZMANN_COLOR, lw=2.0, label="Boltzmann electrons")

    for marker_x in (Z_CENTER, Z_MIRROR_THROAT, Z_SHEATH):
        ax.axvline(marker_x, color="0.35", lw=1.0, ls="--")

    ax.set_xlabel(r"$z$ [m]")
    ax.set_ylabel(r"$E_z = -\nabla \phi$ [V/m]")
    ax.set_title("Electric field from potential gradient")
    ax.grid(alpha=0.3)
    ax.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_dir / "electric-field-vm-overlay-kinetic-vs-boltzmann.pdf")
    plt.close(fig)

    print("[electric-field] Saved V/m figure to", out_dir)


def main(context: dict[str, Any] | None = None) -> None:
    if context is None:
        raise ValueError("This subfile is intended to be called by master.py with context.")
    run(context)
