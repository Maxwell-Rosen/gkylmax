"""Plot normalized electric-field proxy overlay for Boltzmann vs kinetic runs."""

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
    phi_kin_norm = data["phi_kin_norm"]
    z_boltz = data["z_boltz"]
    phi_boltz_norm = data["phi_boltz_norm"]

    e_kin_norm = -np.gradient(phi_kin_norm, z_kin)
    e_boltz_norm = -np.gradient(phi_boltz_norm, z_boltz)

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.plot(z_kin, e_kin_norm, color=KINETIC_COLOR, lw=2.0, label="Kinetic electrons")
    ax.plot(z_boltz, e_boltz_norm, color=BOLTZMANN_COLOR, lw=2.0, label="Boltzmann electrons")

    for marker_x in (Z_CENTER, Z_MIRROR_THROAT, Z_SHEATH):
        ax.axvline(marker_x, color="0.35", lw=1.0, ls="--")

    for marker_x, marker_label in (
        (Z_CENTER, "center"),
        (Z_MIRROR_THROAT, "mirror throat"),
    ):
        ax.axvline(marker_x, color="0.35", lw=1.0, ls="--")
        ax.text(
            marker_x+0.1,
            0.9,
            marker_label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            rotation=90,
            fontsize=12,
        )

    for marker_x, marker_label in (
    (Z_SHEATH, "sheath"),
    ):
        ax.axvline(marker_x, color="0.35", lw=1.0, ls="--")
        ax.text(
            marker_x-0.1,
            0.9,
            marker_label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            rotation=90,
            fontsize=12,
        )

    ax.set_xlabel(r"$z$ [m]")
    ax.set_ylabel(r"$-\nabla (e\phi/T_e)$ [m$^{-1}$]")
    ax.set_xlim(-Z_SHEATH, Z_SHEATH)
    # ax.set_title("Electric field proxy from normalized potential")
    # ax.grid(alpha=0.3)
    ax.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_dir / "electric-field-normalized-overlay-kinetic-vs-boltzmann.pdf")
    plt.close(fig)

    print("[electric-field-normalized] Saved normalized figure to", out_dir)


def main(context: dict[str, Any] | None = None) -> None:
    if context is None:
        raise ValueError("This subfile is intended to be called by master.py with context.")
    run(context)
