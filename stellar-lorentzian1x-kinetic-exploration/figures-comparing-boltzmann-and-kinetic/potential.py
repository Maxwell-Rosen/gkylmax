"""Plot normalized potential overlay for Boltzmann vs kinetic runs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

from field_plot_common import (
    BOLTZMANN_COLOR,
    KINETIC_COLOR,
    Z_CENTER,
    Z_MIRROR_THROAT,
    Z_SHEATH,
    delta_labels,
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

    kin_label, kin_d1, kin_d2 = delta_labels(z_kin, phi_kin_norm, "Kinetic")
    boltz_label, boltz_d1, boltz_d2 = delta_labels(z_boltz, phi_boltz_norm, "Boltzmann")

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.plot(z_kin, phi_kin_norm, color=KINETIC_COLOR, lw=2.0, label="Kinetic electrons")
    ax.plot(z_boltz, phi_boltz_norm, color=BOLTZMANN_COLOR, lw=2.0, label="Boltzmann electrons")

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
    ax.set_ylabel(r"$e\phi / T_e$")
    # ax.set_title("Potential profile at last frame (normalized)")
    # ax.grid(alpha=0.3)
    ax.legend(loc="best")

    annotation = kin_label + "\n\n" + boltz_label
    ax.text(
        0.02,
        0.45,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10.5,
        linespacing=1.25,
        bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.95},
    )
    ax.set_xlim(-Z_SHEATH, Z_SHEATH)

    fig.tight_layout()
    fig.savefig(out_dir / "potential-overlay-kinetic-vs-boltzmann.pdf")
    plt.close(fig)

    print("[potential] Saved figure to", out_dir)
    print(
        "[potential] Te used [eV]: "
        f"kinetic(z=0)={data['kin_te0_ev']:.3f}, boltzmann(master)={data['boltz_te0_ev']:.3f}"
    )
    print("[potential]", kin_label)
    print("[potential]", boltz_label)
    print(
        "[potential] Delta summary "
        f"(center->throat, throat->sheath) [dimensionless]: "
        f"kinetic=({kin_d1:.3f}, {kin_d2:.3f}), "
        f"boltzmann=({boltz_d1:.3f}, {boltz_d2:.3f})"
    )


def main(context: dict[str, Any] | None = None) -> None:
    if context is None:
        raise ValueError("This subfile is intended to be called by master.py with context.")
    run(context)
