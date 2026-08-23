#!/usr/bin/env python3
"""Check that the 2x beam source is the same on every field line as in 1x.

Run after both simulations have written frame-zero source diagnostics, e.g.

  python check_source_equivalence.py --one-x-dir ../1x-beams --two-x-dir .

The local M0/M2 fit is the relevant source-amplitude check.  Jacobian-weighted
power is also reported, but it measures power per unit psi and can vary between
otherwise identical field lines because their flux-tube volumes differ.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

try:
    import postgkyl as pg

    POSTGKYL_V2 = hasattr(pg.GData, "interpolate")
except (ImportError, AttributeError):
    POSTGKYL_V2 = False

if not POSTGKYL_V2:
    from postgkyl.data import GData, GInterpModal


TRAPEZOID = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
M_DEUTERIUM = 2.014 * 1.67262192369e-27


def interpolated_values(path: Path, component: int = 0) -> np.ndarray:
    if POSTGKYL_V2:
        data = pg.GData(str(path), basis_type="ms", poly_order=1)
        nodal = data.select(comp=component).interpolate()
        return np.squeeze(nodal.get_values())
    else:
        data = GData(str(path))
        interp = GInterpModal(data, 1, "ms")
        interp.interpolate(component, overwrite=True)
        return np.squeeze(data.get_values()[..., 0])


def mapped_coordinates(path: Path) -> np.ndarray:
    if POSTGKYL_V2:
        data = pg.GData(str(path), basis_type="ms", poly_order=1)
        return np.asarray(data.interpolate().get_values())
    else:
        data = GData(str(path))
        ndim = data.get_num_dims()
        _, values = GInterpModal(data, 1, "ms").interpolate(tuple(range(ndim)))
        return np.asarray(values)


def fit_multiplier(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Return a such that a*candidate best matches reference in L2."""
    cutoff = max(np.max(np.abs(reference)), np.max(np.abs(candidate))) * 1.0e-12
    keep = (np.abs(reference) > cutoff) | (np.abs(candidate) > cutoff)
    denominator = np.dot(candidate[keep], candidate[keep])
    return float(np.dot(reference[keep], candidate[keep]) / denominator)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--one-x-dir", type=Path, required=True)
    parser.add_argument("--two-x-dir", type=Path, required=True)
    parser.add_argument("--one-x-prefix", default="zzim")
    parser.add_argument("--two-x-prefix", default="zzim")
    parser.add_argument("--frame", type=int, default=0)
    parser.add_argument(
        "--one-x-radial-coordinate",
        choices=("psi", "sqrt-psi"),
        default="psi",
        help="radial coordinate used to construct the 1x geometry Jacobian",
    )
    parser.add_argument(
        "--two-x-radial-coordinate",
        choices=("psi", "sqrt-psi"),
        default="psi",
        help="radial coordinate used to construct the 2x geometry Jacobian",
    )
    parser.add_argument(
        "--psi-1x",
        type=float,
        default=None,
        help="1x flux label; defaults to the outermost interpolated 2x psi",
    )
    args = parser.parse_args()

    one = args.one_x_dir / args.one_x_prefix
    two = args.two_x_dir / args.two_x_prefix
    frame = args.frame

    coords_1x = mapped_coordinates(Path(f"{one}-mc2nu_pos_deflated.gkyl"))
    coords_2x = mapped_coordinates(Path(f"{two}-mc2nu_pos_deflated.gkyl"))
    z_1x = coords_1x[:, 0]
    psi = coords_2x[:, 0, 0]
    z_2x = coords_2x[0, :, 1]

    m0_1x = interpolated_values(Path(f"{one}-ion_source_M0_{frame}.gkyl"))
    m2_1x = interpolated_values(Path(f"{one}-ion_source_M2_{frame}.gkyl"))
    jac_1x = interpolated_values(Path(f"{one}-jacobgeo.gkyl"))

    m0_2x = interpolated_values(Path(f"{two}-ion_source_M0_{frame}.gkyl"))
    m2_2x = interpolated_values(Path(f"{two}-ion_source_M2_{frame}.gkyl"))
    jac_2x = interpolated_values(Path(f"{two}-jacobgeo.gkyl"))

    # Compare local source moments on a common Z grid.  These multipliers answer
    # the actual question: by what factor would the 2x source need to change to
    # reproduce the 1x local source?
    m0_ref = np.interp(z_2x, z_1x, m0_1x)
    m2_ref = np.interp(z_2x, z_1x, m2_1x)
    scale_m0 = np.array([fit_multiplier(m0_ref, line) for line in m0_2x])
    scale_m2 = np.array([fit_multiplier(m2_ref, line) for line in m2_2x])

    psi_1x = float(psi[-1] if args.psi_1x is None else args.psi_1x)

    # J_c is the Jacobian for the radial coordinate selected by fl_coord.  For
    # q=sqrt(psi), J_q = J_psi*dpsi/dq = 2*q*J_psi.  Therefore the raw line
    # integral is dP/dq and must either be integrated over q or divided by 2*q
    # before it is interpreted as dP/dpsi.
    line_power_2x_raw = (
        np.pi * M_DEUTERIUM * TRAPEZOID(m2_2x * jac_2x, x=z_2x, axis=1)
    )
    if args.two_x_radial_coordinate == "sqrt-psi":
        radial_grid = np.sqrt(psi)
        dp_dpsi = line_power_2x_raw / (2.0 * radial_grid)
    else:
        radial_grid = psi
        dp_dpsi = line_power_2x_raw
    total_power_2x = TRAPEZOID(line_power_2x_raw, x=radial_grid)

    line_power_1x_raw = (
        np.pi * M_DEUTERIUM * TRAPEZOID(m2_1x * jac_1x, x=z_1x)
    )
    if args.one_x_radial_coordinate == "sqrt-psi":
        dp_dpsi_1x = line_power_1x_raw / (2.0 * np.sqrt(psi_1x))
    else:
        dp_dpsi_1x = line_power_1x_raw

    line_power_2x = float(np.interp(psi_1x, psi, dp_dpsi))
    kappa_1x_from_2x = line_power_2x / dp_dpsi_1x

    print("Local multiplier to apply to the 2x source (target: 1):")
    print(
        f"  from M0: median={np.median(scale_m0):.8f}, "
        f"range=[{scale_m0.min():.8f}, {scale_m0.max():.8f}]"
    )
    print(
        f"  from M2: median={np.median(scale_m2):.8f}, "
        f"range=[{scale_m2.min():.8f}, {scale_m2.max():.8f}]"
    )
    print()
    print("Power bookkeeping (not a local source normalization):")
    print(
        "  radial coordinates             = "
        f"1x:{args.one_x_radial_coordinate}, 2x:{args.two_x_radial_coordinate}"
    )
    print(f"  P_2x                         = {total_power_2x:.8e} W")
    print(f"  (dP/dpsi)_1x                = {dp_dpsi_1x:.8e} W/psi")
    print(f"  (dP/dpsi)_2x at psi_1x      = {line_power_2x:.8e} W/psi")
    print(f"  psi_1x used                 = {psi_1x:.8e}")
    print(f"  kappa to apply to 1x source = {kappa_1x_from_2x:.8f}")
    print(f"  inverse (apply to 2x)       = {1.0 / kappa_1x_from_2x:.8f}")


if __name__ == "__main__":
    main()
