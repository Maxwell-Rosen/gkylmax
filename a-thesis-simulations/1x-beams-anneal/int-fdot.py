#!/usr/bin/env python3
"""Compute 1D M0 marginals of a raw 1x2v gyrokinetic fdot file.

The input is the native modal ``J*fdot`` written by
``gk_species_write_fdot``.  Choose whether the plotted x axis is physical
``z``, ``vpar``, or ``mu``; the other two phase-space directions are
integrated out.  A configuration-space map and a velocity-space map convert
the computational grid edges to their non-uniform physical coordinates.

The raw array already contains Gkeyll's configuration, magnetic, and
velocity-map Jacobians.  The maps therefore change the displayed coordinate
and the density per unit x, but their Jacobians must not be multiplied into
the values a second time.

The script also constructs the native, 1D, serendipity-p1 M0 field using the
exact contraction in Gkeyll's ``gyrokinetic_M0_1x2v_ser_p1`` kernel.  Its
native full-grid integral is compared with the selected marginal's integral
and both values are printed.
"""

from __future__ import annotations

import argparse

import numpy as np

import postgkyl as pg
from postgkyl.gpython import GkylArray


AXES = {"z": 0, "vpar": 1, "mu": 2}


def _validate_fdot(fdot: pg.GData) -> None:
  """Check the phase-space layout assumed by the generated Gkeyll kernel."""
  if fdot.backend != "gkyl":
    raise ValueError("fdot must be loaded with Postgkyl's native Gkeyll backend")
  # end
  if (
      fdot.num_dims != 3
      or int(fdot.ctx.get("vdim", -1)) != 2
      or str(fdot.ctx.get("basis_type", "")).lower() != "gkhybrid"
      or int(fdot.ctx.get("poly_order", -1)) != 1
  ):
    raise ValueError(
        "expected 1x2v gkhybrid-p1 data with axes (z, vpar, mu)")
  # end
# end


def fdot_m0_1x2v(fdot_file: str | pg.GData) -> pg.GData:
  """Velocity-integrate a raw 1x2v p1 ``J*fdot`` into a native M0 field."""
  fdot = pg.load(fdot_file) if isinstance(fdot_file, str) else fdot_file
  _validate_fdot(fdot)

  cells = np.asarray(fdot.ctx["cells"], dtype=int)
  lower = np.asarray(fdot.ctx["lower"], dtype=float)
  upper = np.asarray(fdot.ctx["upper"], dtype=float)

  mass = float(fdot.ctx["mass"])
  dz = (upper - lower) / cells

  # Gkeyll's 1x2v p1 M0 kernel maps phase coefficients 0 and 1 to the
  # configuration-space p1 coefficients.  J_geo, B, and the velocity-map
  # Jacobian are already included in the raw J*fdot coefficients.
  phase_coeff = fdot.values[..., [0, 1]]
  m0_coeff = (
      np.pi * dz[1] * dz[2] / mass
      * phase_coeff.sum(axis=(1, 2))
  )

  ctx = dict(fdot.ctx)
  ctx.update({
      "basis_type": "serendipity",
      "poly_order": 1,
      "value_form": "modal",
      "cells": cells[:1],
      "lower": lower[:1],
      "upper": upper[:1],
      "num_comps": 2,
      "vdim": 0,
      "Description": "M0 moment of the raw gyrokinetic stage-1 fdot array",
  })

  m0 = pg.GData(ctx=ctx, tag="fdot_M0", label=r"$J_{geo}\,\dot{M}_0$")
  m0.push(
      [np.asarray(fdot.grid[0], dtype=float)],
      GkylArray.from_numpy(m0_coeff),
  )
  return m0
# end


def _phase_cell_m0(fdot: pg.GData) -> np.ndarray:
  """M0 contribution of every complete (z, vpar, mu) phase-space cell."""
  cells = np.asarray(fdot.ctx["cells"], dtype=int)
  lower = np.asarray(fdot.ctx["lower"], dtype=float)
  upper = np.asarray(fdot.ctx["upper"], dtype=float)
  dz = (upper - lower) / cells
  mass = float(fdot.ctx["mass"])

  # Integral of the normalized constant basis function over a 3D reference
  # cell gives 2**(-3/2) after the physical-cell map.  The gyrokinetic
  # velocity measure supplies 2*pi/m.
  factor = 2.0 * np.pi / mass * np.prod(dz) / 2.0 ** 1.5
  return factor * fdot.values[..., 0]
# end


def _mapped_phase_grid(fdot: pg.GData, velocity_map: str,
    configuration_map: str) -> list[np.ndarray]:
  """Map phase-space grid edges using the same verbs as ``load_gk_distf``.

  The configuration map must be the scalar/deflated 1D map (two p1 modal
  components), for example ``geo_corn_mc2nu_pos_deflated.gkyl``.  A full
  three-component Cartesian ``mapc2p`` curve has no unique scalar x axis.
  """
  # ``map`` only changes a NumPy-backed target's grid.  This lightweight
  # target supplies the phase-grid topology; its values are never used.
  target = pg.GData()
  target.push(
      [np.asarray(g, dtype=float) for g in fdot.grid],
      np.zeros(tuple(int(n) for n in fdot.num_cells) + (1,)),
  )

  target = target.map(velocity_map, space="vel",
      basis_type="serendipity", poly_order=1)
  target = target.map(configuration_map, space="conf",
      basis_type="serendipity", poly_order=1)
  return [np.asarray(g, dtype=float) for g in target.grid]
# end


def marginal_profile(fdot: pg.GData, x_axis: str, mapped_grid: list[np.ndarray],
    *, absolute: bool = False) -> tuple[pg.GData, float]:
  """Integrate out the other two axes and return density per physical x.

  With ``absolute=True``, contributions are multiplied by the sign of the
  velocity-integrated M0 in their z cell before reduction.  Consequently all
  three possible marginals integrate to Gkeyll's cellwise-absolute M0 total.
  The z marginal is non-negative; vpar/mu marginals may contain negative
  velocity-bin contributions even though their total is non-negative.
  """
  axis = AXES[x_axis]
  cell_contrib = _phase_cell_m0(fdot)

  if absolute:
    m0_z_cell = cell_contrib.sum(axis=(1, 2))
    cell_contrib = cell_contrib * np.sign(m0_z_cell)[:, None, None]
  # end

  reduce_axes = tuple(d for d in range(3) if d != axis)
  marginal_bin = cell_contrib.sum(axis=reduce_axes)

  edges = mapped_grid[axis]
  widths = np.diff(edges)
  if edges.ndim != 1 or len(edges) != len(marginal_bin) + 1:
    raise ValueError(
        f"mapped {x_axis} grid has shape {edges.shape}; expected a separable "
        f"1D edge grid of length {len(marginal_bin) + 1}")
  # end
  if np.any(widths <= 0.0):
    raise ValueError(f"mapped {x_axis} coordinate is not strictly increasing")
  # end

  # Store a cell-average density so Postgkyl's integrate_axis multiplies it
  # by the non-uniform physical cell widths and recovers marginal_bin.sum().
  values = (marginal_bin / widths)[:, None]
  mode = "absolute-M0 contribution" if absolute else "signed M0"
  profile = pg.GData(
      ctx={
          "time": fdot.ctx.get("time"),
          "frame": fdot.ctx.get("frame"),
          "Description": f"{mode} marginal of stage-1 gyrokinetic fdot",
      },
      tag=f"fdot_M0_{x_axis}",
      label=rf"$d\dot{{M}}_0/d\,{x_axis}$",
  )
  profile.push([edges], values)

  # Exercise Postgkyl's public non-uniform-axis reduction for the check.
  integrated = float(profile.integrate_axis(0).values.squeeze())
  return profile, integrated
# end


def jdiagnostic_abs_total(m0: pg.GData) -> float:
  """Cellwise-absolute M0 total used by the integrated diagnostic.

  Gkeyll first integrates M0 over each configuration cell and then takes its
  absolute value, preventing cancellation between z cells.
  """
  dz0 = np.diff(np.asarray(m0.grid[0], dtype=float))
  m0_cell = dz0 / np.sqrt(2.0) * m0.values[:, 0]
  return float(np.abs(m0_cell).sum())
# end


def _parse_args() -> argparse.Namespace:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--fdot_file", help="<name>-<species>_fdot_<frame>.gkyl", default = "zzim-ion_fdot_84.gkyl")
  parser.add_argument("--x-axis", choices=tuple(AXES), default="z",
      help="physical coordinate retained in the plotted marginal (default: z)")
  parser.add_argument("--velocity-map",
      help="species mapc2p_vel.gkyl file",
      default = "zzim-ion_mapc2p_vel.gkyl")
  parser.add_argument("--configuration-map",
      help="scalar 1D configuration map, e.g. mc2nu_pos_deflated.gkyl",
      default = "zzim-geo_corn_mc2nu_pos_deflated.gkyl")
  parser.add_argument("--output", help="write the native 1D M0 field to this .gkyl file")
  parser.add_argument("--plot", help="write a PNG/PDF plot to this path")
  parser.add_argument("--absolute", action="store_true",
      help="decompose the built-in-diagnostic-style cellwise-absolute M0 total")
  parser.add_argument("--show", action="store_true", help="show the plot interactively")
  return parser.parse_args()
# end


def main() -> None:
  args = _parse_args()
  fdot = pg.load(args.fdot_file)
  _validate_fdot(fdot)
  m0 = fdot_m0_1x2v(fdot)
  mapped_grid = _mapped_phase_grid(
      fdot, args.velocity_map, args.configuration_map)
  profile, profile_total = marginal_profile(
      fdot, args.x_axis, mapped_grid, absolute=args.absolute)

  # This independent terminal reduction stays in Gkeyll's native DG backend.
  signed_total = float(m0.integrate())
  abs_total = diagnostic_abs_total(m0)
  expected_total = abs_total if args.absolute else signed_total
  error = profile_total - expected_total
  scale = max(abs_total, abs(signed_total), 1.0)

  print(f"x axis:                                {args.x_axis}")
  print(f"signed total (native M0 integral):     {signed_total:.16e}")
  print(f"absolute M0 total (cellwise abs):      {abs_total:.16e}")
  print(f"selected marginal integral:            {profile_total:.16e}")
  print(f"marginal minus selected total:         {error:.16e}")
  print(f"normalized validation error:           {abs(error) / scale:.3e}")

  tolerance = 128.0 * np.finfo(float).eps * scale
  if abs(error) > tolerance:
    raise RuntimeError(
        f"marginal validation failed: |{error}| > tolerance {tolerance}")
  # end

  if args.output:
    print(f"wrote {m0.save(args.output)}")
  # end

  if args.plot or args.show:
    profile.plot(saveas=args.plot, show=args.show)
  # end
# end


if __name__ == "__main__":
  main()
# end
