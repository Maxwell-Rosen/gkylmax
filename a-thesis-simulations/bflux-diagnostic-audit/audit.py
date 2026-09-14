"""Reconstruct 1x2v P1 wall and MPI-face fluxes from saved nodal surface flux.

Run with the user's pgkyl environment. No simulation files are modified.
The eight spatial ranks and 400 cells match these three saved runs.
"""
from pathlib import Path
import csv
import numpy as np
import postgkyl as pg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
RUNS = ["1x-beams", "1x-maxw", "1x-maxw/no-crazy-time-dilation"]
# Nodal-to-modal coefficient Ghat[0] in dg_gyrokinetic_boundary_surfx_1x2v_ser_p1.
WEIGHTS = np.array([5/18, 4/9, 5/18] * 2)
rows = []
fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

for ax, run in zip(axes, RUNS):
    prefix = ROOT / run
    source = pg.GData(str(prefix / "zzim-ion_source_integrated_moms.gkyl"))
    diag = pg.GData(str(prefix / "zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl"))
    times = diag.get_grid()[0]
    s_times = source.get_grid()[0]
    s_values = source.get_values()[:, 0]
    diag_ratio = 2 * diag.get_values()[:, 0] / np.interp(times, s_times, s_values)
    current = []
    for frame in range(10, 45):
        path = prefix / f"zzim-ion_collisionless_surf_flux_{frame}.gkyl"
        if not path.exists():
            continue
        data = pg.GData(str(path))
        flux = data.get_values()
        assert flux.shape == (400, 64, 32, 10)
        grid = data.get_grid()
        dv = grid[1][1] - grid[1][0]
        dmu = grid[2][1] - grid[2][0]
        # Combining the boundary DG kernel, M0 kernel and configuration integral:
        # outward lower-face Gamma = -pi*dv*dmu/m * sum_velocity Ghat[0].
        # The saved flux already includes mapped-velocity Jacobians.
        gamma = -np.pi * dv * dmu / data.ctx["mass"] * np.sum(
            flux[..., :6] * WEIGHTS, axis=(1, 2, 3))
        t = data.ctx["time"]
        i = int(np.argmin(abs(times - t)))
        s = float(np.interp(t, s_times, s_values))
        rank_ratios = 2 * gamma[::50] / s
        row = dict(run=run, frame=frame, time=t, sample=i,
                   wall_ratio=rank_ratios[0],
                   internal_ratio=rank_ratios[1:].sum(),
                   reconstructed_diag_ratio=rank_ratios.sum(),
                   saved_diag_ratio=diag_ratio[i],
                   snapshot_minus_rk_diag=rank_ratios.sum()-diag_ratio[i])
        row.update({f"rank_{r}_lower_ratio": value for r, value in enumerate(rank_ratios)})
        rows.append(row)
        current.append(row)
    ax.plot(np.arange(len(times)), diag_ratio, color="#bf5939", lw=1.5,
            label="Saved integrated diagnostic")
    ax.plot([r["sample"] for r in current], [r["wall_ratio"] for r in current],
            "o-", color="#006f9f", ms=3, lw=1.5,
            label="Actual lower wall, doubled (snapshots)")
    tau_fdp = .005 if run == "1x-maxw" else .000015
    for cycle in range(1, 5):
        lo = cycle * (.1 + tau_fdp)
        start = np.searchsorted(times, lo)
        stop = np.searchsorted(times, lo + .1)
        ax.axvspan(start, stop, color="gray", alpha=.12)
    ax.axhline(1, color="black", lw=.7, ls=":")
    ax.set_title(run, loc="left")
    ax.set_ylabel(r"$2\Gamma_\mathrm{lower}/S$")
    ax.grid(alpha=.2)
axes[0].set_ylim(-1, 21)
axes[1].set_ylim(-1, 16)
axes[2].set_ylim(-1, 36)
axes[0].legend(loc="upper right", fontsize=9)
axes[-1].set_xlim(1000, 4505)
axes[-1].set_xlabel("Integrated diagnostic sample index; gray bands: orbit-average phases")
fig.tight_layout()
fig.savefig(OUT / "wall_vs_diagnostic.png", dpi=160)
fig.savefig(OUT / "wall_vs_diagnostic.pdf")
with (OUT / "flux_comparison.csv").open("w") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

for row in rows:
    if row["frame"] == 44:
        print(row)
# The saved flux is an instantaneous RK-stage snapshot, while the dynvec is
# an RK-weighted flux. They agree best away from switches and fast transients.
stable = [r for r in rows if r["frame"] % 10 in (1, 2, 3, 4)]
print("Max reconstruction error inside OAPs:",
      max(abs(r["snapshot_minus_rk_diag"]) for r in stable))
