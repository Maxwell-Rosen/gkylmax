"""Shared moment-profile comparison for the z, vpar, and mu scans."""
from pathlib import Path
import argparse
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import postgkyl as pg


def main(scan_dir, axis):
    scan_dir = Path(scan_dir).resolve()
    parser = argparse.ArgumentParser(description=f"Compare {axis} resolution profiles at one saved frame.")
    parser.add_argument("--frame", type=int, default=65)
    parser.add_argument("--baseline", type=Path, default=scan_dir.parent.parent / "1x-beams")
    parser.add_argument("--no-baseline", action="store_true")
    parser.add_argument("--output", type=Path, default=scan_dir / f"res-scan-{axis}.pdf")
    args = parser.parse_args()
    if args.frame < 0:
        parser.error("--frame must be nonnegative")
    if not args.no_baseline and not args.baseline.is_dir():
        parser.error(f"Baseline directory does not exist: {args.baseline}")

    symbol = {"z": r"N_z", "vpar": r"N_{v_\parallel}", "mu": r"N_\mu"}[axis]
    paths = sorted((p for p in scan_dir.iterdir() if p.is_dir() and p.name.isdigit()),
                   key=lambda p: int(p.name))
    candidates = [(p, f"${symbol}={int(p.name)}$") for p in paths]
    if not args.no_baseline:
        candidates.append((args.baseline, "Original (1x-beams)"))
    datasets, labels, colors = [], [], []
    end_time = None
    palette = plt.get_cmap("tab10")
    for i, (folder, label) in enumerate(candidates):
        moment = folder / f"zzim-ion_BiMaxwellianMoments_{args.frame}.gkyl"
        mapping = folder / "zzim-geo_corn_mc2nu_pos_deflated.gkyl"
        missing = [p.name for p in (moment, mapping) if not p.is_file()]
        if missing:
            print(f"SKIP {folder}: missing {', '.join(missing)}", file=sys.stderr)
            continue
        data = pg.load(str(moment))
        time = float(data.ctx["time"])
        if not math.isfinite(time):
            parser.error(f"{moment}: non-finite frame time")
        if end_time is not None and not math.isclose(time, end_time, rel_tol=0, abs_tol=1e-11):
            parser.error(f"{moment}: frame time {time} differs from {end_time}")
        end_time = time
        data = data.interpolate().map(str(mapping))
        # Convert the saved ion temperature (velocity squared) to keV.
        data[..., 2:4] *= 1.602176634e-27 * 2.014 / 1.602176634e-19 / 1000
        datasets.append(data)
        labels.append(label)
        colors.append("black" if label.startswith("Original") else palette(i % 10))
        print(f"Loaded {folder} (t={time:g} s)")

    if len(datasets) < 2:
        parser.exit(1, "Need at least two runs with moment and mapping output for comparison.\n")
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    direction = {"z": r"$z$", "vpar": r"$v_\parallel$", "mu": r"$\mu$"}[axis]
    pg.plot(
        *datasets,
        color=colors,
        figure=0,
        legend_labels=labels,
        legend_subplot=0,
        legend_loc="best",
        title=f"Resolution Scan along {direction} (t={end_time:g} s)",
        xlabel="z [m]",
        figsize="12,10",
        split_linear_log=True,
        split_point=0.0,
        split_log_side="right",
        split_width_ratios=(1.0, 1.0),
        split_gap=0.0,
        split_legend_side="log",
        split_log_nonpositive="mask",
        split_linear_ylim=[
            (0.0, 3e19),
            (-1.4e6, 1.4e6),
            (-2, 14.0),
            (0.0, 27.0),
        ],
        split_log_ylim=[
            (1e13, 4e19),
            (1e2, 2e6),
            (1e-1, 20.0),
            (1e-1, 50.0),
        ],
        subplot_ylabels=(
            r"Density [$\mathrm{m}^{-3}$],"
            r"$U_\parallel$ [m/s],"
            r"$T_\parallel$ [keV],"
            r"$T_\perp$ [keV]"
        ),
        saveas=str(output),
        no_show=True,
    )
    print(f"Wrote {output}")
