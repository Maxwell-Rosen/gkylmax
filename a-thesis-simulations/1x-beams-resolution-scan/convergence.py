#!/usr/bin/env python3
"""Resolution convergence for this 1x2v, p=1, Boltzmann-electron scan.

See z-scan/CONVERGENCE.md for diagnostic definitions and normalization.
Only moment/geometry files are read; distribution functions are not loaded.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import re
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import postgkyl as pg

ROOT = Path(__file__).resolve().parent / 'z-scan'
AXES = {'z': ('Nz', r'N_z'), 'vpar': ('Nvpar', r'N_{v_\parallel}'), 'mu': ('Nmu', r'N_\mu')}
COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#332288', '#882255']
# Integrated quantities retain Gkeyll's per-unit transverse-coordinate measure.
INVENTORY = [
    ('Ni', 'Integrated ion number', 'Number [simulation normalization]', 1),
    ('Ntot', 'Integrated ion + electron number', 'Number [simulation normalization]', 1),
    ('Ei', 'Integrated ion kinetic energy', 'Energy [simulation normalization]', 1),
    ('Ee', 'Integrated electron kinetic energy', 'Energy [simulation normalization]', 1),
    ('Etot', 'Total particle kinetic energy', 'Energy [simulation normalization]', 1),
    ('n0', 'Central ion density', r'$n_i(z=0)$ [$\mathrm{m}^{-3}$]', 1),
]
WALL = [
    (f'{key}_{side}', f'{title}: {label} wall', unit, scale)
    for key, title, unit, scale in [
        ('gamma', 'Outward ion particle flux', r'$\Gamma_i$ [$\mathrm{m}^{-2}\,\mathrm{s}^{-1}$]', 1),
        ('qi', 'Ion wall energy flux', r'$q_{i,\mathrm{wall}}$ [MW m$^{-2}$]', 1e-6),
        ('qtot', 'Total wall heat flux (model)', r'$q_{i+e,\mathrm{wall}}$ [MW m$^{-2}$]', 1e-6),
    ] for side, label in [('lo', 'lower'), ('hi', 'upper')]]
SUMMARY = [INVENTORY[0], INVENTORY[4], INVENTORY[5],
           ('loss', 'Particle loss, both walls', 'Rate [simulation normalization]', 1),
           ('Pi', 'Ion wall power, both walls', 'Power [simulation normalization]', 1),
           ('Ptot', 'Total wall power, both walls (model)', 'Power [simulation normalization]', 1)]


def read(path):
    if not path.is_file():
        raise FileNotFoundError(2, 'Missing diagnostic', str(path))
    d = pg.GData(str(path))
    if not np.isfinite(d.get_values()).all():
        raise ValueError(f'Non-finite values in {path}')
    return d


def series(path, components):
    d = read(path)
    t, v = np.asarray(d.get_grid()[0]), d.get_values()
    if v.ndim != 2 or v.shape != (len(t), components):
        raise ValueError(f'Unexpected diagnostic shape in {path}: {v.shape}')
    if np.any(np.diff(t) < 0):
        raise ValueError(f'Time goes backwards in {path}; check appended restart output')
    # Keep the last copy of an exactly duplicated restart/phase-boundary time.
    keep = np.r_[np.diff(t) != 0, True]
    return t[keep], v[keep], d.ctx


def sample(s, times):
    t, v = s
    times = np.asarray(times)
    if np.min(times) < t[0] - 1e-12 or np.max(times) > t[-1] + 1e-12:
        raise ValueError('Attempted to extrapolate a diagnostic in time')
    return np.interp(times, t, v)


def modal(path):
    d = read(path)
    if d.ctx.get('poly_order') != 1 or d.ctx.get('basis_type') != 'serendipity':
        raise ValueError(f'This script requires the scan\'s p=1 serendipity basis: {path}')
    v = d.get_values()
    if v.ndim != 2 or v.shape[1] % 2:
        raise ValueError(f'Expected 1D p=1 modal data: {path}')
    return d, v.reshape(v.shape[0], -1, 2)


def evaluate(coeff, xi):
    """Orthonormal p=1 basis on [-1,1], last coefficient dimension."""
    return coeff[..., 0, None] / np.sqrt(2) + coeff[..., 1, None] * np.sqrt(1.5) * np.asarray(xi)


def product_integral(a, b, grid):
    # Exact integral of a product of two p=1 modal expansions.
    return np.sum(np.diff(grid) * np.sum(a * b, axis=-1) / 2)


def input_signature(path, allow_mapping=False, axis='z'):
    # Only the scanned resolution may differ. Refuse accidental mixed physics.
    s = path.read_text()
    s = re.sub(r'//[^\n]*|/\*.*?\*/', '', s, flags=re.S)
    key = AXES[axis][0]
    s = re.sub(rf'\bint\s+{key}\s*=\s*\d+\s*;', f'int {key} = RESOLUTION;', s)
    if allow_mapping:
        # Only these three reviewed grid parameters may differ for the extra run.
        for key in ['maximum_slope_at_min_B', 'maximum_slope_at_max_B', 'gaussian_std']:
            s = re.sub(rf'(\.{key}\s*=\s*)[0-9.eE+\-]+', rf'\g<1>MAPPING', s)
    return re.sub(r'\s+', '', s)


def read_resolution(path, axis, nz):
    """Use saved cell counts, falling back to input for velocity resolutions."""
    if axis == 'z':
        return nz
    key = AXES[axis][0]
    stats = path / 'zzim-stat.json'
    if stats.is_file():
        # Gkeyll's stat file is not strict JSON (unquoted keys, trailing commas).
        matches = re.findall(r'ion_cells\s*:\s*\[([^]]+)\]', stats.read_text())
        if matches:
            cells = [int(v) for v in re.findall(r'\d+', matches[-1])]
            if len(cells) == 3 and cells[0] == nz:
                return cells[{'vpar': 1, 'mu': 2}[axis]]
    source = re.sub(r'//[^\n]*|/\*.*?\*/', '', (path / 'sim.c').read_text(), flags=re.S)
    match = re.search(rf'\bint\s+{key}\s*=\s*(\d+)\s*;', source)
    if match is None or int(match.group(1)) < 1:
        raise ValueError(f'{path}/sim.c: cannot determine positive {key}')
    warnings.warn(f'{path}: using {key} from sim.c; saved velocity cell counts unavailable')
    return int(match.group(1))


def run_label(run):
    return f"${run['symbol']}={run['resolution']}$" + (' (new mapping)' if run['mapping'] == 'nunif' else '')


def line_style(run):
    return dict(lw=2.0, ls='--') if run['mapping'] == 'nunif' else dict(lw=1.2)


def load_run(path, frame, te_ev, final_time=None):
    prefix = path / 'zzim'
    def f(s):
        return Path(f'{prefix}-{s}.gkyl')
    final, _ = modal(f(f'ion_BiMaxwellianMoments_{frame}'))
    end = float(final.ctx['time'])
    if final_time is not None and not np.isclose(end, final_time, rtol=0, atol=1e-11):
        raise ValueError(f'{path}: frame {frame} has a different time ({end})')
    n = int(final.ctx['cells'][0])
    ti, ion, ctx = series(f('ion_integrated_moms'), 3)
    te, elc, ectx = series(f('elc_integrated_moms'), 4)
    charge, mass, emass = ctx['charge'], ctx['mass'], ectx['mass']
    temp = te_ev * abs(ectx['charge'])
    out = {'Ni': (ti, ion[:, 0]), 'Hi': (ti, ion[:, 2]),
           'Ne': (te, elc[:, 0]), 'Ee': (te, .5 * emass * (elc[:, 2] + elc[:, 3]))}
    out['Ntot'] = (ti, ion[:, 0] + sample(out['Ne'], ti))
    # Surface files pack lower-face values, then an upper-boundary copy.
    # They are NOT volume-modal data despite the generic file metadata.
    area_data = read(f('geo_surf0_lenr')).get_values()
    if area_data.shape != (n, 2):
        raise ValueError(f'{path}: unexpected surface-area-factor layout')
    areas = [float(area_data[0, 0]), float(area_data[-1, 1])]
    if min(areas) <= 0:
        raise ValueError(f'{path}: nonpositive wall area factor')
    for side, edge, area in zip(['lo', 'hi'], ['lower', 'upper'], areas):
        t, b, _ = series(f(f'ion_bflux_x{edge}_integrated_HamiltonianMoments'), 3)
        # Both particle/energy fluxes are already outward positive at both ends.
        out[f'loss_{side}'] = (t, b[:, 0])
        out[f'Pi_{side}'] = (t, b[:, 2])
        out[f'Ptot_{side}'] = (t, b[:, 2] + 2 * temp * b[:, 0])
        out[f'gamma_{side}'] = (t, b[:, 0] / area)
        out[f'qi_{side}'] = (t, b[:, 2] / area)
        out[f'qe_{side}'] = (t, 2 * temp * b[:, 0] / area)
        out[f'qtot_{side}'] = (t, (b[:, 2] + 2 * temp * b[:, 0]) / area)
    for key in ['loss', 'Pi', 'Ptot']:
        t, v = out[f'{key}_lo']
        out[key] = (t, v + sample(out[f'{key}_hi'], t))

    geom, jac = modal(f('geo_int_jacobgeo'))
    _, mapping = modal(f('geo_corn_mc2nu_pos_deflated'))
    grid = geom.get_grid()[0]
    profile_times, central, energy = [], [], []
    number_check, hamiltonian_check = [], []
    profile = None
    for k in range(frame + 1):
        d, bimax = modal(f(f'ion_BiMaxwellianMoments_{k}'))
        _, m0 = modal(f(f'ion_M0_{k}'))
        _, m2 = modal(f(f'ion_M2_{k}'))
        _, phi = modal(f(f'field_{k}'))
        t = float(d.ctx['time'])
        # Physical z=0 can lie on a DG interface: average both one-sided traces.
        xi = -mapping[:, 0, 0] / (np.sqrt(3) * mapping[:, 0, 1])
        cells = np.flatnonzero(abs(xi) <= 1 + 1e-10)
        if not len(cells):
            raise ValueError(f'{path}: physical z=0 is outside the position map')
        n0 = np.mean([evaluate(bimax[j, 0], [np.clip(xi[j], -1, 1)])[0] for j in cells])
        ei = .5 * mass * product_integral(m2[:, 0], jac[:, 0], grid)
        ni = product_integral(m0[:, 0], jac[:, 0], grid)
        # Independent reconstruction checks the Hamiltonian component convention.
        xq, wq = np.polynomial.legendre.leggauss(3)
        potential = charge * np.sum(np.diff(grid)[:, None] / 2 * wq *
                    evaluate(m0[:, 0], xq) * evaluate(jac[:, 0], xq) * evaluate(phi[:, 0], xq))
        number_check.append((ni - sample(out['Ni'], t)) / max(abs(sample(out['Ni'], t)), 1))
        hamiltonian_check.append((ei + potential - sample(out['Hi'], t)) /
                                max(abs(sample(out['Hi'], t)), 1))
        profile_times.append(t)
        central.append(n0)
        energy.append(ei)
        if k == frame:
            xp = np.array([-.5, .5])
            z = evaluate(mapping[:, 0], xp).ravel()
            vals = evaluate(bimax, xp).transpose(0, 2, 1).reshape(-1, 4)
            vals[:, 2:] *= mass / charge / 1000
            profile = (z, vals)
    ts = np.array(profile_times)
    if np.any(np.diff(ts) <= 0):
        raise ValueError(f'{path}: frame times are not strictly increasing')
    out['n0'] = (ts, np.array(central))
    out['Ei'] = (ts, np.array(energy))
    out['Etot'] = (ts, np.array(energy) + sample(out['Ee'], ts))
    # No extrapolation: every curve must cover the selected comparison frame.
    for s in out.values():
        sample(s, end)
    checks = {'max_relative_number_reconstruction': float(np.max(abs(np.array(number_check)))),
              'max_relative_H_reconstruction': float(np.max(abs(np.array(hamiltonian_check))))}
    checks['negative_final_Tpar_samples'] = int(np.sum(profile[1][:, 2] < 0))
    return dict(path=path, nz=n, end=end, data=out, profile=profile, areas=areas,
                checks=checks, changeset=ctx.get('changeset', 'unknown'),
                physical_edges=evaluate(mapping[:, 0], [-1, 1]))


def window_stats(s, start, end):
    t, _ = s
    times = np.r_[start, t[(t > start) & (t < end)], end]
    values = sample(s, times)
    mean = np.trapezoid(values, times) / (end - start)
    return mean, np.min(values), np.max(values)


def save(fig, output, name):
    fig.savefig(output / f'{name}.pdf', bbox_inches='tight')
    fig.savefig(output / f'{name}.png', dpi=170, bbox_inches='tight')
    plt.close(fig)


def histories(runs, metrics, output, name, end, zoom_start=None):
    fig, axes = plt.subplots(3, 2, figsize=(12, 10), layout='constrained')
    for ax, (key, title, unit, scale) in zip(axes.flat, metrics):
        for run, color in zip(runs, (COLORS[i % len(COLORS)] for i in range(len(runs)))):
            if run['mapping'] == 'nunif':
                color = '#6A3D9A'
            t, v = run['data'][key]
            lower = max(t[0], zoom_start if zoom_start is not None else t[0])
            # Explicit endpoint interpolation, without extrapolating sparse frames.
            tt = np.r_[lower, t[(t > lower) & (t < end)], end]
            xx = (tt - zoom_start) * 1e6 if zoom_start is not None else tt
            ax.plot(xx, sample((t, v), tt) * scale, color=color,
                    label=run_label(run), **line_style(run))
        ax.set(title=title, ylabel=unit,
               xlabel=r'Time in final relaxation [$\mu$s]' if zoom_start is not None else 'Simulation time [s]')
        ax.grid(alpha=.25)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 4), useMathText=True)
    axes.flat[0].legend(fontsize=9)
    fig.suptitle('Resolution scan: ' + ('final relaxation' if zoom_start is not None else 'full history (OAP/FDP simulation clock)'), fontsize=14)
    save(fig, output, name)


def main(root=ROOT, axis='z'):
    root = Path(root).resolve()
    resolution_key, symbol = AXES[axis]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=root)
    parser.add_argument('--baseline', type=Path, default=root.parent.parent / '1x-beams',
                        help='Baseline used in compare.py; default: ../../1x-beams')
    parser.add_argument('--no-baseline', action='store_true')
    parser.add_argument('--exclude', '--exclude-resolutions', type=int, nargs='+',
                        action='extend', default=[], metavar='N',
                        help=f'Exclude these {resolution_key} values, including matching baseline/extra runs; may be repeated')
    parser.add_argument('--output', type=Path,
                        help='Default: convergence-plots, or convergence-plots-with-nunif with --nonuniform-run')
    parser.add_argument('--nonuniform-run', type=Path,
                        help='Additional mapped-grid run; retained as a separate point at its actual Nz')
    parser.add_argument('--frame', type=int, default=65, help='Common comparison frame (default: completed run, 65)')
    parser.add_argument('--window-us', type=float, default=15, help='Trailing averaging window, microseconds')
    parser.add_argument('--zoom-us', type=float, default=60, help='Final relaxation zoom, microseconds')
    parser.add_argument('--electron-temperature-ev', type=float, default=940)
    args = parser.parse_args()
    if args.frame < 1 or min(args.window_us, args.zoom_us, args.electron_temperature_ev) <= 0:
        parser.error('Frame, time windows and electron temperature must be positive')
    if any(n < 1 for n in args.exclude):
        parser.error('Excluded resolutions must be positive integers')
    if not args.no_baseline and not args.baseline.is_dir():
        parser.error(f'Baseline directory does not exist: {args.baseline}')
    if args.nonuniform_run and axis != 'z':
        parser.error('--nonuniform-run is only supported for the z scan')
    candidates = sorted((p for p in args.root.iterdir() if p.is_dir() and p.name.isdigit()), key=lambda p: int(p.name))
    if not args.no_baseline and args.baseline.is_dir():
        candidates.append(args.baseline)
    runs, skipped = [], []
    excluded = set(args.exclude)

    def is_excluded(path):
        if not excluded:
            return False
        last = path / f'zzim-ion_BiMaxwellianMoments_{args.frame}.gkyl'
        if not last.is_file():
            return False
        moment, _ = modal(last)
        resolution = read_resolution(path, axis, int(moment.ctx['cells'][0]))
        if resolution in excluded:
            skipped.append(f'{path}: excluded {resolution_key}={resolution}')
            return True
        return False

    if args.nonuniform_run is not None and is_excluded(args.nonuniform_run):
        args.nonuniform_run = None
    signature = None
    for path in candidates:
        if is_excluded(path):
            continue
        last = path / f'zzim-ion_BiMaxwellianMoments_{args.frame}.gkyl'
        if not last.exists():
            skipped.append(f'{path}: frame {args.frame} is not available')
            continue
        current = input_signature(path / 'sim.c', axis=axis)
        if signature is not None and current != signature:
            raise ValueError(f'{path}/sim.c differs in more than {resolution_key}; check physics compatibility')
        print(f'Reading {path}', flush=True)
        try:
            run = load_run(path, args.frame, args.electron_temperature_ev, runs[0]['end'] if runs else None)
        except FileNotFoundError as exc:
            skipped.append(f'{path}: missing required diagnostic: {exc.filename}')
            continue
        signature = current
        run['mapping'] = 'original'
        run['resolution'] = read_resolution(path, axis, run['nz'])
        run['symbol'] = symbol
        run['baseline'] = not args.no_baseline and path.resolve() == args.baseline.resolve()
        runs.append(run)
    for message in skipped:
        print('SKIP:', message)
    if len(runs) < 2:
        parser.error('At least two resolutions with the requested frame and diagnostics are required')
    runs.sort(key=lambda r: r['resolution'])
    reference = runs[-1]  # The extra mapping must never replace the scan reference.
    if args.nonuniform_run is not None:
        path = args.nonuniform_run.resolve()
        if path in {r['path'].resolve() for r in runs}:
            parser.error('The extra run is already in the original scan')
        if input_signature(path / 'sim.c', True) != input_signature(reference['path'] / 'sim.c', True):
            raise ValueError(f'{path}/sim.c differs in more than Nz and the three allowed mapping parameters')
        print(f'Reading new mapping: {path}', flush=True)
        run = load_run(path, args.frame, args.electron_temperature_ev, reference['end'])
        run['mapping'] = 'nunif'
        run['resolution'] = run['nz']
        run['symbol'] = symbol
        runs.append(run)
    if axis != 'z' and len({r['nz'] for r in runs}) > 1:
        warnings.warn('Spatial cell counts differ across these velocity-scan runs; '
                      'the baseline comparison also includes a spatial-resolution change. See run_report.txt.')
    end = runs[0]['end']
    start = end - args.window_us * 1e-6
    zoom_start = end - args.zoom_us * 1e-6
    if min(start, zoom_start) < 0:
        parser.error('Requested time window extends before the simulation start')
    output = args.output or args.root / ('convergence-plots-with-nunif' if args.nonuniform_run else 'convergence-plots')
    if args.nonuniform_run and output.resolve() == (args.root / 'convergence-plots').resolve():
        parser.error('Use a separate output folder for the extra mapping')
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    for metrics, name in [(INVENTORY, 'inventories'), (WALL, 'wall_fluxes')]:
        histories(runs, metrics, output, name + '_history', end)
        histories(runs, metrics, output, name + '_final_relaxation', end, zoom_start)

    # Every individual metric, not just the six overview quantities, is exported.
    rows = []
    for run in runs:
        for key, s in run['data'].items():
            mean, low, high = window_stats(s, start, end)
            ref_end = float(sample(reference['data'][key], end))
            ref_mean = window_stats(reference['data'][key], start, end)[0]
            rows.append(dict({resolution_key: run['resolution']}, run=str(run['path']), mapping=run['mapping'], metric=key, final=float(sample(s, end)), mean=mean,
                             window_min=low, window_max=high,
                             final_difference_percent=100 * (sample(s, end) - ref_end) / abs(ref_end) if ref_end else np.nan,
                             mean_difference_percent=100 * (mean - ref_mean) / abs(ref_mean) if ref_mean else np.nan))
        columns = list(run['data'])
        # Frame-aligned values; original high-cadence diagnostics remain untouched.
        times = run['data']['n0'][0]
        suffix = '_nunif' if run['mapping'] == 'nunif' else ('_baseline' if run.get('baseline') else '')
        with (output / f'frame_diagnostics_{resolution_key}{run["resolution"]}{suffix}.csv').open('w') as fp:
            writer = csv.writer(fp)
            writer.writerow(['time_s'] + columns)
            writer.writerows(zip(times, *(sample(run['data'][k], times) for k in columns)))
    with (output / 'convergence_summary.csv').open('w') as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout='constrained')
    errfig, erraxes = plt.subplots(2, 3, figsize=(14, 8), layout='constrained')
    for ax, ex, (key, title, unit, scale) in zip(axes.flat, erraxes.flat, SUMMARY):
        selected = [r for r in rows if r['metric'] == key and r['mapping'] == 'original']
        nz = np.array([r[resolution_key] for r in selected])
        means = np.array([r['mean'] for r in selected])
        low = np.array([r['window_min'] for r in selected])
        high = np.array([r['window_max'] for r in selected])
        ax.errorbar(nz, means * scale, yerr=np.array([means-low, high-means]) * scale,
                    fmt='o-', capsize=4, label=f'Last {args.window_us:g} µs mean and range')
        ax.plot(nz, [r['final'] * scale for r in selected], 's--', label='Final value')
        ex.plot(nz[:-1], [r['final_difference_percent'] for r in selected[:-1]], 's--', label='Final value')
        ex.plot(nz[:-1], [r['mean_difference_percent'] for r in selected[:-1]], 'o-', label='Window mean')
        extra = [r for r in rows if r['metric'] == key and r['mapping'] == 'nunif']
        for row in extra:
            ax.errorbar(row[resolution_key], row['mean'] * scale,
                        yerr=np.array([[row['mean']-row['window_min']], [row['window_max']-row['mean']]]) * scale,
                        fmt='D', color='#6A3D9A', markersize=8, capsize=5,
                        label='New mapping: mean and range', zorder=5)
            ax.plot(row[resolution_key], row['final'] * scale, '*', color='#6A3D9A',
                    ms=12, label='New mapping: final', zorder=6)
            ex.plot(row[resolution_key], row['mean_difference_percent'], 'D', color='#6A3D9A',
                    ms=8, label='New mapping: mean', zorder=5)
            ex.plot(row[resolution_key], row['final_difference_percent'], '*', color='#6A3D9A',
                    ms=12, label='New mapping: final', zorder=6)
        ex.axhline(0, color='0.5', lw=.7)
        for a in [ax, ex]:
            a.set(title=title, xlabel=f'${symbol}$')
            a.set_xticks(sorted(set((nz if a is ax else nz[:-1]).tolist() + [r[resolution_key] for r in extra])))
            a.grid(alpha=.25)
        ax.set_ylabel(unit)
        ex.set_ylabel('Signed difference [% of reference]')
    axes.flat[0].legend(fontsize=8)
    erraxes.flat[0].legend(fontsize=8)
    fig.suptitle(f'Convergence at t = {end:.9g} s; bars show temporal range, not uncertainty')
    errfig.suptitle(f'Differences relative to original-mapping ${symbol}={reference["resolution"]}$ (finite-resolution reference)')
    save(fig, output, 'convergence_vs_resolution')
    save(errfig, output, 'relative_differences')

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    titles = ['Ion density', 'Parallel velocity', 'Parallel temperature', 'Perpendicular temperature']
    units = [r'$n_i$ [m$^{-3}$]', r'$u_\parallel$ [m s$^{-1}$]', r'$T_\parallel$ [keV]', r'$T_\perp$ [keV]']
    for run, color in zip(runs, (COLORS[i % len(COLORS)] for i in range(len(runs)))):
        if run['mapping'] == 'nunif':
            color = '#6A3D9A'
        z, v = run['profile']
        for i, ax in enumerate(axes.flat):
            ax.plot(z, v[:, i], color=color, label=run_label(run), **line_style(run))
            ax.set(title=titles[i], xlabel='Physical z [m]', ylabel=units[i])
            ax.grid(alpha=.25)
    axes.flat[0].set_yscale('log')
    axes.flat[0].legend()
    fig.suptitle(f'Ion profiles at common frame {args.frame}, t = {end:.9g} s')
    save(fig, output, 'final_profiles')

    if args.nonuniform_run:
        fig, ax = plt.subplots(figsize=(10, 5), layout='constrained')
        for run, color in zip(runs, (COLORS[i % len(COLORS)] for i in range(len(runs)))):
            edges = run['physical_edges']
            widths = edges[:, 1] - edges[:, 0]
            if np.any(widths <= 0):
                raise ValueError(f'{run["path"]}: nonpositive physical cell widths')
            ax.plot(edges.mean(axis=1), widths,
                    color='#6A3D9A' if run['mapping'] == 'nunif' else color,
                    label=run_label(run), **line_style(run))
        ax.set(xlabel='Physical z [m]', ylabel='Physical cell width [m]',
               title='Grid spacing: original scan and new mapping', yscale='log')
        ax.grid(alpha=.25)
        ax.legend()
        save(fig, output, 'grid_spacing')

    lines = [f'Comparison frame: {args.frame}; time: {end:.12g} s',
             f'Window: [{start:.12g}, {end:.12g}] s; reference {resolution_key}={reference["resolution"]}',
             'Integrated numbers, energies and powers retain the simulation transverse normalization.',
             'Wall flux densities use the surface area factor J |grad(z_comp)|.',
             f'Model wall heat: ion Hamiltonian flux + 2 Te Gamma_i, Te={args.electron_temperature_ev:g} eV.',
             'Assumptions: grounded wall, ambipolar losses, Maxwellian Boltzmann electrons.',
             'Total stored energy = ion + electron kinetic energy, including bulk flow.',
             'Field-energy diagnostic excluded: this Boltzmann-field diagnostic is a phi-squared norm.',
             'Range bars measure time variation; no continuum error or convergence order is inferred.',
             'Sparse frame quantities are interpolated linearly in time for window statistics.', '']
    for run in runs:
        spatial = f', Nz={run["nz"]}' if axis != 'z' else ''
        lines.append(f'{resolution_key}={run["resolution"]}{spatial}, mapping={run["mapping"]}: {run["path"]}; code={run["changeset"]}; wall factors={run["areas"]}; checks={run["checks"]}')
        if max(v for k, v in run['checks'].items() if k.startswith('max_relative')) > .01:
            warnings.warn(f'{resolution_key}={run["resolution"]}: reconstruction differs by >1%; inspect report before interpreting energy.')
    lines.extend(['', *['SKIP: ' + s for s in skipped]])
    (output / 'run_report.txt').write_text('\n'.join(lines) + '\n')
    print(f'Wrote {8 if args.nonuniform_run else 7} figures (PDF + PNG), CSV data and run_report.txt to {output}')


if __name__ == '__main__':
    main()
