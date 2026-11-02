#!/usr/bin/env python3
"""Five publication figures combining all three convergence scans.

Uses the diagnostic definitions and raw-data checks in ../diagnostics.py. Default
execution reads current Gkeyll output; --reuse-data explicitly replots the saved
snapshot. No pre-existing convergence-plots CSV is used.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/publication-convergence-mpl')
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import diagnostics as cv
from publication_style import ABS_ROWS, style, finish, paper_axis, history_values
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1] / '1x-beams-resolution-scan'
AXES = ('z', 'vpar', 'mu')
COLORS = ('#0072B2', '#D55E00', '#009E73')
RUN_COLORS = ('#0072B2', '#D55E00', '#009E73', '#CC79A7', '#6A51A3', '#E69F00')
# Restart the history/profile palette in each column, in increasing resolution.
# The production baseline is fourth (purple) in all three scans.
SCAN_COLORS = ('#0072B2', '#D55E00', '#009E73', '#6A3D9A', '#8C564B')
MARKERS = ('o', 's', '^')
OVERVIEW_STYLES = tuple(
    dict(color=color, marker=marker, ls=linestyle, ms=5, lw=1.3)
    for color, marker, linestyle in zip(
        RUN_COLORS, ('o', 's', '^', 'D', 'v', 'P'),
        ('-', '--', '-.', (0, (7, 2)), (0, (5, 1, 1, 1, 1, 1)), ':'))
)


def write_csv(path, rows):
    with path.open('w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def cell_counts(path, nz):
    return [nz, cv.read_resolution(path, 'vpar', nz), cv.read_resolution(path, 'mu', nz)]


def collect(args):
    """Load each physical run once, with scan memberships recorded separately."""
    runs, scans, missing = {}, {}, []
    end = None
    for axis in AXES:
        folder = args.root / f'{axis}-scan'
        paths = sorted((p for p in folder.iterdir() if p.is_dir() and p.name.isdigit()),
                       key=lambda p: int(p.name)) + [args.baseline]
        scans[axis] = []
        signature = cv.input_signature(args.baseline / 'sim.c', axis=axis)
        for path in paths:
            path = path.resolve()
            if cv.input_signature(path / 'sim.c', axis=axis) != signature:
                raise ValueError(f'{path}: physics differs beyond the scanned resolution')
            if not (path / f'zzim-ion_BiMaxwellianMoments_{args.frame}.gkyl').is_file():
                missing.append(dict(axis=axis, resolution=int(path.name), path=str(path),
                                    reason=f'No diagnostics at comparison frame {args.frame}'))
                print(f'UNAVAILABLE: {path}', flush=True)
                continue
            uid = 'baseline' if path == args.baseline else f'{axis}_{path.name}'
            if uid not in runs:
                print(f'Reading {path}', flush=True)
                run = cv.load_run(path, args.frame, args.electron_temperature_ev, end)
                if max(v for k, v in run['checks'].items() if k.startswith('max_relative')) > .01:
                    raise ValueError(f'{path}: moment reconstruction differs by more than 1%')
                end = run['end']
                run['cells'] = cell_counts(path, run['nz'])
                run['baseline'] = uid == 'baseline'
                run['uid'] = uid
                runs[uid] = run
            run = runs[uid]
            resolution = run['cells'][AXES.index(axis)]
            if path != args.baseline and resolution != int(path.name):
                raise ValueError(f'{path}: directory label and saved resolution disagree')
            scans[axis].append(dict(uid=uid, resolution=resolution))
        scans[axis].sort(key=lambda r: r['resolution'])
        if len(scans[axis]) < 2:
            raise ValueError(f'{axis}: at least two completed runs required')
        if len({r['resolution'] for r in scans[axis]}) != len(scans[axis]):
            raise ValueError(f'{axis}: duplicate resolution; distinguish runs explicitly')
        fixed = [tuple(v for j, v in enumerate(runs[s['uid']]['cells'])
                       if j != AXES.index(axis)) for s in scans[axis]]
        if len(set(fixed)) != 1:
            raise ValueError(f'{axis}: unscanned cell counts differ')
    manifest = dict(frame=args.frame, electron_temperature_ev=args.electron_temperature_ev,
                    root=str(args.root), baseline=str(args.baseline),
                    end=end, scans=scans, unavailable=missing, runs={})
    arrays = {}
    for uid, run in runs.items():
        manifest['runs'][uid] = dict(path=str(run['path']), cells=run['cells'],
                                    baseline=run['baseline'], checks=run['checks'],
                                    changeset=run['changeset'], areas=run['areas'],
                                    metrics=list(run['data']))
        for key, (t, v) in run['data'].items():
            arrays[f'{uid}__{key}'] = np.vstack([t, v])
        arrays[f'{uid}__profile_z'], arrays[f'{uid}__profile_values'] = run['profile']
    np.savez_compressed(args.output / 'diagnostic_snapshot.npz', **arrays)
    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest, runs


def restore(args):
    manifest = json.loads((args.output / 'manifest.json').read_text())
    for key in ('frame', 'electron_temperature_ev', 'root', 'baseline'):
        expected = str(getattr(args, key)) if key in ('root', 'baseline') else getattr(args, key)
        if manifest[key] != expected:
            raise ValueError(f'Snapshot {key} differs; rerun without --reuse-data')
    runs = {}
    with np.load(args.output / 'diagnostic_snapshot.npz', allow_pickle=False) as data:
        for uid, meta in manifest['runs'].items():
            runs[uid] = dict(meta, uid=uid, end=manifest['end'],
                            data={k: tuple(data[f'{uid}__{k}']) for k in meta['metrics']},
                            profile=(data[f'{uid}__profile_z'], data[f'{uid}__profile_values']))
    return manifest, runs


def statistics(args, manifest, runs):
    end = manifest['end']
    start = end - args.window_us * 1e-6
    rows, lookup = [], {}
    for axis in AXES:
        scan = manifest['scans'][axis]
        reference = runs[scan[-1]['uid']]
        for member in scan:
            run = runs[member['uid']]
            for key, data in run['data'].items():
                mean, low, high = cv.window_stats(data, start, end)
                final = float(cv.sample(data, end))
                refmean = cv.window_stats(reference['data'][key], start, end)[0]
                reffinal = float(cv.sample(reference['data'][key], end))
                if not np.isfinite([mean, low, high, final, refmean, reffinal]).all():
                    raise ValueError(f'Non-finite statistics: {axis}/{member}/{key}')
                if not low - abs(mean)*1e-12 <= mean <= high + abs(mean)*1e-12:
                    raise ValueError('Time mean lies outside temporal range')
                row = dict(axis=axis, resolution=member['resolution'], uid=member['uid'],
                           run=str(run['path']), metric=key, mean=mean, final=final,
                           window_min=low, window_max=high,
                           reference_resolution=scan[-1]['resolution'], reference_mean=refmean,
                           mean_difference_percent=100*(mean-refmean)/abs(refmean) if refmean else np.nan,
                           final_difference_percent=100*(final-reffinal)/abs(reffinal) if reffinal else np.nan,
                           window_start_s=start, window_end_s=end)
                rows.append(row)
                lookup[axis, member['uid'], key] = row
    write_csv(args.output / 'all_diagnostics.csv', rows)
    return lookup


def overview_reference(axis, scan):
    """Keep Nmu=32 as the overview reference even when finer runs exist."""
    if axis == 'mu':
        for member in scan:
            if member['resolution'] == 32:
                return member
        raise ValueError('Convergence overview requires a completed Nmu=32 reference run')
    return scan[-1]


def main_figure(args, manifest, stats):
    fig, axs = plt.subplots(2, 3, figsize=(13, 9))
    fig.subplots_adjust(left=.10, right=.985, bottom=.10, top=.935, hspace=.25, wspace=.35)
    fig.supylabel(r'$(Q-Q_{\rm ref})/|Q_{\rm ref}|$ [\%]', x=.018, y=.5175,
                  fontsize=plt.rcParams['axes.labelsize'])
    labels = (r'$N_i$', r'$E_i$', r'$n_i(z=0)$',
              r'$\Gamma_i\,/\!\int S_i\,d^3x\,d^3v$', r'$q_i$')
    handles = [Line2D([], [], **style, label=label)
               for style, label in zip(OVERVIEW_STYLES, labels)]
    fig.legend(handles=handles, ncol=len(labels), loc='upper center', bbox_to_anchor=(.5, 1.005),
               columnspacing=1.6, handlelength=2.8)
    exported = []
    for col, axis in enumerate(AXES):
        scan = manifest['scans'][axis]
        reference = overview_reference(axis, scan)
        reference_index = scan.index(reference)
        x = np.array([s['resolution'] for s in scan])
        values = final_values(axis, scan, stats)
        for row in range(2):
            ax = axs[row, col]
            paper_axis(ax)
            ax.axhline(0, color='#60666b', lw=.8)
            for j, (key, _, _, _) in enumerate(ABS_ROWS[row*3:row*3+3]):
                final = values[key]
                ref = final[reference_index]
                if ref == 0:
                    raise ValueError(f'{axis}/{key}: zero final reference for relative difference')
                y = 100 * (final-ref) / abs(ref)
                ax.plot(x, y, **OVERVIEW_STYLES[row*3+j], zorder=3)
                for member, value, difference in zip(scan, final, y):
                    exported.append(dict(axis=axis, resolution=member['resolution'],
                                         uid=member['uid'], metric=key,
                                         final_time_s=manifest['end'], final=float(value),
                                         reference_resolution=reference['resolution'],
                                         reference_final=float(ref),
                                         final_difference_percent=float(difference)))
            ax.set_xticks(x, [str(n) for n in x])
            ax.margins(x=.09, y=.15)
            if row == 1:
                ax.set_xlabel(f'${cv.AXES[axis][1]}$ [cells]')
                ax.set_yscale('linear')
    write_csv(args.output / '01_final_state_convergence.csv', exported)
    finish(fig, args, '01_convergence_overview')


def final_values(axis, scan, stats):
    """Use identical ion diagnostics in all final-state scalar figures."""
    return {key: np.array([
        sum(weight * stats[axis, member['uid'], metric]['final'] for metric, weight in metrics)
        for member in scan]) for key, _, _, metrics in ABS_ROWS}


def absolute_figure(args, manifest, stats):
    fig, axs = plt.subplots(len(ABS_ROWS), 3, figsize=(13, 3 * len(ABS_ROWS)), sharey='row')
    fig.subplots_adjust(left=.15, right=.985, bottom=.06, top=.94, hspace=.36, wspace=.17)
    exported = []
    for col, axis in enumerate(AXES):
        scan = manifest['scans'][axis]
        x = [s['resolution'] for s in scan]
        values = final_values(axis, scan, stats)
        for row, (output_key, ylabel, scale, _) in enumerate(ABS_ROWS):
            ax = axs[row, col]
            # Read final values only; no temporal averages, ranges or electron terms.
            ax.plot(x, values[output_key]*scale, color=COLORS[col], marker='o',
                    ms=4, lw=1.3)
            paper_axis(ax)
            if col == 0:
                ax.set_ylabel(ylabel)
            # Completed runs only, with plain numeric tick labels.
            ax.set_xticks(x, [str(n) for n in x])
            ax.margins(x=.09, y=.15)
            if row == len(ABS_ROWS) - 1:
                ax.set_xlabel(f'${cv.AXES[axis][1]}$ [cells]')
        for i, s in enumerate(scan):
            exported.append(dict(axis=axis, resolution=s['resolution'], uid=s['uid'],
                                 final_time_s=manifest['end'],
                                 **{key: float(v[i]) for key, v in values.items()},
                                 **{'Gamma_i_s-1': sum(stats[axis, s['uid'], key]['final']
                                                      for key in ('loss_lo', 'loss_hi')),
                                    'source_i_s-1': stats[axis, s['uid'], 'source']['final']}))
    # Show source-loss balance explicitly on the shared ratio axis.
    ratios = [r['Gamma_i_over_source'] for r in exported]
    lower = min(.4, min(ratios) - .05)
    upper = max(1.05, max(ratios) + .05)
    ratio_ax = axs[3, 0]
    ratio_ax.set_ylim(lower, upper)
    ticks = np.arange(np.ceil(lower / .2), np.floor(upper / .2) + 1) * .2
    ratio_ax.set_yticks(ticks)
    write_csv(args.output / '02_final_state_data.csv', exported)
    finish(fig, args, '02_absolute_diagnostics')


def particle_balance_figure(args, manifest, stats):
    """Isolate the absolute source-loss ratio and compare every run with unity."""
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.4), sharey=True)
    fig.subplots_adjust(left=.09, right=.89, bottom=.20, top=.85, wspace=.17)
    balance_style = dict(color='#60666b', ls='--', lw=1.2)
    baseline_style = dict(color='#20272c', marker='D', mfc='white',
                          mew=1.3, ms=7, ls='none')
    fig.legend(handles=[
        Line2D([], [], **balance_style, label=r'Expected ratio $=1$'),
        Line2D([], [], **baseline_style, label=r'Baseline $(400,64,32)$'),
    ], loc='upper center', bbox_to_anchor=(.49, 1), ncol=2,
       columnspacing=2, handlelength=2.8)
    exported = []
    for col, axis in enumerate(AXES):
        ax = axs[col]
        scan = manifest['scans'][axis]
        x = np.array([member['resolution'] for member in scan])
        ratios = final_values(axis, scan, stats)['Gamma_i_over_source']
        paper_axis(ax)
        ax.axhline(1, **balance_style)
        ax.plot(x, ratios, color=COLORS[col], marker='o', ms=5, lw=1.3, zorder=3)
        for member, ratio in zip(scan, ratios):
            if manifest['runs'][member['uid']]['baseline']:
                ax.plot(member['resolution'], ratio, **baseline_style, zorder=4)
            exported.append(dict(
                axis=axis, resolution=member['resolution'], uid=member['uid'],
                final_time_s=manifest['end'], Gamma_i_over_source=float(ratio),
                deviation_from_unity_percent=float(100 * (ratio - 1)),
                **{'Gamma_i_s-1': sum(stats[axis, member['uid'], key]['final']
                                     for key in ('loss_lo', 'loss_hi')),
                   'source_i_s-1': stats[axis, member['uid'], 'source']['final']}))
        ax.set_xticks(x, [str(n) for n in x])
        ax.set_xlabel(f'${cv.AXES[axis][1]}$ [cells]')
        ax.margins(x=.09)
    ylabel = next(label for key, label, _, _ in ABS_ROWS if key == 'Gamma_i_over_source')
    axs[0].set_ylabel(ylabel)
    # Pool every scan and include unity, so no run is clipped or rescaled by
    # a numerical reference. The right axis is the same data in percent.
    ratios = [row['Gamma_i_over_source'] for row in exported]
    lower, upper = min(.4, min(ratios) - .05), max(1.05, max(ratios) + .05)
    ticks = np.arange(np.ceil(lower / .1), np.floor(upper / .1) + 1) * .1
    axs[0].set_ylim(lower, upper)
    axs[0].set_yticks(ticks)
    deviation = axs[-1].secondary_yaxis(
        'right', functions=(lambda ratio: 100 * (ratio - 1),
                            lambda percent: 1 + percent / 100))
    deviation.set_yticks(np.round(100 * (ticks - 1), 10))
    deviation.set_ylabel(r'Deviation from unity [\%]')
    deviation.spines['right'].set_visible(True)
    write_csv(args.output / '05_particle_balance_data.csv', exported)
    finish(fig, args, '05_particle_balance')


def run_style(manifest, axis, member):
    scan = manifest['scans'][axis]
    if len(scan) > len(SCAN_COLORS):
        raise ValueError('Extend SCAN_COLORS to give every resolution within a scan a distinct color')
    i = scan.index(member)
    baseline = manifest['runs'][member['uid']]['baseline']
    return dict(color=SCAN_COLORS[i], ls=('-', '--', '-.', ':', '-')[i % 5],
                lw=1.7 if baseline else 1.3)


def run_legend(manifest, axis):
    return [Line2D([], [], **run_style(manifest, axis, member),
                   label=f'${cv.AXES[axis][1]}={member["resolution"]}$')
            for member in manifest['scans'][axis]]


def histories_figure(args, manifest, runs):
    fig, axs = plt.subplots(len(ABS_ROWS), 3, figsize=(13, 3 * len(ABS_ROWS)), sharex=True, sharey='row')
    fig.subplots_adjust(left=.13, right=.985, bottom=.06, top=.93, hspace=.28, wspace=.17)
    short_labels = {'Ni': r'$N_i$ [$10^{19}$]', 'Ei_J': r'$E_i$ [kJ]'}
    end = manifest['end']
    start = end-args.zoom_us*1e-6
    for col, axis in enumerate(AXES):
        scan = manifest['scans'][axis]
        center = (axs[0, col].get_position().x0 + axs[0, col].get_position().x1) / 2
        fig.legend(handles=run_legend(manifest, axis), loc='upper center',
                   bbox_to_anchor=(center, .99), ncol=2,
                   columnspacing=.8, handlelength=1.8, handletextpad=.5)
        for row, (key, ylabel, scale, metrics) in enumerate(ABS_ROWS):
            ax = axs[row, col]
            paper_axis(ax)
            for member in scan:
                times, values = history_values(runs[member['uid']]['data'], metrics, start, end)
                ax.plot((times-end)*1e6, values*scale, **run_style(manifest, axis, member))
            if key == 'Gamma_i_over_source':
                ax.axhline(1, color='#60666b', ls='--', lw=1.2)
            ax.set_xlim(-args.zoom_us, 0)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.margins(y=.15)
            if col == 0:
                ax.set_ylabel(short_labels.get(key, ylabel))
            if row == len(ABS_ROWS) - 1:
                ax.set_xlabel(r'$t-t_{\rm end}$ [$\mu$s]')
    finish(fig, args, '03_final_relaxation')


def profiles_figure(args, manifest, runs):
    fig = plt.figure(figsize=(16, 14))
    grid = fig.add_gridspec(4, 3, left=.09, right=.95, bottom=.08, top=.90,
                           hspace=.24, wspace=.45)
    axs = np.empty((4, 3, 2), dtype=object)
    labels = (r'$n_i$ [m$^{-3}$]', r'$u_\parallel$ [$10^6$ m s$^{-1}$]',
              r'$T_\parallel$ [keV]', r'$T_\perp$ [keV]')
    scales = (1, 1e-6, 1, 1)
    throat_style = dict(color='#888888', ls='--', lw=1, zorder=10)
    extent = 1.015 * max(np.max(np.abs(run['profile'][0])) for run in runs.values())
    for col, axis in enumerate(AXES):
        for row in range(4):
            pair = grid[row, col].subgridspec(1, 2, wspace=0)
            for side in range(2):
                ax = fig.add_subplot(pair[0, side],
                                     sharex=axs[0, col, side] if row else None,
                                     sharey=axs[row, 0, side] if col else None)
                axs[row, col, side] = ax
                paper_axis(ax)
                if side:
                    ax.set_yscale('log', nonpositive='mask')
                    ax.yaxis.tick_right()
                    ax.spines['right'].set_visible(True)
                    ax.spines['left'].set_visible(False)
                else:
                    ax.spines['right'].set_visible(True)
                    ax.spines['right'].set_color('#c3c8cc')
                for member in manifest['scans'][axis]:
                    z, v = runs[member['uid']]['profile']
                    values = v[:, row] * scales[row]
                    keep = (z >= 0) & (values > 0) if side else z < 0
                    # Mask, rather than delete, samples so omitted values and DG
                    # cell faces remain gaps instead of acquiring connecting lines.
                    values = np.where(keep, values, np.nan)
                    zz = np.column_stack([z.reshape(-1, 2),
                                          np.full(len(z)//2, np.nan)]).ravel()
                    yy = np.column_stack([values.reshape(-1, 2),
                                          np.full(len(z)//2, np.nan)]).ravel()
                    ax.plot(zz, yy, **run_style(manifest, axis, member))
                ax.axvline(.98 if side else -.98, **throat_style)
                if side:
                    ax.set_xlim(0, extent)
                    ax.set_xticks([1, 2])
                else:
                    ax.set_xlim(-extent, 0)
                    ax.set_xticks([-2, -1, 0])
                    if row in (2, 3):
                        ax.axhline(0, color='#8e3030', ls=':', lw=.8)
                    if col == 0:
                        ax.set_ylabel(labels[row])
                ax.tick_params(axis='x', labelbottom=row == 3)
                ax.margins(y=.08)
            if row == 3:
                axs[row, col, 1].set_xlabel(r'$z$ [m]')
                axs[row, col, 1].xaxis.set_label_coords(0, -.18)
        center = (axs[0, col, 0].get_position().x0 + axs[0, col, 1].get_position().x1) / 2
        handles = run_legend(manifest, axis)
        if col == 0:
            handles.append(Line2D([], [], **throat_style, label=r'$Z=\pm Z_m$'))
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(center, .985),
                   ncol=2, columnspacing=.8, handlelength=1.8, handletextpad=.5)
    axs[0, 0, 0].set_ylim(bottom=0)
    # Estimate display bounds from unique runs, so the baseline is not counted
    # once per scan. Keep all plotted samples; only the viewport clips outliers.
    for row, scale in enumerate(scales):
        positive = np.concatenate([v[(z >= 0) & (v[:, row] > 0), row] * scale
                                   for z, v in (run['profile'] for run in runs.values())])
        if not positive.size:
            raise ValueError(f'Profile component {row} has no positive right-half samples')
        if row == 3:
            # Include the full positive perpendicular-temperature exhaust tail,
            # with padding in log space.
            log_pad = (positive.max() / positive.min()) ** .08
            log_low, log_high = positive.min() / log_pad, positive.max() * log_pad
        else:
            log_low = 10.0 ** np.floor(np.log10(np.quantile(positive, .05)))
            log_high = positive.max() * (positive.max() / log_low) ** .08
        axs[row, 0, 1].set_ylim(log_low, log_high)
        if row in (2, 3):
            left = np.concatenate([v[z < 0, row] * scale
                                   for z, v in (run['profile'] for run in runs.values())])
            low, high = min(0, np.quantile(left, .01)), left.max()
            pad = .08 * (high - low)
            axs[row, 0, 0].set_ylim(low - pad, high + pad)
    # Retain red shading for the visible negative-temperature region.
    for ax in axs[2:, :, 0].flat:
        low, high = ax.get_ylim()
        if low < 0:
            ax.axhspan(low, 0, color='#f5e5e3', zorder=0)
            ax.set_ylim(low, high)
    finish(fig, args, '04_final_profiles')


def report(args, manifest, runs, stats):
    lines = [
        'PUBLICATION CONVERGENCE FIGURES',
        f'Frame {args.frame}; comparison time {manifest["end"]:.12g} s.',
        f'all_diagnostics.csv averaging window: final {args.window_us:g} microseconds. '
        'Figures 1, 2, 4 and 5 use only the final comparison state.',
        f'{len(runs)} unique completed simulations; the baseline is shared by all scans.',
        'Diagnostic definitions, interpolation, reconstruction checks and normalizations are from figures/diagnostics.py.',
        'Data read from raw Gkeyll output, not from older convergence-plots CSV files.',
        'Snapshot reuse is explicit (--reuse-data); omit it to incorporate new simulation output.',
        '', 'FIGURE CAPTIONS',
        f'1. Final-state ion convergence at frame {args.frame}, t={manifest["end"]:.9g} s. '
        'The five quantities and their definitions match figure 2: Ni, Ei, ni(0), '
        'Gamma_i / integral(S_i d^3x d^3v), and total wall heat flux qi. '
        'Gamma_i is the sum of the outward ion loss rates through both walls. '
        'All points show 100*(Q_final-Q_reference_final)/abs(Q_reference_final), without '
        'temporal averages or error bars. Each column varies Nz, Nvpar or Nmu, with the '
        'other cell counts fixed at the baseline values (400, 64, 32). References are the '
        'finest available runs for Nz and Nvpar, and Nmu=32 for the mu scan, not continuum solutions; zero at the reference '
        'is definitional. The top row shows Ni, Ei and ni(0) on independent linear scales. '
        'The bottom row shows Gamma_i/source and qi on independent linear scales. '
        'Lines guide the eye; no convergence order is fitted. '
        'Only completed runs appear; the deleted Nz=128 directory is omitted. '
        '01_final_state_convergence.csv records each final value, final reference and percent difference. '
        'All five figures share typography and clean axes, '
        'with no panel annotations, subtitles, footnotes or daggers. The overview has no title, '
        'one centered shared y-label, and a compact one-row, five-column legend; each quantity '
        'has a distinct color, line style and marker. There are no master titles.',
        f'2. Final-state ion diagnostics at frame {args.frame}, t={manifest["end"]:.9g} s. '
        'Rows show ion number Ni, ion kinetic energy Ei, central ion density ni(0), boundary '
        'flux/source Gamma_i / integral(S_i d^3x d^3v), and total ion wall heat flux qi. '
        'qi sums the saved ion Hamiltonian energy fluxes integrated over the lower and upper walls, '
        'without division by wall area. The snapshot reader stores the two integrated fluxes '
        'under Pi_lo and Pi_hi; their sum appears only as qi in the figures and qi_W in plotted CSVs. '
        'All points are final values, without averaging in time or error bars. '
        'Gamma_i=loss_lower+loss_upper sums the outward ion particle losses through both walls. '
        'The denominator is component 0 of ion_source_integrated_moms, the phase-space integral '
        'of the ion source. Their ratio is dimensionless; unity denotes source-loss balance '
        'and is explicitly included in the axis limits and ticks. '
        'The SI normalization is specified by the simulation author: Ni is a particle count, '
        'Ei is in J, central density is in m^-3, and qi is in W. '
        'Axis labels give the display multipliers (including kJ and kW). '
        'Only ion quantities appear; electron and ion-plus-electron curves are omitted. '
        'Columns vary Nz, Nvpar and Nmu respectively; unscanned cell counts are fixed at '
        'the baseline values (400, 64, 32). Rows share vertical scales. '
        'All completed simulations are included, with the shared baseline in each scan. '
        'The deleted Nz=128 directory is not plotted. '
        '02_final_state_data.csv records every plotted value without display scaling, along with '
        'the separate particle-loss numerator and source denominator. '
        'Ni=integral(n_i d^3x), and Ei=(m_i/2) integral(M2_i d^3x), where '
        'M2_i=integral[(v_parallel^2+2 mu B/m_i) f_i d^3v]. '
        'Thus Ei=integral[(m_i v_parallel^2/2+mu B) f_i d^3x d^3v], including '
        'bulk parallel-flow energy and parallel/perpendicular thermal energy. '
        'The corresponding continuum moment expression is '
        'integral[n_i(m_i u_parallel^2/2+T_parallel/2+T_perp) d^3x], with temperatures in J. '
        'The plotted Ei is reconstructed from ion_M2 and geo_int_jacobgeo using the exact '
        'modal product integral for the saved p=1 expansions. Ni uses component 0 of '
        'ion_integrated_moms and is independently checked against the M0/Jacobian integral. '
        'The configured Hamiltonian diagnostic has components '
        '[integral(f_i), integral(m_i v_parallel f_i), '
        'integral((m_i v_parallel^2/2+mu B+e_i phi) f_i)] over phase space, '
        'so its third component is Ei+integral(e_i n_i phi d^3x). '
        'That potential-energy term is excluded from the plotted kinetic energy. '
        'Here e_i denotes the ion charge, to distinguish it from the heat flux qi.',
        f'3. Final-relaxation histories over the last {args.zoom_us:g} microseconds, with the same '
        'five absolute ion quantities, row order and display units as figure 2: Ni in 10^19, '
        'Ei in kJ, ni(0) in 10^19 m^-3, Gamma_i / integral(S_i d^3x d^3v), '
        'and total wall heat flux qi in kW. '
        'The source-loss ratio sums the outward ion losses through both walls; a dashed line '
        'marks unity. There is no normalization by another simulation or temporal mean, '
        'and no averaging-window shading. Time is t-t_end in microseconds, with zero at the '
        'simulation end. Rows share vertical scales across scans. '
        'High-cadence diagnostics retain their native samples; energy and central density '
        'use linear interpolation between saved times. Combined wall diagnostics use the '
        'union of both sample grids. Compact mathematical resolution legends and colors '
        'match figure 4, with the shared baseline purple. There are no titles, subtitles, '
        'panel annotations, daggers or bottom explanatory text. The OAP/FDP clock '
        'is not interpreted as an uninterrupted physical transient. The full histories remain '
        'available in diagnostic_snapshot.npz and the original per-scan plots.',
        '4. Final density, parallel flow, parallel temperature and perpendicular temperature versus '
        'mapped physical position. The two interior samples per DG cell match diagnostics.py. '
        'Cell traces are drawn separately. Every quantity is split at z=0: the left half '
        '(z<0) uses a linear vertical scale and the right half (z>=0) a logarithmic scale. '
        'Nonpositive samples are masked on the logarithmic halves without connecting across gaps. '
        'The perpendicular-temperature logarithmic halves show the full positive data range '
        'with 8% padding in log space. Other logarithmic lower limits use the fifth percentile '
        'of positive right-half samples, rounded down to a power of ten. '
        'Linear temperature lower limits use the first '
        'percentile of left-half samples (or zero, if smaller), with 8% range padding. '
        'These bounds pool each distinct run once and are shared across scans; outliers '
        'are clipped by the display limits, not removed from the data. Upper bounds include '
        'all positive peaks with 8% padding in the corresponding axis scale. '
        'Density is in m^-3, flow in 10^6 m s^-1 and temperatures in keV. '
        'Only the final comparison frame is shown, without temporal averaging or error bars. '
        'Scan legends show completed resolutions, restarting the same palette in increasing '
        'resolution order in each column, with the production baseline purple in all columns. Grey dashed '
        'mirror-throat markers at z=+/-0.98 m are drawn above all profile curves, with '
        'one shared Z=+/-Zm entry in the Nz legend. There are no daggers, unavailable entries, '
        'figure title or inline throat labels. '
        'Negative temperatures in either component are retained and shaded on the linear halves; the numerical '
        'checks remain in quality_checks.csv. Positivity between sampled positions is not '
        'certified by this check.',
        f'5. Final-state ion source-loss balance at frame {args.frame}, t={manifest["end"]:.9g} s. '
        'A single row isolates the Gamma_i / integral(S_i d^3x d^3v) row of figure 2, '
        'with columns varying Nz, Nvpar and Nmu and identical shared vertical limits. '
        'All completed scan members are shown at their actual cell counts, including '
        'Nmu=40; the shared baseline (400, 64, 32) is marked by an open diamond in '
        'each column. The dashed horizontal line marks the expected ratio of 1. '
        'The right-hand axis expresses the same ratio B as the signed percentage '
        'deviation 100*(B-1), with zero denoting source-loss balance. '
        'The numerator sums outward ion losses through both walls and the denominator '
        'is the phase-space integral of the ion source, exactly as in figure 2. '
        'Values are taken at the final comparison time, without temporal averaging, '
        'error bars, or normalization by a finite-grid reference. Lines guide the eye. '
        '05_particle_balance_data.csv records all ratios, signed deviations from unity, '
        'particle-loss numerators and source denominators.',
        '', 'SCAN MEMBERSHIP AND FINEST-PAIR DIFFERENCES',
    ]
    for axis in AXES:
        scan = manifest['scans'][axis]
        lines.append(f'{axis}: ' + ', '.join(str(s['resolution']) for s in scan)
                     + f'; overview reference = {overview_reference(axis, scan)["resolution"]}'
                     + f'; all_diagnostics.csv averaging reference = {scan[-1]["resolution"]}.')
        previous = scan[-2]
        for key, final in final_values(axis, scan, stats).items():
            difference = 100 * (final[-2]-final[-1]) / abs(final[-1]) if final[-1] else np.nan
            lines.append(f'  {key}: {previous["resolution"]} vs {scan[-1]["resolution"]}, '
                         f'final difference {difference:+.6g}%.')
    lines.extend(['', 'RUN CHECKS (two sampled positions per cell)'])
    quality = []
    for uid, run in runs.items():
        z, v = run['profile']
        row = dict(uid=uid, path=str(run['path']), Nz=run['cells'][0],
                   Nvpar=run['cells'][1], Nmu=run['cells'][2],
                   samples=len(z), min_Tpar_keV=float(v[:, 2].min()),
                   min_Tperp_keV=float(v[:, 3].min()), **run['checks'])
        row['negative_final_Tperp_samples'] = int(np.sum(v[:, 3] < 0))
        quality.append(row)
        lines.append(f'{uid}: cells={run["cells"]}, min Tpar={row["min_Tpar_keV"]:.8g} keV; '
                     f'{run["checks"]["negative_final_Tpar_samples"]}/{len(z)} negative Tpar samples; '
                     f'min Tperp={row["min_Tperp_keV"]:.8g} keV, '
                     f'{row["negative_final_Tperp_samples"]}/{len(z)} negative Tperp samples; '
                     f'code={run["changeset"]}.')
    write_csv(args.output / 'quality_checks.csv', quality)
    lines.extend(['', 'UNAVAILABLE RUNS'])
    for entry in manifest['unavailable']:
        lines.append(f'{entry["axis"]}={entry["resolution"]}: {entry["reason"]}; {entry["path"]}')
    lines.extend(['', 'REPRODUCTION', 'python publication_convergence.py',
                  'python publication_convergence.py --reuse-data  # explicitly use the saved snapshot',
                  'Outputs: vector PDF and SVG, 300 dpi PNG (configurable with --dpi).',
                  f'all_diagnostics.csv includes all {len(next(iter(runs.values()))["data"])} exported diagnostic quantities for every scan member.',
                  'manifest.json records input paths, grid sizes, code changesets, checks and missing runs.'])
    lines.extend(['', 'GKEYLL SOURCE VERIFICATION',
                  'Checked local /home/mr1884/gkeyll at 1dd609d0006c5bc2a0bc1bddb888d446ef2ab98c, '
                  'during the original resolution analysis. Current per-run changesets are recorded in manifest.json; '
                  'the current set can contain multiple revisions.',
                  'gyrokinetic/ker/gyrokinetic_mom/gyrokinetic_int_mom_1x2v_ser_p1.c: '
                  'lines 2-8 (integrated M0), 40-53 (integrated M2), 111-125 (Hamiltonian moments).',
                  'gyrokinetic/ker/gyrokinetic_mom/gyrokinetic_mom_1x2v_ser_p1.c: '
                  'lines 16-29 (saved M2 kernel).',
                  'gyrokinetic/apps/gk_species_moment.c: lines 26-38 and 89-96 '
                  '(geometric Jacobian removed from nonintegrated diagnostics).',
                  'gyrokinetic/apps/gk_species.c: lines 898-922 '
                  '(calculate cell integrals, sum over the domain and MPI ranks).',
                  'gyrokinetic/apps/gk_species_source.c: lines 146-175 '
                  '(same moment integration applied to the phase-space source).',
                  'gyrokinetic/zero/gk_maxwellian_moments.c: lines 157-193 '
                  '(parallel/perpendicular temperature definitions and bulk-flow subtraction).',
                  'figures/diagnostics.py: load_run (kinetic-energy reconstruction and independent Hamiltonian check).'])
    (args.output / 'captions_and_methods.txt').write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--baseline', type=Path, default=ROOT.parent / '1x-beams')
    parser.add_argument('--output', type=Path, default=SCRIPT_DIR / 'publication-figures')
    parser.add_argument('--frame', type=int, default=75)
    parser.add_argument('--window-us', type=float, default=15)
    parser.add_argument('--zoom-us', type=float, default=60)
    parser.add_argument('--electron-temperature-ev', type=float, default=940)
    parser.add_argument('--dpi', type=int, default=300)
    parser.add_argument('--reuse-data', action='store_true')
    parser.add_argument('--figures', type=int, nargs='+', choices=(1, 2, 3, 4, 5),
                        default=[1, 2, 3, 4, 5], help='Figures to render (default: all five)')
    args = parser.parse_args()
    if args.frame < 1 or min(args.window_us, args.zoom_us, args.electron_temperature_ev, args.dpi) <= 0:
        parser.error('Frame, windows, electron temperature and dpi must be positive')
    if args.window_us > args.zoom_us:
        parser.error('Averaging window must fit inside the displayed relaxation window')
    args.root, args.baseline, args.output = args.root.resolve(), args.baseline.resolve(), args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest, runs = restore(args) if args.reuse_data else collect(args)
    stats = statistics(args, manifest, runs)
    style()
    if 1 in args.figures:
        main_figure(args, manifest, stats)
    if 2 in args.figures:
        absolute_figure(args, manifest, stats)
    if 3 in args.figures:
        histories_figure(args, manifest, runs)
    if 4 in args.figures:
        profiles_figure(args, manifest, runs)
    if 5 in args.figures:
        particle_balance_figure(args, manifest, stats)
    report(args, manifest, runs, stats)
    print(f'Completed figures {sorted(set(args.figures))} for {len(runs)} distinct simulations in {args.output}')


if __name__ == '__main__':
    main()
