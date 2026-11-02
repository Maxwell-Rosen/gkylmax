"""Publication figures for alpha, cycle/FDP-length, and time-dilation studies.

Times are saved simulation-clock times. Final values use each run's own endpoint;
incomplete alpha cases are excluded, while the long continuation is explicitly
marked provisional. No conversion to a physical time is inferred from a CFL cap.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import operator
import re
from datetime import datetime, timezone
from pathlib import Path

from publication_style import ABS_ROWS, diagnostic_panels, finish, history_values, paper_axis, style
import diagnostics as cv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np

SCHEMA = 1
PARAMETERS = ('alpha_oap', 'alpha_fdp', 'tau_oap', 'tau_fdp', 'tau_fdp_extra',
              'num_cycles', 'num_frames_oap', 'num_frames_fdp',
              'num_frames_fdp_extra', 'cfl_factor_times_omega_max', 'Z_m')
COLORS = ('#6A3D9A', '#0072B2', '#D55E00', '#009E73', '#8C564B')
LINESTYLES = ('-', '--', '-.', ':', '-')
SHORT_LABELS = { 'Ni': r'$N_i$ [$10^{19}$]', 'Ei_J': r'$E_i$ [kJ]' }
# Reviewed, comment/whitespace-independent phase-assembly blocks. Checking the
# entire block retains the damping and multiplier checks despite different C
# loop/index bookkeeping in the annealing input. Schedule-code changes require
# renewed review; scalar durations/counts outside these blocks are parsed below.
PHASE_ASSEMBLY_HASHES = {
    'standard': '154f8fddcb24fb198c282fa4b09ed39f554cdab8cf8293f4aa8de914d83d0891',
    'annealing': 'c4e16581a039c7e5e45ef52ba47474c9ca0913a6c68602714e0cbea2e4aab07f',
}


def write_csv(path, rows):
    if not rows:
        raise ValueError(f'No rows for {path}')
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def arithmetic(expression):
    """Evaluate numeric C-style scalar expressions without executing input code."""
    binary = {ast.Add: operator.add, ast.Sub: operator.sub,
              ast.Mult: operator.mul, ast.Div: operator.truediv}
    unary = {ast.UAdd: operator.pos, ast.USub: operator.neg}

    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return node.value
        if isinstance(node, ast.BinOp) and type(node.op) in binary:
            return binary[type(node.op)](visit(node.left), visit(node.right))
        if isinstance(node, ast.UnaryOp) and type(node.op) in unary:
            return unary[type(node.op)](visit(node.operand))
        raise ValueError(f'Unsupported scalar expression: {expression}')

    result = float(visit(ast.parse(expression, mode='eval').body))
    if not np.isfinite(result):
        raise ValueError(f'Nonfinite scalar: {expression}')
    return result


def source_text(path):
    return re.sub(r'//[^\n]*|/\*.*?\*/', '', path.read_text(), flags=re.S)


def read_parameters(path):
    source = source_text(path / 'sim.c')
    parameters = {}
    for key in PARAMETERS:
        match = re.search(rf'\b(?:double|int)\s+{key}\s*=\s*([^;]+);', source)
        if match is None:
            raise ValueError(f'{path}: missing phase parameter {key}')
        parameters[key] = arithmetic(match[1])
    for key in PARAMETERS:
        if key.startswith('num_'):
            value = parameters[key]
            if value < 0 or not value.is_integer():
                raise ValueError(f'{path}: invalid count {key}={value}')
            parameters[key] = int(value)
    p = parameters
    if min(p['tau_oap'], p['tau_fdp']) < 0 or p['tau_fdp_extra'] <= 0:
        raise ValueError(f'{path}: invalid phase durations')
    if p['alpha_fdp'] != 1:
        raise ValueError(f'{path}: expected alpha_fdp=1')
    p['expected_frame'] = p['num_cycles'] * (p['num_frames_oap'] + p['num_frames_fdp']) + p['num_frames_fdp_extra']
    p['expected_end_s'] = p['num_cycles'] * (p['tau_oap'] + p['tau_fdp']) + p['tau_fdp_extra']
    p['final_phase_start_s'] = p['num_cycles'] * (p['tau_oap'] + p['tau_fdp'])
    p['continuous_fdp_start_s'] = p['final_phase_start_s'] - (p['tau_fdp'] if p['num_cycles'] else 0)
    p['final_phase_start_frame'] = p['num_cycles'] * (p['num_frames_oap'] + p['num_frames_fdp'])
    p['continuous_fdp_start_frame'] = p['final_phase_start_frame'] - (p['num_frames_fdp'] if p['num_cycles'] else 0)
    if 'annealing_alphas' in source:
        arrays = {}
        for key in ('annealing_alphas', 'annealing_durations'):
            match = re.search(rf'\b{key}\[\]\s*=\s*\{{([^}}]+)\}};', source)
            if match is None:
                raise ValueError(f'{path}: missing {key} array')
            expressions = re.sub(r'\btau_oap\b', str(p['tau_oap']), match[1])
            arrays[key] = [arithmetic(v.strip()) for v in expressions.split(',') if v.strip()]
        alphas, durations = arrays['annealing_alphas'], arrays['annealing_durations']
        if not alphas or len(alphas) != len(durations) or min(alphas + durations) <= 0:
            raise ValueError(f'{path}: invalid annealing stages')
        start, first = p['final_phase_start_s'], p['final_phase_start_frame']
        stages = []
        for alpha, duration in zip(alphas, durations):
            stages.append(dict(kind='OAP', alpha=alpha, duration_s=duration,
                               start_s=start, end_s=start + duration,
                               start_frame=first, end_frame=first + p['num_frames_oap']))
            start += duration
            first += p['num_frames_oap']
        p.update(final_phase_start_s=start, continuous_fdp_start_s=start,
                 final_phase_start_frame=first, continuous_fdp_start_frame=first,
                 expected_end_s=start + p['tau_fdp_extra'],
                 expected_frame=first + p['num_frames_fdp_extra'])
        stages.append(dict(kind='FDP', alpha=p['alpha_fdp'], duration_s=p['tau_fdp_extra'],
                           start_s=start, end_s=p['expected_end_s'],
                           start_frame=first, end_frame=p['expected_frame']))
        p['annealing_stages'] = stages
    imported = re.search(r'\.file_name\s*=\s*"([^"]+)"', source)
    p['imported_distribution'] = imported[1] if imported else None
    return p


def physics_signature(path, study):
    """Allow only the reviewed input changes within each comparison."""
    source = source_text(path / 'sim.c')
    source = re.sub(r'(\.filename_psi\s*=\s*)"([^"]+)"',
                    lambda m: f'{m[1]}"{Path(m[2]).name}"', source)
    allowed = ['alpha_oap'] if study == 'alpha-scan' else []
    if study == 'extra-cycles-fdp-length':
        allowed += ['num_cycles', 'tau_oap', 'tau_fdp', 'tau_fdp_extra', 'num_frames_oap',
                    'num_frames_fdp', 'num_frames_fdp_extra']
        # File-imported ions replace only the initial ion projection.
        source = re.sub(r'\.(?:projection|init_from_file)\s*=\s*\{[^{}]*\}',
                        '.initialization = {REVIEWED}', source)
        assembly = re.search(r'(?:double low_pass_damping_rate\s*=|double t_end\s*=).*?'
                             r'(?=double write_phase_freq\s*=)', source, re.S)
        kind = 'annealing' if 'annealing_alphas' in source else 'standard'
        if assembly is None or hashlib.sha256(re.sub(r'\s+', '', assembly[0]).encode()).hexdigest() != PHASE_ASSEMBLY_HASHES[kind]:
            raise ValueError(f'{path}: unreviewed {kind} phase assembly')
        source = source[:assembly.start()] + 'REVIEWED_PHASE_ASSEMBLY;' + source[assembly.end():]
    if study == 'extra-time-dilation':
        allowed += ['tau_fdp', 'tau_fdp_extra', 'cfl_factor_times_omega_max']
        # The OAP collisionless multiplier stays at 0.1 in both inputs.
        source = source.replace('.cfl_factor_times_omega_max = 1/10.0,',
                                '.cfl_factor_times_omega_max = cfl_factor_times_omega_max,')
    for key in allowed:
        source = re.sub(rf'(\b(?:int|double)\s+{key}\s*=)[^;]+;', rf'\g<1> REVIEWED;', source)
    return re.sub(r'\s+', '', source)


def candidates(root, study):
    baseline = dict(uid='baseline', path=root / '1x-beams', label='Baseline', short_label='Baseline', provisional=False)
    if study == 'alpha-scan':
        runs = [baseline]
        for path in sorted((root / '1x-beams-alpha-scan' / 'alpha-scan').glob('alpha-*')):
            if path.is_dir() and path.resolve() != baseline['path'].resolve():
                runs.append(dict(uid=path.name.replace('.', 'p'), path=path,
                                 label='', short_label='', provisional=False))
        return sorted(runs, key=lambda r: read_parameters(r['path'])['alpha_oap'])
    if study == 'extra-cycles-fdp-length':
        return [baseline,
                dict(uid='ten_cycles', path=root / '1x-beams-ecycles', label='10 cycles', short_label='10 cycles', provisional=False),
                dict(uid='fdp_300us', path=root / '1x-beams-fdp-continuation',
                     label=r'Baseline + 300 $\mu$s FDP', short_label=r'+300 $\mu$s', provisional=False),
                dict(uid='fdp_ongoing', path=root / '1x-beams-fdp-continuation-extra',
                     label='', short_label='Ongoing', provisional=True),
                dict(uid='annealing', path=root / '1x-beams-anneal',
                     label='Annealing', short_label='Annealing', provisional=False)]
    if study == 'extra-time-dilation':
        baseline.update(label=r'$c_{\omega}=10^{-1}$', short_label=r'$10^{-1}$')
        return [baseline, dict(uid='extra_dilated', path=root / '1x-beams-extra-dilated',
                               label=r'$c_{\omega}=10^{-3}$', short_label=r'$10^{-3}$', provisional=False)]
    raise ValueError(study)


def complete_frames(path):
    frames = []
    for file in path.glob('zzim-ion_BiMaxwellianMoments_*.gkyl'):
        frame = int(file.stem.rsplit('_', 1)[-1])
        if all((path / f'zzim-{key}_{frame}.gkyl').is_file()
               for key in ('ion_M0', 'ion_M2', 'field')):
            frames.append(frame)
    return sorted(frames)


def input_inventory(path, frame):
    names = ['sim.c', 'zzim-ion_integrated_moms.gkyl', 'zzim-elc_integrated_moms.gkyl',
             'zzim-ion_source_integrated_moms.gkyl', 'zzim-geo_surf0_lenr.gkyl',
             'zzim-geo_int_jacobgeo.gkyl', 'zzim-geo_corn_mc2nu_pos_deflated.gkyl']
    names += [f'zzim-ion_bflux_x{side}_integrated_HamiltonianMoments.gkyl' for side in ('lower', 'upper')]
    names += [f'zzim-{key}_{k}.gkyl' for k in range(frame + 1)
              for key in ('ion_BiMaxwellianMoments', 'ion_M0', 'ion_M2', 'field')]
    return {name: dict(size=(path / name).stat().st_size, mtime_ns=(path / name).stat().st_mtime_ns)
            for name in names}


def collect(args):
    runs, unavailable = [], []
    retained, previous = {}, None
    specs = candidates(args.root, args.study)
    if args.refresh_runs:
        unknown = set(args.refresh_runs) - {r['uid'] for r in specs}
        if unknown:
            raise ValueError(f'Unknown run IDs: {sorted(unknown)}')
        previous, old_runs = restore(args, require_current=False)
        refreshed_paths = {str(s['path'].resolve()) for s in specs if s['uid'] in args.refresh_runs}
        unavailable = [r for r in previous['unavailable'] if r['path'] not in refreshed_paths]
        retained = {r['uid']: r for r in old_runs if r['uid'] not in args.refresh_runs}
        for run in retained.values():
            run.setdefault('collected_utc', previous['collected_utc'])
        if 'baseline' in retained and retained['baseline']['input_sha256'] != hashlib.sha256((args.root / '1x-beams' / 'sim.c').read_bytes()).hexdigest():
            raise ValueError('Retained baseline input has changed; refresh baseline as well')
    reference_signature = physics_signature(args.root / '1x-beams', args.study)
    for spec in specs:
        if args.refresh_runs and spec['uid'] not in args.refresh_runs:
            if spec['uid'] not in retained:
                raise ValueError(f'{spec["uid"]}: no retained snapshot; include in --refresh-runs')
            runs.append(retained[spec['uid']])
            print(f'Keeping recorded snapshot for {spec["uid"]} at frame {retained[spec["uid"]]["frame"]}', flush=True)
            continue
        path = spec['path'].resolve()
        params = read_parameters(path)
        if physics_signature(path, args.study) != reference_signature:
            raise ValueError(f'{path}: input physics differs beyond the reviewed parameters')
        frames = complete_frames(path)
        frame = params['expected_frame']
        if spec['provisional'] and frames:
            frame = max(frames)
        if frame not in frames or frame < 1:
            unavailable.append(dict(path=str(path), reason='No complete endpoint diagnostics',
                                    expected_frame=params['expected_frame'], latest_frame=max(frames, default=None)))
            print(f'UNAVAILABLE: {path} (expected frame {params["expected_frame"]})', flush=True)
            continue
        if frames[:frame + 1] != list(range(frame + 1)):
            raise ValueError(f'{path}: missing earlier moment frames')
        before = input_inventory(path, frame)
        print(f'Reading {path} at frame {frame}', flush=True)
        run = cv.load_run(path, frame, 940)
        if args.study == 'alpha-scan':
            first = params['continuous_fdp_start_frame']
            n_oap, _ = cv.modal(path / f'zzim-ion_M0_{first}.gkyl')
            n_fdp, _ = cv.modal(path / f'zzim-ion_M0_{frame}.gkyl')
            mapping, _ = cv.modal(path / 'zzim-geo_corn_mc2nu_pos_deflated.gkyl')
            for data, selected_frame in ((n_oap, first), (n_fdp, frame)):
                if not np.isclose(data.ctx['time'], run['data']['Ei'][0][selected_frame], rtol=0, atol=1e-12):
                    raise ValueError(f'{path}: density and moment frame times disagree')
            summary, profile = cv.density_comparison(n_oap, n_fdp, mapping)
            run['density_comparison'] = dict(summary, oap_frame=first, fdp_frame=frame,
                                             oap_time_s=float(n_oap.ctx['time']),
                                             fdp_time_s=float(n_fdp.ctx['time']))
            run['density_comparison_profile'] = profile
        if input_inventory(path, frame) != before:
            raise ValueError(f'{path}: input changed while reading; rerun to obtain a consistent snapshot')
        if max(v for k, v in run['checks'].items() if k.startswith('max_relative')) > .01:
            raise ValueError(f'{path}: moment reconstruction differs by more than 1%')
        complete = frame == params['expected_frame']
        if complete and not np.isclose(run['end'], params['expected_end_s'], rtol=0, atol=1e-9):
            raise ValueError(f'{path}: saved endpoint disagrees with configured schedule')
        if run['end'] > params['expected_end_s'] + 1e-9:
            raise ValueError(f'{path}: endpoint exceeds configured schedule')
        spec = dict(spec, path=str(path), provisional=not complete)
        if args.study == 'alpha-scan':
            spec['label'] = rf'$\alpha_{{\rm OAP}}={params["alpha_oap"] * 1e3:g}\times10^{{-3}}$'
            spec['short_label'] = f'{params["alpha_oap"] * 1e3:g}'
            label_alpha = float(path.name.removeprefix('alpha-')) if spec['uid'] != 'baseline' else params['alpha_oap']
            if not np.isclose(label_alpha, params['alpha_oap']):
                raise ValueError(f'{path}: directory alpha and input disagree')
        if spec['uid'] == 'fdp_ongoing':
            spec['label'] = rf'+300 $\mu$s + {run["end"] * 1e3:.3g} ms FDP'
            spec['label'] += ' (ongoing)' if not complete else ''
            spec['short_label'] = 'Ongoing' if not complete else 'Extended'
        first = params['continuous_fdp_start_frame']
        if first >= frame:
            raise ValueError(f'{path}: final FDP has not yet been sampled')
        observed_start = run['data']['Ei'][0][first]
        # Boundary frames are written just after the step crossing the boundary.
        tolerance = max(1e-7, params['tau_fdp_extra'] * 1e-5)
        if abs(observed_start - params['continuous_fdp_start_s']) > tolerance:
            raise ValueError(f'{path}: phase boundary frame disagrees with input schedule')
        for stage in params.get('annealing_stages', []):
            for edge in ('start', 'end'):
                observed = float(run['data']['Ei'][0][stage[f'{edge}_frame']])
                # The step entering a new OAP can overshoot by ~0.5 us.
                if abs(observed - stage[f'{edge}_s']) > 1e-6:
                    raise ValueError(f'{path}: annealing stage boundary disagrees with saved frame')
                stage[f'observed_{edge}_s'] = observed
        available_start = max(observed_start, *(s[0][0] for s in run['data'].values()))
        if available_start > run['end'] - max(args.zoom_us, args.window_us) * 1e-6:
            raise ValueError(f'{path}: requested window extends before the sampled FDP interval')
        run.update(spec, parameters=params, frame=frame, complete=complete,
                   phase_start=available_start, files=before,
                   collected_utc=datetime.now(timezone.utc).isoformat(),
                   input_sha256=hashlib.sha256((path / 'sim.c').read_bytes()).hexdigest(),
                   cells=[run['nz'], cv.read_resolution(path, 'vpar', run['nz']),
                          cv.read_resolution(path, 'mu', run['nz'])])
        run['checks']['negative_final_Tperp_samples'] = int(np.sum(run['profile'][1][:, 3] < 0))
        runs.append(run)
    if len(runs) < 2 or not any(r['uid'] == 'baseline' for r in runs):
        raise ValueError('Need the baseline and at least one comparison run')
    if len(runs) > len(COLORS):
        raise ValueError('Extend COLORS and LINESTYLES to give every run a distinct style')
    baseline = next(r for r in runs if r['uid'] == 'baseline')
    if len({r['parameters']['alpha_oap'] for r in runs}) != len(runs) and args.study == 'alpha-scan':
        raise ValueError('Duplicate alpha values')
    for i, run in enumerate(runs):
        run['color'] = COLORS[0] if run['uid'] == 'baseline' else COLORS[1 + sum(r['uid'] != 'baseline' for r in runs[:i])]
        run['linestyle'] = LINESTYLES[0] if run['uid'] == 'baseline' else LINESTYLES[1 + sum(r['uid'] != 'baseline' for r in runs[:i])]
        if run['cells'] != baseline['cells']:
            raise ValueError(f'{run["path"]}: resolution differs from the baseline')
        # Two physical positions per p=1 cell also determine both cell edges,
        # and are retained in snapshots for selective refreshes.
        if not np.allclose(run['profile'][0], baseline['profile'][0], rtol=1e-10, atol=1e-12):
            raise ValueError(f'{run["path"]}: physical geometry mapping differs')
    # Verify that continuation inputs actually point to the compared parent state.
    by_uid = {r['uid']: r for r in runs}
    for uid, parent in [('fdp_300us', 'baseline'), ('fdp_ongoing', 'fdp_300us')]:
        if uid in by_uid:
            run, origin = by_uid[uid], by_uid[parent]
            expected = Path(origin['path']) / f'zzim-ion_{origin["frame"]}.gkyl'
            if Path(run['parameters']['imported_distribution']).resolve() != expected.resolve():
                raise ValueError(f'{uid}: continuation does not import the expected parent endpoint')
    if args.study == 'extra-time-dilation':
        continuation = args.root / '1x-beams-extra-dilated-fdp-continuation'
        frames = complete_frames(continuation)
        unavailable.append(dict(path=str(continuation), reason='Separate continuation; not part of the completed five-cycle comparison',
                                expected_frame=read_parameters(continuation)['expected_frame'], latest_frame=max(frames, default=None)))
    manifest = dict(schema=SCHEMA, study=args.study, root=str(args.root),
                    collected_utc=datetime.now(timezone.utc).isoformat(), unavailable=unavailable,
                    time_convention='saved simulation clock; no inferred physical-time rescaling',
                    runs=[])
    if previous:
        manifest['refresh_runs'] = args.refresh_runs
        manifest['retained_snapshot_utc'] = previous['collected_utc']
    arrays = {}
    for run in runs:
        meta = {key: value for key, value in run.items()
                if key not in ('data', 'profile', 'physical_edges', 'density_comparison_profile')}
        meta['metrics'] = list(run['data'])
        manifest['runs'].append(meta)
        uid = run['uid']
        for key, (t, values) in run['data'].items():
            arrays[f'{uid}__{key}'] = np.vstack((t, values))
        arrays[f'{uid}__profile_z'], arrays[f'{uid}__profile_values'] = run['profile']
        if 'density_comparison_profile' in run:
            arrays[f'{uid}__density_comparison_profile'] = run['density_comparison_profile']
    np.savez_compressed(args.output / 'diagnostic_snapshot.npz', **arrays)
    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest, runs


def restore(args, require_current=True):
    manifest = json.loads((args.output / 'manifest.json').read_text())
    if any(manifest[key] != value for key, value in [('schema', SCHEMA), ('study', args.study), ('root', str(args.root))]):
        raise ValueError('Snapshot configuration differs; rerun without --reuse-data')
    runs = []
    with np.load(args.output / 'diagnostic_snapshot.npz', allow_pickle=False) as data:
        for meta in manifest['runs']:
            uid = meta['uid']
            runs.append(dict(meta, data={key: tuple(data[f'{uid}__{key}']) for key in meta['metrics']},
                             profile=(data[f'{uid}__profile_z'], data[f'{uid}__profile_values'])))
            if 'density_comparison' in meta:
                runs[-1]['density_comparison_profile'] = data[f'{uid}__density_comparison_profile']
    if require_current and args.study == 'alpha-scan' and 7 in args.figures and any(run.get('density_comparison', {}).get('version') != 2 for run in runs):
        raise ValueError('Snapshot lacks relative/absolute OAP/FDP density data; rerun without --reuse-data')
    if require_current and args.study == 'extra-cycles-fdp-length' and not any('annealing_stages' in r['parameters'] for r in runs):
        raise ValueError('Snapshot lacks annealing data; rerun without --reuse-data')
    return manifest, runs


def run_style(run):
    return dict(color=run['color'], ls=run['linestyle'], lw=1.7 if run['uid'] == 'baseline' else 1.35)


def legend(fig, runs, y=.99, ncol=None):
    if ncol is None:
        ncol = 2 if len(runs) > 2 else len(runs)
    fig.legend(handles=[Line2D([], [], **run_style(r), label=r['label']) for r in runs],
               loc='upper center', bbox_to_anchor=(.52, y), ncol=ncol,
               columnspacing=1.6, handlelength=2.6)


def combined(run, metrics, start, end):
    if not start < end:
        raise ValueError('Time window must have positive duration')
    return history_values(run['data'], metrics, start, end)


def statistics(args, runs):
    rows, sampling = [], []
    for run in runs:
        end = run['end']
        for window, start in [('last_15us' if args.window_us == 15 else 'final_window', end - args.window_us * 1e-6),
                              ('last_zoom', end - args.zoom_us * 1e-6), ('last_continuous_fdp', run['phase_start'])]:
            if start < run['phase_start'] - 1e-12:
                raise ValueError(f'{run["uid"]}: {window} begins before the final FDP interval')
            for key, _, _, metrics in ABS_ROWS:
                t, values = combined(run, metrics, start, end)
                mean, low, high = cv.window_stats((t, values), start, end)
                if not np.isfinite([mean, low, high, values[-1]]).all():
                    raise ValueError('Nonfinite statistics')
                rows.append(dict(uid=run['uid'], metric=key, window=window, complete=run['complete'],
                                 start_s=start, end_s=end, final=float(values[-1]), mean=float(mean),
                                 minimum=float(low), maximum=float(high),
                                 range_over_abs_mean_percent=float(100 * (high - low) / abs(mean)) if mean else None))
                for metric, _ in metrics:
                    times = run['data'][metric][0]
                    selected = times[(times >= start) & (times <= end)]
                    lo = max(0, np.searchsorted(times, start, side='right') - 1)
                    hi = min(len(times) - 1, np.searchsorted(times, end))
                    gaps = np.diff(times[lo:hi + 1])
                    sampling.append(dict(uid=run['uid'], window=window, metric=metric,
                                         native_samples_in_window=len(selected),
                                         maximum_bracketing_gap_us=float(gaps.max() * 1e6) if gaps.size else 0,
                                         window_duration_us=(end - start) * 1e6,
                                         interpolation_resolves_window=bool(gaps.size and gaps.max() <= (end - start) / 2)))
    write_csv(args.output / 'window_statistics.csv', rows)
    write_csv(args.output / 'sampling.csv', sampling)
    return rows, sampling


def endpoint_figure(args, runs):
    fig, axs = diagnostic_panels()
    is_alpha_scan = args.study == 'alpha-scan'
    fig.subplots_adjust(left=.13, right=.97, bottom=.12, top=.97 if is_alpha_scan else .88,
                        hspace=.40, wspace=.40)
    rows = ABS_ROWS
    if is_alpha_scan:
        # Alpha tick labels identify the runs; put particle balance across the bottom.
        rows = [ABS_ROWS[i] for i in (0, 1, 2, 4, 3)]
    else:
        legend(fig, runs)
    x = [r['parameters']['alpha_oap'] * 1e3 for r in runs] if args.study == 'alpha-scan' else np.arange(len(runs))
    baseline = next(r for r in runs if r['uid'] == 'baseline')
    exported = []
    for ax, (key, label, scale, metrics) in zip(axs.flat, rows):
        paper_axis(ax)
        values = [sum(w * float(cv.sample(r['data'][k], r['end'])) for k, w in metrics) for r in runs]
        ref = sum(w * float(cv.sample(baseline['data'][k], baseline['end'])) for k, w in metrics)
        if args.study == 'alpha-scan':
            ax.plot(x, np.array(values) * scale, color='#8a9095', lw=.9, zorder=1)
        for xx, run, value in zip(x, runs, values):
            ax.plot(xx, value * scale, marker='D' if run['provisional'] else 'o',
                    color=run['color'], mfc='white' if run['provisional'] else run['color'], ms=6)
            exported.append(dict(uid=run['uid'], metric=key, frame=run['frame'], end_s=run['end'],
                                 complete=run['complete'], final=value, baseline_final=ref,
                                 difference_from_baseline_percent=100 * (value - ref) / abs(ref) if ref else None))
        ax.set_ylabel(SHORT_LABELS.get(key, label))
        ax.set_xticks(x, [r['short_label'] for r in runs])
        ax.tick_params(axis='x', labelsize=12)
        ax.margins(x=.12, y=.15)
        if key == 'Gamma_i_over_source':
            ax.axhline(1, color='#60666b', ls='--', lw=1)
    if args.study == 'alpha-scan':
        for ax in axs[-1:]:
            ax.set_xlabel(r'$\alpha_{\rm OAP}$ [$10^{-3}$]')
    elif args.study == 'extra-time-dilation':
        for ax in axs[-1:]:
            ax.set_xlabel(r'$c_{\omega}$')
    write_csv(args.output / '01_endpoint_diagnostics.csv', exported)
    finish(fig, args, '01_endpoint_diagnostics')


def plot_history(ax, run, metrics, start, end, scale, time_scale):
    times, values = combined(run, metrics, start, end)
    gaps = []
    native_times = []
    for metric, _ in metrics:
        t = run['data'][metric][0]
        lo = max(0, np.searchsorted(t, start, side='right') - 1)
        hi = min(len(t) - 1, np.searchsorted(t, end))
        gaps.extend(np.diff(t[lo:hi + 1]))
        native_times.append(t[(t >= start) & (t <= end)])
    resolved = bool(gaps) and max(gaps) <= (end - start) / 2
    if resolved:
        ax.plot((times - end) * time_scale, values * scale, **run_style(run))
    else:
        # Do not manufacture a short transient from a long interpolation gap.
        # For sums, use only timestamps actually shared by all components.
        times = native_times[0]
        for t in native_times[1:]:
            times = np.intersect1d(times, t)
        values = sum(w * cv.sample(run['data'][k], times) for k, w in metrics) if times.size else np.array([])
        ax.plot((times - end) * time_scale, values * scale, ls='none',
                marker='o', ms=4, color=run['color'], clip_on=False)
    # Show the actual stored times of the two spatially reconstructed histories.
    # A line spanning a sparse interval is interpolation, not additional data.
    if resolved and len(metrics) == 1 and metrics[0][0] in ('Ei', 'n0'):
        t, v = run['data'][metrics[0][0]]
        keep = (t >= start) & (t <= end)
        ax.plot((t[keep] - end) * time_scale, v[keep] * scale, ls='none',
                marker='o', ms=2.6, color=run['color'], clip_on=False)
    return times, values


def histories_figure(args, runs, duration_us, number):
    fig, axs = diagnostic_panels(sharex=True)
    fig.subplots_adjust(left=.13, right=.97, bottom=.09, top=.88, hspace=.30, wspace=.40)
    legend(fig, runs)
    exported = []
    for ax, (key, label, scale, metrics) in zip(axs.flat, ABS_ROWS):
        paper_axis(ax)
        for run in runs:
            end = run['end']
            times, values = plot_history(ax, run, metrics, end - duration_us * 1e-6, end, scale, 1e6)
            exported.extend(dict(uid=run['uid'], metric=key, time_s=float(t),
                                 relative_time_us=float((t - end) * 1e6), value=float(v)) for t, v in zip(times, values))
        ax.set_ylabel(SHORT_LABELS.get(key, label))
        ax.set_xlim(-duration_us, 0)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.margins(y=.12)
        if key == 'Gamma_i_over_source':
            ax.axhline(1, color='#60666b', ls='--', lw=1)
    for ax in axs[-1:]:
        ax.set_xlabel(r'$t-t_{\rm end}$ [$\mu$s]')
    name = f'{number:02d}_last_{duration_us:g}us'
    write_csv(args.output / f'{name}.csv', exported)
    finish(fig, args, name)


def last_fdp_column_figure(args, runs):
    """Compare the last 90 microseconds on one shared time axis."""
    duration_us = 90
    fig, axs = plt.subplots(len(ABS_ROWS), 1, figsize=(8, 3 * len(ABS_ROWS) - 1), sharex=True)
    fig.subplots_adjust(left=.24, right=.98, bottom=.065, top=.90, hspace=.28)
    fig.legend(handles=[Line2D([], [], **run_style(r), label=r['label']) for r in runs],
               loc='upper center', bbox_to_anchor=(.55, .995), ncol=1, handlelength=2.6)
    exported = []
    for ax, (key, label, scale, metrics) in zip(axs, ABS_ROWS):
        paper_axis(ax)
        boundaries = []
        for run in runs:
            end = run['end']
            # The first saved FDP sample can be just after the nominal boundary.
            start = max(end - duration_us * 1e-6, run['phase_start'])
            t, values = plot_history(ax, run, metrics, start, end, scale, 1e6)
            exported.extend(dict(uid=run['uid'], metric=key, time_s=float(tt),
                                 relative_time_s=float(tt - end), value=float(v)) for tt, v in zip(t, values))
            boundary_us = (run['parameters']['final_phase_start_s'] - end) * 1e6
            if -duration_us < boundary_us < 0 and not any(np.isclose(boundary_us, b) for b in boundaries):
                ax.axvline(boundary_us, color='#90979b', lw=.8, ls=':')
                boundaries.append(boundary_us)
        ax.set_ylabel(SHORT_LABELS.get(key, label))
        ax.set_xlim(-duration_us, 0)
        ax.set_xticks([-90, -60, -30, 0])
        if key == 'Gamma_i_over_source':
            ax.axhline(1, color='#60666b', ls='--', lw=1)
    axs[-1].set_xlabel(r'$t-t_{\rm end}$ [$\mu$s]')
    write_csv(args.output / '04_last_full_dynamics.csv', exported)
    finish(fig, args, '04_last_full_dynamics')


def full_fdp_figure(args, runs):
    if args.study in ('alpha-scan', 'extra-cycles-fdp-length'):
        return last_fdp_column_figure(args, runs)
    # Each column uses its own true time extent; no fractional-phase rescaling.
    n = len(runs)
    fig, axs = plt.subplots(len(ABS_ROWS), n, figsize=(4 * n + 1.5, 3 * len(ABS_ROWS) - 1),
                            sharey='row', squeeze=False)
    fig.subplots_adjust(left=.12, right=.98, bottom=.065, top=.94, hspace=.28, wspace=.18)
    exported = []
    for col, run in enumerate(runs):
        end, start = run['end'], run['phase_start']
        duration = end - start
        time_scale, unit = (1e3, 'ms') if duration >= .002 else (1e6, r'\mu\mathrm{s}')
        center = (axs[0, col].get_position().x0 + axs[0, col].get_position().x1) / 2
        fig.legend(handles=[Line2D([], [], **run_style(run), label=run['label'])],
                   loc='upper center', bbox_to_anchor=(center, .99), handlelength=1.7)
        for row, (key, label, scale, metrics) in enumerate(ABS_ROWS):
            ax = axs[row, col]
            paper_axis(ax)
            t, values = plot_history(ax, run, metrics, start, end, scale, time_scale)
            exported.extend(dict(uid=run['uid'], metric=key, time_s=float(tt),
                                 relative_time_s=float(tt - end), value=float(v)) for tt, v in zip(t, values))
            boundary = run['parameters']['final_phase_start_s']
            if start < boundary < end:
                ax.axvline((boundary - end) * time_scale, color='#90979b', lw=.8, ls=':')
            if key == 'Gamma_i_over_source':
                ax.axhline(1, color='#60666b', ls='--', lw=1)
            if col == 0:
                ax.set_ylabel(SHORT_LABELS.get(key, label))
            ax.set_xlim(-duration * time_scale, 0)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.tick_params(axis='x', labelbottom=row == len(ABS_ROWS) - 1)
            if row == len(ABS_ROWS) - 1:
                ax.set_xlabel(rf'$t-t_{{\rm end}}$ [${unit}$]')
    write_csv(args.output / '04_last_full_dynamics.csv', exported)
    finish(fig, args, '04_last_full_dynamics')


def profiles_figure(args, runs):
    fig = plt.figure(figsize=(10, 13))
    grid = fig.add_gridspec(4, 2, left=.13, right=.90, bottom=.08, top=.89, hspace=.28, wspace=0)
    labels = (r'$n_i$ [m$^{-3}$]', r'$u_\parallel$ [$10^6$ m s$^{-1}$]',
              r'$T_\parallel$ [keV]', r'$T_\perp$ [keV]')
    scales = (1, 1e-6, 1, 1)
    legend(fig, runs)
    exported = []
    extent = max(np.max(abs(r['profile'][0])) for r in runs) * 1.015
    for row, (label, scale) in enumerate(zip(labels, scales)):
        left, right = fig.add_subplot(grid[row, 0]), fig.add_subplot(grid[row, 1])
        for side, ax in enumerate((left, right)):
            paper_axis(ax)
            if side:
                ax.set_yscale('log', nonpositive='mask')
                ax.yaxis.tick_right()
                ax.spines['right'].set_visible(True)
                ax.spines['left'].set_visible(False)
                ax.set_xlim(0, extent)
                ax.set_xticks([1, 2])
            else:
                ax.spines['right'].set_visible(True)
                ax.spines['right'].set_color('#c3c8cc')
                ax.set_xlim(-extent, 0)
                ax.set_xticks([-2, -1, 0])
                ax.set_ylabel(label)
            for run in runs:
                z, values = run['profile']
                y = values[:, row] * scale
                keep = (z >= 0) & (y > 0) if side else z < 0
                zz = np.column_stack((z.reshape(-1, 2), np.full(len(z) // 2, np.nan))).ravel()
                yy = np.column_stack((np.where(keep, y, np.nan).reshape(-1, 2), np.full(len(z) // 2, np.nan))).ravel()
                ax.plot(zz, yy, **run_style(run))
            ax.axvline(runs[0]['parameters']['Z_m'] * (1 if side else -1), color='#888888', ls='--', lw=1)
            ax.tick_params(axis='x', labelbottom=row == 3)
            ax.margins(y=.10)
        # Preserve the entire linear range, including negative temperatures.
        if row == 0:
            left.set_ylim(bottom=min(0, left.get_ylim()[0]))
        if row in (2, 3):
            low, high = left.get_ylim()
            if low < 0:
                left.axhspan(low, 0, color='#f5e5e3', zorder=0)
                left.axhline(0, color='#8e3030', ls=':', lw=.8)
                left.set_ylim(low, high)
        if row == 3:
            right.set_xlabel(r'$z$ [m]')
            right.xaxis.set_label_coords(0, -.18)
    for run in runs:
        z, values = run['profile']
        exported.extend(dict(uid=run['uid'], frame=run['frame'], end_s=run['end'],
                             z_m=float(zz), density_m3=float(v[0]), upar_ms=float(v[1]),
                             Tpar_keV=float(v[2]), Tperp_keV=float(v[3])) for zz, v in zip(z, values))
    write_csv(args.output / '05_endpoint_profiles.csv', exported)
    finish(fig, args, '05_endpoint_profiles')


def balance_figure(args, runs):
    fig, ax = plt.subplots(figsize=(9, 5))
    fig.subplots_adjust(left=.13, right=.87, bottom=.17, top=.74 if len(runs) > 4 else .80)
    legend(fig, runs)
    paper_axis(ax)
    exported = []
    for run in runs:
        end = run['end']
        times, values = plot_history(ax, run, [('particle_balance', 1)],
                                    end - args.window_us * 1e-6, end, 1, 1e6)
        exported.extend(dict(uid=run['uid'], relative_time_us=float((t-end)*1e6), ratio=float(v),
                             deviation_from_unity_percent=float(100*(v-1))) for t,v in zip(times,values))
    ax.axhline(1, color='#60666b', ls='--', lw=1)
    ax.set_xlim(-args.window_us, 0)
    ax.set_ylabel(r'$\Gamma_i\,/\!\int S_i\,d^3x\,d^3v$')
    ax.set_xlabel(r'$t-t_{\rm end}$ [$\mu$s]')
    deviation = ax.secondary_yaxis('right', functions=(lambda x: 100*(x-1), lambda x: 1+x/100))
    deviation.set_ylabel(r'Deviation from unity [\%]')
    deviation.spines['right'].set_visible(True)
    write_csv(args.output / '06_particle_balance.csv', exported)
    finish(fig, args, '06_particle_balance')


def annealing_figure(args, runs):
    """Resolve the ramp and final FDP on separate simulation-clock axes."""
    run = next(r for r in runs if r['uid'] == 'annealing')
    baseline = next(r for r in runs if r['uid'] == 'baseline')
    stages = run['parameters']['annealing_stages']
    rows = [ABS_ROWS[i] for i in (0, 1, 3)]
    fig, axs = plt.subplots(len(rows), len(stages), figsize=(16, 10), sharey='row', squeeze=False)
    fig.subplots_adjust(left=.10, right=.98, bottom=.11, top=.85, hspace=.28, wspace=.15)
    fig.legend(handles=[Line2D([], [], **run_style(run), label='Annealing'),
                        Line2D([], [], color=baseline['color'], ls='--', lw=1,
                               label='Baseline final value')],
               loc='upper center', bbox_to_anchor=(.54, .99), ncol=2)
    exported, schedule = [], []
    for col, stage in enumerate(stages):
        start, end = stage['observed_start_s'], stage['observed_end_s']
        time_scale, unit = (1e3, 'ms') if stage['duration_s'] >= .001 else (1e6, r'\mu\mathrm{s}')
        title = rf'{stage["kind"]}: $\alpha={stage["alpha"]:g}$'
        axs[0, col].set_title(title + '\n' + rf'{stage["duration_s"] * time_scale:g} ${unit}$', fontsize=16)
        schedule.append(dict(stage=col + 1, **stage))
        for row, (key, label, scale, metrics) in enumerate(rows):
            ax = axs[row, col]
            paper_axis(ax)
            t, values = plot_history(ax, run, metrics, start, end, scale, time_scale)
            exported.extend(dict(uid=run['uid'], stage=col + 1, kind=stage['kind'], alpha=stage['alpha'],
                                 metric=key, time_s=float(tt), relative_to_stage_end_s=float(tt - end),
                                 value=float(v)) for tt, v in zip(t, values))
            ref = sum(w * float(cv.sample(baseline['data'][k], baseline['end'])) for k, w in metrics)
            ax.axhline(ref * scale, color=baseline['color'], ls='--', lw=1)
            if key == 'Gamma_i_over_source':
                ax.axhline(1, color='#60666b', ls=':', lw=.8)
            if col == 0:
                ax.set_ylabel(SHORT_LABELS.get(key, label))
            ax.set_xlim((start - end) * time_scale, 0)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.tick_params(axis='x', labelbottom=row == len(rows) - 1)
            ax.margins(y=.12)
            if row == len(rows) - 1:
                ax.set_xlabel(rf'$t-t_{{\rm stage,end}}$ [${unit}$]', fontsize=17)
    write_csv(args.output / '07_annealing_stages.csv', exported)
    write_csv(args.output / '07_annealing_schedule.csv', schedule)
    finish(fig, args, '07_annealing_stages')


def density_comparison_figure(args, runs):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.3))
    fig.subplots_adjust(left=.10, right=.97, bottom=.17, top=.86, wspace=.34)
    legend(fig, runs, ncol=len(runs))
    frames = {(r['density_comparison']['oap_frame'], r['density_comparison']['fdp_frame']) for r in runs}
    if len(frames) != 1:
        raise ValueError('Alpha density comparison requires matching phase frames')
    profiles, summary = [], []
    for ax in axs.flat:
        paper_axis(ax)
    for run in runs:
        values = run['density_comparison_profile']
        # NaNs preserve independent one-sided DG traces at cell interfaces.
        def broken(v):
            return np.column_stack((v.reshape(-1, 2), np.full(len(v) // 2, np.nan))).ravel()
        axs[0].plot(broken(values[:, 0]), broken(values[:, 3] * 1e-17), **run_style(run))
        common = dict(uid=run['uid'], alpha_oap=run['parameters']['alpha_oap'])
        summary.append(dict(common, **run['density_comparison']))
        profiles.extend(dict(common, cell=i // 2, trace='left' if i % 2 == 0 else 'right',
                             z_m=float(v[0]), n_oap_m3=float(v[1]), n_fdp_m3=float(v[2]),
                             delta_n_m3=float(v[3]), relative_delta_n=float(v[4]))
                        for i, v in enumerate(values))
    axs[0].set_xlabel(r'$z$ [m]')
    axs[0].set_xlim(summary[0]['z_min_m'], summary[0]['z_max_m'])
    for sign in (-1, 1):
        axs[0].axvline(sign * runs[0]['parameters']['Z_m'], color='#888888', ls='--', lw=.8)
    axs[0].set_ylabel(r'$\Delta n_i$ [$10^{17}$ m$^{-3}$]')
    axs[0].yaxis.set_major_locator(MaxNLocator(nbins=5))
    axs[0].axhline(0, color='#60666b', ls=':', lw=1)
    ax = axs[1]
    x = [row['alpha_oap'] * 1e3 for row in summary]
    y = [row['absolute_integral_m2'] * 1e-16 for row in summary]
    ax.plot(x, y, color='#8a9095', lw=.9, zorder=1)
    for xx, yy, run in zip(x, y, runs):
        ax.plot(xx, yy, marker='o', color=run['color'], ms=6)
    ax.axhline(0, color='#60666b', ls=':', lw=1)
    ax.set_xticks(x, [r['short_label'] for r in runs])
    ax.set_xlabel(r'$\alpha_{\rm OAP}$ [$10^{-3}$]')
    ax.set_ylabel(r'$\int |\Delta n|\,dz$ [$10^{16}$ m$^{-2}$]')
    ax.margins(x=.12, y=.15)
    ax.set_ylim(bottom=0)
    write_csv(args.output / '07_oap_fdp_density_comparison.csv', summary)
    write_csv(args.output / '07_oap_fdp_density_profiles.csv', profiles)
    finish(fig, args, '07_oap_fdp_density_comparison')


def report(args, manifest, runs, sampling):
    fdp_caption = [
        '4. The entire available last continuous FDP interval, in separate columns with shared vertical scales.',
        '   Each column uses its true duration and stated units. No time normalization or stitching across restart clocks is applied.',
    ]
    if args.study == 'alpha-scan':
        fdp_caption = [
            f'4. The full final 90-microsecond continuous FDP interval for all {len(runs)} alpha cases, overlaid in a single column of five diagnostic panels.',
            '   All panels share the same endpoint-aligned time axis; each curve starts at its first saved FDP sample.',
            '   Dots mark saved energy and central-density samples. The dotted vertical line marks the start of the configured extra FDP.',
        ]
    if args.study == 'extra-cycles-fdp-length':
        fdp_caption = [
            f'4. The last up to 90 microseconds of all {len(runs)} runs, overlaid in a single column of five diagnostic panels.',
            '   All panels share the same endpoint-aligned time axis. Sparse energy and central-density histories show native samples only.',
            '   Annealing has only 75 microseconds of continuous final FDP; its curve begins there, excluding the preceding OAP.',
            '   The historical output filename 04_last_full_dynamics is retained; this view is limited to the last 90 microseconds.',
        ]
    lines = [f'PUBLICATION FIGURES: {args.study}',
             f'Snapshot assembled: {manifest["collected_utc"]}', '',
             'DIAGNOSTICS AND TIME CONVENTION',
             'Ni is the particle count from ion_integrated_moms component 0.',
             'Ei=(mi/2) integral(M2_i d^3x), reconstructed with the exact p=1 modal product and geometric Jacobian.',
             'Ei excludes electrostatic potential energy. Number and Hamiltonian reconstructions are independently checked.',
             'n0 averages the two one-sided DG traces at physical z=0 when z=0 is a cell interface.',
             'Gamma_i sums outward ion particle loss through both walls. The denominator is ion_source_integrated_moms component 0.',
             'qi is the sum of the saved ion Hamiltonian energy fluxes integrated over the lower and upper walls, without division by wall area.',
             'The snapshot reader stores these integrated wall fluxes under Pi_lo and Pi_hi; the figures show their sum only as qi.',
             'Units follow the author-specified SI normalization: Ni is a count, Ei J, n0 m^-3, qi W (displayed in kW). Plotted CSVs use qi_W.',
             'All clocks are the saved simulation clock. A time-dilation/CFL factor is not a global physical-time conversion.',
             'Each curve is aligned to its own selected endpoint. Different phase schedules are not synchronized physical transients.',
             'The full-FDP plot starts at the saved boundary after the last OAP, including adjacent FDP phases with unchanged settings.',
             'A dotted vertical line marks the start of the configured extra FDP. For a continuation, only that continuation segment is shown.',
             'The first boundary sample can be slightly after the nominal boundary; no diagnostic is extrapolated.',
             '', 'FIGURE CAPTIONS',
             '1. Endpoint ion diagnostics. Markers are individual saved endpoints, not temporal means. Open diamonds denote an ongoing run.',
             '   Alpha points are joined only to guide the eye. Other comparisons are categorical. CSV percent differences use the baseline endpoint; they are not numerical errors.',
             f'2. Final {args.window_us:g} microseconds of the five ion diagnostics, aligned to the selected endpoint.',
             f'3. Final {args.zoom_us:g} microseconds of the same diagnostics.',
             '   Curves retain native integrated-diagnostic samples. Energy and central density use piecewise-linear interpolation; dots show their actual saved times.',
             '   If any bracketing sample gap exceeds half the displayed window, only native samples are plotted, without an interpolated curve.',
             '   For combined diagnostics in such a window, only timestamps shared by all contributing components are plotted.',
             '   Window statistics still use piecewise-linear interpolation and must be read together with sampling.csv; a sparse window mean does not establish relaxation.',
             *fdp_caption,
             '5. Endpoint ion density, parallel velocity, parallel temperature, and perpendicular temperature.',
             '   Negative z uses linear axes; positive z uses logarithmic axes. Nonpositive right-half samples are masked, never replaced by absolute values.',
             '   DG cell interfaces remain gaps. Linear panels retain all samples, including negative temperatures (red background). Dashed lines mark mirror throats.',
             f'6. Source-loss ratio over the final {args.window_us:g} microseconds; the right axis gives 100*(ratio-1). Unity denotes source-loss balance.',
             '']
    if args.study == 'alpha-scan' and all(r.get('density_comparison', {}).get('version') == 2 for r in runs):
        lines += ['7. Ion density difference n_OAP-n_FDP versus z in units of 10^17 m^-3 (left) and integral |n_OAP-n_FDP| dz versus alpha (right).',
                  '   Frames are derived from the phase schedule (45 and 75 for the standard runs; 90 us apart). Saved boundary times are slightly after the nominal OAP end.',
                  '   Density is ion_M0. Postgkyl performs native modal subtraction for the plotted density difference. The additional relative profile in the CSV uses weak division by n_FDP; the denominator must be positive throughout each cell, with no density floor imposed.',
                  '   For the absolute integral, exact means of |delta n| are computed from each p=1 polynomial, splitting at any within-cell zero crossing. These modal means are weakly multiplied by dz/dz_comp from the local DG derivative of geo_corn_mc2nu_pos_deflated and integrated through Postgkyl.',
                  '   The native integrate(op="abs") shortcut is not used: its cell-mean absolute value can undercount within-cell sign changes.',
                  '   The full physical z domain is used. This is a line integral in m^-2, not the volume-integrated particle-number difference; no flux-tube area/Jacobian is included.',
                  '   The density difference profile retains its sign. The absolute integral accumulates all local density differences without cancellation; it is not the absolute value of the signed integral.',
                  '   Both one-sided cell-end traces are shown in the profiles; they are not used for numerical quadrature. Dashed vertical lines mark the mirror throats.',
                  '   The summary CSV includes both signed and absolute integrals and their percentages of the FDP line density. The profile CSV adds the weak relative difference. Alpha-point connections only guide the eye.', '']
    if args.study == 'extra-cycles-fdp-length':
        lines += ['7. Annealing-stage evolution of ion number, kinetic energy, and source-loss ratio after four regular OAP/FDP cycles.',
                  '   Columns show OAP alpha=0.001 for 100 ms, alpha=0.01 for 10 ms, alpha=0.1 for 1 ms, then the 75-us final FDP.',
                  '   Each column uses its own simulation-clock axis, relative to its saved stage-end frame. Shared row scales compare the stages; dots mark stored energy samples.',
                  '   Dashed purple lines give the baseline final value for reference; the gray dotted line marks source-loss balance.',
                  '   07_annealing_schedule.csv records nominal and observed boundaries (frames 40, 45, 50, 55, 80). Boundary frames are written just after crossing a phase boundary.',
                  '   The source allocates one unused zero-duration phase; saved output confirms the nonzero stages and frame 80 at 0.511135 s. Stale input comments are not used as phase definitions.',
                  '   The archived 1x-beams-anneal/256 data are excluded from this fixed-resolution comparison; the main annealing run matches the baseline 400 x 64 x 32 grid.', '']
    lines += ['RUNS AND CHECKS']
    for run in runs:
        p = run['parameters']
        lines += [f'{run["uid"]}: {run["path"]}',
                  f'  raw-data collection: {run.get("collected_utc", manifest["collected_utc"])}',
                  f'  frame={run["frame"]}/{p["expected_frame"]}; endpoint={run["end"]:.12g} s; completed={run["complete"]}; changeset={run["changeset"]}',
                  f'  cycles={p["num_cycles"]}; alpha_oap={p["alpha_oap"]:g}; c_omega={p["cfl_factor_times_omega_max"]:g}; tau_fdp={p["tau_fdp"]:g} s; tau_extra={p["tau_fdp_extra"]:g} s',
                  f'  observed full-FDP interval: {run["phase_start"]:.12g} to {run["end"]:.12g} s ({(run["end"]-run["phase_start"])*1e6:.9g} us)',
                  f'  reconstruction / temperature checks: {run["checks"]}']
        if 'density_comparison' in run:
            lines.append(f'  OAP/FDP density comparison: {run["density_comparison"]}')
        if 'annealing_stages' in p:
            lines.append(f'  annealing stages: {p["annealing_stages"]}')
    if len({r['changeset'] for r in runs}) > 1:
        lines += ['Saved solver revisions differ across this comparison; per-run revisions are recorded above. Input and diagnostic checks do not establish code-revision equivalence.']
    if args.study == 'extra-time-dilation':
        lines += ['', 'The extra-dilation input changes c_omega from 0.1 to 0.001 AND tau_fdp from 15 us to 5 ms and tau_extra from 75 us to 25 ms.',
                  'It is therefore a joint time-dilation / duration sensitivity comparison, not an isolated scan of c_omega.',
                  'c_omega denotes the input cfl_factor_times_omega_max. The OAP collisionless multiplier remains at 0.1.',
                  'The final continuous FDP intervals are 90 us and 30 ms, respectively.']
    lines += ['', 'SAMPLING LIMITATIONS (see sampling.csv for all diagnostics/windows)']
    for row in sampling:
        if not row['interpolation_resolves_window']:
            lines.append(f'{row["uid"]}/{row["window"]}/{row["metric"]}: {row["native_samples_in_window"]} stored samples in window; maximum bracketing gap {row["maximum_bracketing_gap_us"]:.6g} us.')
    lines += ['', 'UNAVAILABLE / EXCLUDED']
    lines += [json.dumps(row) for row in manifest['unavailable']] or ['None.']
    if manifest.get('refresh_runs'):
        lines += ['', 'SELECTIVE REFRESH',
                  f'Only {", ".join(manifest["refresh_runs"])} were reread when assembling this snapshot.',
                  'Other runs retain their recorded arrays, endpoints, collection times, and source inventories.',
                  'Use --reuse-data to reproduce this selection without reading newer simulation outputs.']
    lines += ['', 'REPRODUCTION', 'python publication_figures.py',
              'python publication_figures.py --reuse-data  # explicitly replot the recorded snapshot',
              'python publication_figures.py --refresh-runs annealing  # cycle/FDP study: refresh only annealing in an existing snapshot',
              'Use --root to select the simulation parent directory, --output to change the destination, and --figures to select figures.',
              'Default: raw diagnostics; PDF and SVG vectors plus 300 dpi PNG. Snapshots, plotted CSVs, window statistics, and input provenance accompany the figures.']
    (args.output / 'captions_and_methods.txt').write_text('\n'.join(lines) + '\n')


def main(study, script_dir):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=script_dir.parents[1])
    parser.add_argument('--output', type=Path, default=script_dir / 'publication-figures')
    parser.add_argument('--window-us', type=float, default=15)
    parser.add_argument('--zoom-us', type=float, default=60)
    parser.add_argument('--dpi', type=int, default=300)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--reuse-data', action='store_true')
    mode.add_argument('--refresh-runs', nargs='+', metavar='UID',
                      help='Read only these run IDs; preserve other runs from the output snapshot')
    figure_numbers = list(range(1, 8 if study in ('alpha-scan', 'extra-cycles-fdp-length') else 7))
    parser.add_argument('--figures', type=int, nargs='+', choices=figure_numbers, default=figure_numbers)
    args = parser.parse_args()
    if not np.isfinite([args.window_us, args.zoom_us, args.dpi]).all() or min(args.window_us, args.zoom_us, args.dpi) <= 0:
        parser.error('Window lengths and dpi must be finite and positive')
    args.root, args.output, args.study = args.root.resolve(), args.output.resolve(), study
    args.output.mkdir(parents=True, exist_ok=True)
    manifest, runs = restore(args) if args.reuse_data else collect(args)
    _, sampling = statistics(args, runs)
    style()
    if 1 in args.figures:
        endpoint_figure(args, runs)
    if 2 in args.figures:
        histories_figure(args, runs, args.window_us, 2)
    if 3 in args.figures:
        histories_figure(args, runs, args.zoom_us, 3)
    if 4 in args.figures:
        full_fdp_figure(args, runs)
    if 5 in args.figures:
        profiles_figure(args, runs)
    if 6 in args.figures:
        balance_figure(args, runs)
    if 7 in args.figures:
        if study == 'extra-cycles-fdp-length':
            annealing_figure(args, runs)
        else:
            density_comparison_figure(args, runs)
    report(args, manifest, runs, sampling)
    print(f'Completed figures {sorted(set(args.figures))} for {len(runs)} runs in {args.output}', flush=True)
