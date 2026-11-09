"""Shared, checked p=1 Gkeyll diagnostics for the publication analyses.

Extracted from the resolution study's convergence.py. Only moment and geometry
files are read; full distribution functions are not loaded.
"""
from __future__ import annotations
import re
import warnings
from pathlib import Path
import numpy as np
import postgkyl as pg

AXES = {'z': ('Nz', r'N_z'), 'vpar': ('Nvpar', r'N_{v_\parallel}'), 'mu': ('Nmu', r'N_\mu')}

def read(path):
    if not path.is_file():
        raise FileNotFoundError(2, 'Missing diagnostic', str(path))
    d = pg.GData(str(path))
    if not np.isfinite(d.get_values()).all():
        raise ValueError(f'Non-finite values in {path}')
    return d


def series(path, components):
    """Read append-only diagnostics, letting later restarts replace old tails."""
    d = read(path)
    t, v = np.asarray(d.get_grid()[0]), d.get_values()
    if t.ndim != 1 or not len(t) or not np.isfinite(t).all():
        raise ValueError(f'Expected nonempty, finite diagnostic times in {path}')
    if v.ndim != 2 or v.shape != (len(t), components):
        raise ValueError(f'Unexpected diagnostic shape in {path}: {v.shape}')
    # A restart appends samples from an earlier time without removing the old
    # tail. Keep a sample only if it precedes every later-written sample. This
    # discards superseded tails even when restart timestamps differ slightly,
    # and keeps the last copy of an exactly duplicated phase-boundary time.
    # Sorting would incorrectly interleave old and restarted trajectories.
    future_min = np.minimum.accumulate(t[::-1])[::-1]
    keep = np.r_[t[:-1] < future_min[1:], True]
    restarts = np.count_nonzero(np.diff(t) < 0)
    if restarts:
        warnings.warn(f'{path}: reconciled {restarts} appended restart(s); '
                      f'discarded {np.count_nonzero(~keep)} superseded samples',
                      stacklevel=2)
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


def absolute_p1_cell_means(data):
    """Modal constant field with each cell's exact mean of |data|.

    For f=A+B*xi, a root r=-A/B inside [-1,1] gives mean(|f|)
    = |B|*(1+r*r)/2. Keeping this mean suffices for integration against
    the cellwise constant derivative of the p=1 coordinate map.
    """
    coefficients = data.get_values()
    a, b = coefficients[:, 0] / np.sqrt(2), coefficients[:, 1] * np.sqrt(1.5)
    mean = np.abs(a)
    crossing = np.abs(a) < np.abs(b)
    root = -a[crossing] / b[crossing]
    mean[crossing] = np.abs(b[crossing]) * (1 + root**2) / 2
    values = np.column_stack((mean * np.sqrt(2), np.zeros_like(mean)))
    result = pg.GData(ctx=dict(data.ctx))
    result.push(data.get_grid(), pg.gpython.GkylArray.from_numpy(values))
    return result


def density_comparison(n_oap, n_fdp, mapping):
    """Weak relative density change and absolute/signed physical z integrals.

    Use Postgkyl weak subtraction/division and native modal integration.
    Absolute values use exact p=1 cell means, including within-cell roots.
    The coordinate map supplies dz/dz_comp, not the flux-tube volume Jacobian.
    Returned profiles use both one-sided cell endpoints only for plotting.
    """
    for data in (n_oap, n_fdp, mapping):
        if (data.backend != 'gkyl' or data.ctx.get('value_form', 'modal') != 'modal'
                or data.ctx.get('basis_type') != 'serendipity'
                or data.ctx.get('poly_order') != 1
                or data.get_values().shape != (int(data.ctx['cells'][0]), 2)):
            raise ValueError('Density comparison requires scalar native 1D p=1 modal data')
        if not np.isfinite(data.get_values()).all():
            raise ValueError('Nonfinite density comparison input')
    edges = evaluate(mapping.get_values(), [-1, 1])
    if (np.any(np.diff(edges, axis=1) <= 0)
            or not np.allclose(edges[:-1, 1], edges[1:, 0], rtol=0, atol=1e-10)):
        raise ValueError('Physical z mapping must be increasing and continuous')
    dz_dcomp = mapping.differentiate(direction=0)
    difference = n_oap - n_fdp
    if np.any(evaluate(n_fdp.get_values(), [-1, 1]) <= 0):
        raise ValueError('Relative density comparison requires positive FDP density throughout each cell')
    relative = difference / n_fdp
    if not np.isfinite(relative.get_values()).all():
        raise ValueError('Nonfinite weak relative density difference')
    signed = float((difference * dz_dcomp).integrate())
    absolute = float((absolute_p1_cell_means(difference) * dz_dcomp).integrate())
    oap = float((n_oap * dz_dcomp).integrate())
    fdp = float((n_fdp * dz_dcomp).integrate())
    # Independent exact p=1 check using physical cell widths and mean modes.
    expected = float(np.sum(np.diff(edges, axis=1)[:, 0]
                            * difference.get_values()[:, 0] / np.sqrt(2)))
    scale = max(abs(oap), abs(fdp), 1.)
    if (not np.isfinite([signed, absolute, oap, fdp, expected]).all() or fdp <= 0
            or absolute < abs(signed) - 1e-12 * scale
            or abs(signed - expected) > 1e-12 * scale
            or abs(signed - (oap - fdp)) > 1e-12 * scale):
        raise ValueError('Physical density integral failed its independent checks')
    summary = dict(version=2, absolute_integral_m2=absolute,
                   absolute_relative_to_fdp_percent=100 * absolute / fdp,
                   signed_integral_m2=signed, oap_integral_m2=oap,
                   fdp_integral_m2=fdp, relative_to_fdp_percent=100 * signed / fdp,
                   z_min_m=float(edges[0, 0]), z_max_m=float(edges[-1, 1]),
                   integration_check_relative=abs(signed - expected) / scale)
    profiles = np.column_stack([edges.ravel(),
                               *(evaluate(d.get_values(), [-1, 1]).ravel()
                                 for d in (n_oap, n_fdp, difference, relative))])
    return summary, profiles


def input_signature(path, allow_mapping=False, axis='z'):
    # Only the scanned resolution may differ. Refuse accidental mixed physics.
    s = path.read_text()
    s = re.sub(r'//[^\n]*|/\*.*?\*/', '', s, flags=re.S)
    # Runs copied between machines retain the same equilibrium file basename.
    # Match the compatibility check already used by the scan-local scripts.
    s = re.sub(r'(\.filename_psi\s*=\s*)"([^"]+)"',
               lambda m: f'{m[1]}"{Path(m[2]).name}"', s)
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
    source_t, source, _ = series(f('ion_source_integrated_moms'), 3)
    if np.any(source[:, 0] <= 0):
        raise ValueError(f'{path}: particle source must be positive for flux normalization')
    out['source'] = (source_t, source[:, 0])
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
        # Compare integrated rates in the same transverse normalization.
        out[f'particle_balance_{side}'] = (t, 2 * b[:, 0] / sample(out['source'], t))
        out[f'Pi_{side}'] = (t, b[:, 2])
        out[f'Ptot_{side}'] = (t, b[:, 2] + 2 * temp * b[:, 0])
        out[f'gamma_{side}'] = (t, b[:, 0] / area)
        out[f'qi_{side}'] = (t, b[:, 2] / area)
        out[f'qe_{side}'] = (t, 2 * temp * b[:, 0] / area)
        out[f'qtot_{side}'] = (t, (b[:, 2] + 2 * temp * b[:, 0]) / area)
    for key in ['loss', 'Pi', 'Ptot']:
        t, v = out[f'{key}_lo']
        out[key] = (t, v + sample(out[f'{key}_hi'], t))
    t, loss = out['loss']
    out['particle_balance'] = (t, loss / sample(out['source'], t))

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


