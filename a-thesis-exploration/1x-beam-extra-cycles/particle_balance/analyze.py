"""Audit saved FDP particle balances without modifying simulation outputs.

Run: MPLCONFIGDIR=/tmp/mpl-particle python particle_balance/analyze.py

The output metadata identifies two different filter implementations:
844defc6166c: f_t = a*(L+C+S-g), fbar_t = g
5ca5b5f29a02: f_t = a*(L+C+S-g), fbar_t = a*g
where g = rate*(f-fbar), and a is the cellwise time-rate multiplier.
Both satisfy integral(f_t/a + g) = S_tot - outward_boundary_flux
in FDP (collisionless scale = 1, no collisionless-only multiplier).
This is an instantaneous weighted rate, not d/dt integral(f/a) when a varies.

fdot phase-space files contain the following step's stage-1 RHS, whereas
integrated fdot and boundary flux are averages over the preceding RK step.
Small corrected residuals therefore include a sampling-time mismatch.
The saved multiplier is from the preceding step; this is another small mismatch.
Only late FDP frames are analyzed; OAP needs its own modified balance.
"""

import csv
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/mpl-particle')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import postgkyl as pg

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
FILTER_MULTIPLIED = {'844defc6166c': False, '5ca5b5f29a02': True}


def read(path, quantity):
    return pg.GData(str(path / f'zzim-ion_{quantity}.gkyl'))


def audit(path, frames, label):
    series = {}
    for name in ['source_integrated_moms', 'integrated_moms',
                 'fdot_integrated_moms', 'bflux_xlower_integrated_HamiltonianMoments',
                 'bflux_xupper_integrated_HamiltonianMoments', 'fbardot_integrated_moms']:
        if (path / f'zzim-ion_{name}.gkyl').exists():
            d = read(path, name)
            series[name] = (np.asarray(d.get_grid()[0]), d.get_values()[:, 0])
    t, source = series['source_integrated_moms']
    for name, (tx, _) in series.items():
        assert np.array_equal(t, tx), f'Unaligned diagnostic: {name}'
    flux = (series['bflux_xlower_integrated_HamiltonianMoments'][1]
            + series['bflux_xupper_integrated_HamiltonianMoments'][1])
    ndot = series['fdot_integrated_moms'][1]
    rows = []
    for frame in frames:
        fdata = read(path, str(frame))
        f = fdata.get_values()[..., 0].copy()
        tm = fdata.ctx['time']
        version = fdata.ctx['changeset']
        bar_multiplied = FILTER_MULTIPLIED[version]
        assert fdata.ctx['basis_type'] == 'gkhybrid'
        assert fdata.ctx['poly_order'] == 1 and f.ndim == 3
        # Exact M0 modal integral: J and velocity-map Jacobian are already in f.
        # gyrokinetic_int_M0_1x2v_ser_p1: sqrt(8)*(pi/4)*prod(dx)/mass*f[0].
        factor = (2*np.pi / fdata.ctx['mass'] / np.sqrt(8)
                  * np.prod([x[1]-x[0] for x in fdata.get_grid()]))
        i = int(np.argmin(abs(t-tm)))
        assert abs(t[i]-tm) < 1e-12
        N = factor * np.sum(f, dtype=np.longdouble)
        np.testing.assert_allclose(N, series['integrated_moms'][1][i], rtol=2e-13)
        S = source[i]

        def integral(x):
            return float(factor * np.sum(x, dtype=np.longdouble) / S)

        a = read(path, f'fdot_multiplier_{frame}').get_values()[..., 0].copy()
        rate = read(path, f'damping_rate_{frame}').get_values()[..., 0].copy()
        bar = read(path, f'damping_fbar_{frame}').get_values()[..., 0].copy()
        assert a.shape == rate.shape == bar.shape == f.shape
        assert np.all(a > 0) and np.all(a <= 1)
        g = rate * (f - bar)
        bar_rhs = a*g if bar_multiplied else g
        row = dict(run=label, frame=frame, time=tm, changeset=version,
                   source=S, N_over_S=float(N/S), flux_over_S=flux[i]/S,
                   fdot_diagnostic_over_S=ndot[i]/S,
                   fbardot_rhs_over_S=integral(bar_rhs),
                   unweighted_filter_exchange_over_S=integral(g),
                   multiplier_min=float(a.min()), multiplier_max=float(a.max()),
                   residual_with_bar_estimate=1-(flux[i]+2*ndot[i])/S,
                   residual_with_bar_rhs=1-(flux[i]+ndot[i])/S-integral(bar_rhs),
                   fbardot_diagnostic_over_S=None,
                   fdot_rhs_over_S=None, weighted_fdot_over_S=None,
                   multiplier_correction_over_S=None, corrected_residual=None,
                   fdot_rhs_minus_diagnostic_over_S=None)
        if 'fbardot_integrated_moms' in series:
            row['fbardot_diagnostic_over_S'] = series['fbardot_integrated_moms'][1][i]/S
        if (path / f'zzim-ion_fdot_{frame}.gkyl').exists():
            d = read(path, f'fdot_{frame}').get_values()[..., 0].copy()
            row.update(fdot_rhs_over_S=integral(d), weighted_fdot_over_S=integral(d/a),
                       multiplier_correction_over_S=integral(d/a-d+g-bar_rhs),
                       corrected_residual=1-flux[i]/S-integral(d/a+g),
                       fdot_rhs_minus_diagnostic_over_S=integral(d)-ndot[i]/S)
        rows.append(row)
    bfinal = (series['fbardot_integrated_moms'][1][-1]
              if 'fbardot_integrated_moms' in series else ndot[-1])
    dtdata = pg.GData(str(path / 'zzim-dt.gkyl'))
    dts = dtdata.get_values()[:, 0]
    # Use ordinary late steps, excluding the final truncated step.
    dtlate = float(np.median(dts[-100:-1]))
    summary = dict(path=str(path.resolve()), final_time=float(t[-1]),
                   final_flux_over_S=float(flux[-1]/source[-1]),
                   final_fdot_over_S=float(ndot[-1]/source[-1]),
                   final_fbardot_over_S=float(bfinal/source[-1]),
                   final_fbardot_is_estimate='fbardot_integrated_moms' not in series,
                   final_requested_residual=float(1-(flux[-1]+ndot[-1]+bfinal)/source[-1]),
                   typical_late_dt=dtlate,
                   epsilon_N_over_dtS=float(np.finfo(float).eps*rows[-1]['N_over_S']/dtlate),
                   max_abs_corrected_residual=max(abs(r['corrected_residual']) for r in rows
                                                  if r['corrected_residual'] is not None))
    return rows, summary


def main():
    runs = [(RUN, list(range(100, 115, 2))+[115], 'extra cycles'),
            (RUN.parent.parent / 'a-thesis-simulations/1x-beams',
             list(range(50, 75, 4))+[75], 'comparison')]
    all_rows, summaries = [], {}
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex='col', layout='constrained')
    for col, (path, frames, label) in enumerate(runs):
        rows, summary = audit(path, frames, label)
        all_rows.extend(rows)
        summaries[label] = summary
        start = rows[0]['time']
        x = np.array([r['time']-start for r in rows])*1e6
        ax = axes[0, col]
        ax.plot(x, [r['residual_with_bar_estimate']*100 for r in rows], 'o-',
                label=r'$\dot{\bar N}\approx\dot N$')
        ax.plot(x, [r['residual_with_bar_rhs']*100 for r in rows], 's-',
                label='Actual filter RHS')
        ax.axhline(0, color='black', lw=.7)
        ax.set(title=label, ylabel='Unweighted residual (% of source)')
        ax.legend(fontsize=9)
        good = [r for r in rows if r['corrected_residual'] is not None]
        axes[1, col].semilogy([(r['time']-start)*1e6 for r in good],
                              [abs(r['corrected_residual']) for r in good], 'o-')
        axes[1, col].set(xlabel='Time since first sampled late FDP frame (µs)',
                         ylabel='Absolute corrected residual / source')
        for ax in axes[:, col]:
            ax.grid(alpha=.25)
    fig.suptitle('Particle balance: account for the phase-space time-rate multiplier')
    fig.savefig(HERE / 'particle_balance.png', dpi=180)
    fig.savefig(HERE / 'particle_balance.pdf')
    with (HERE / 'frame_balance.csv').open('w') as out:
        writer = csv.DictWriter(out, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    (HERE / 'summary.json').write_text(json.dumps(summaries, indent=2)+'\n')
    print(json.dumps(summaries, indent=2))
    for row in all_rows:
        print(row['run'], row['frame'], 'bar RHS/S', row['fbardot_rhs_over_S'],
              'unweighted', row['residual_with_bar_rhs'], 'corrected', row['corrected_residual'])


if __name__ == '__main__':
    main()
