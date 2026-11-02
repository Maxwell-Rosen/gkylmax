"""Common visual style and ion diagnostic definitions for publication figures."""
from __future__ import annotations
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/thesis-publication-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import diagnostics as cv
ROOT = Path(__file__).resolve().parent

def style():
    plt.style.use(ROOT / 'plotting.mplstyle')
    plt.rcParams.update({
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': .7, 'axes.edgecolor': '#687078',
        'legend.frameon': False,
        'lines.linewidth': 1.5, 'savefig.facecolor': 'white',
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'axes.formatter.use_mathtext': True, 'axes.formatter.useoffset': False,
    })


def finish(fig, args, name):
    for extension in ('pdf', 'svg', 'png'):
        fig.savefig(args.output / f'{name}.{extension}', dpi=args.dpi,
                    bbox_inches='tight', pad_inches=.15)
    plt.close(fig)
    print(f'Wrote {name} (PDF, SVG, PNG)', flush=True)


def paper_axis(ax):
    """Shared uncluttered axes for the publication figures."""
    ax.set_axisbelow(True)
    ax.grid(axis='y', color='#dce1e5', lw=.6)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))


def diagnostic_panels(*, sharex=False):
    """Lay out diagnostic pairs, spanning the final odd panel across both columns."""
    fig = plt.figure(figsize=(12, 10))
    grid = fig.add_gridspec((len(ABS_ROWS) + 1) // 2, 2)
    axes = []
    for i in range(len(ABS_ROWS)):
        cell = grid[i // 2, :] if i == len(ABS_ROWS) - 1 and i % 2 == 0 else grid[i // 2, i % 2]
        axes.append(fig.add_subplot(cell, sharex=axes[0] if sharex and axes else None))
    if sharex:
        for ax in axes[:-1]:
            ax.tick_params(axis='x', labelbottom=False)
    return fig, np.asarray(axes)


# Paper figure: output key, SI axis label, display scale, and diagnostic weights.
# Number is a particle count; only the local density has units of m^-3.
# Gamma_i is the summed outward loss rate; the displayed ratio divides by
# the zeroth integrated moment of the phase-space source.
# q_i is the total wall-integrated ion heat flux, in W (displayed in kW).
# Pi_lo/Pi_hi are the saved reader keys for the two integrated wall fluxes;
# qi_lo/qi_hi are area-normalized quantities and are not plotted here.
ABS_ROWS = [
    ('Ni', r'$N_i=\int n_i\,d^3x$ [$10^{19}$]',
     1e-19, [('Ni', 1)]),
    ('Ei_J', r'$E_i=\frac{m_i}{2}\int M_{2i}\,d^3x$ [kJ]',
     1e-3, [('Ei', 1)]),
    ('n0_m-3', r'$n_i(z=0)$ [$10^{19}$ m$^{-3}$]',
     1e-19, [('n0', 1)]),
    ('Gamma_i_over_source', r'$\Gamma_i\,/\!\int S_i\,d^3x\,d^3v$',
     1, [('particle_balance', 1)]),
    ('qi_W', r'$q_i$ [kW]',
     1e-3, [('Pi_lo', 1), ('Pi_hi', 1)]),
]


def history_values(data, metrics, start, end):
    """Combine the same diagnostics as figure 2, preserving native sample times."""
    times = np.unique(np.concatenate([
        [start, end],
        *(data[key][0][(data[key][0] > start) & (data[key][0] < end)]
          for key, _ in metrics),
    ]))
    values = sum(weight * cv.sample(data[key], times) for key, weight in metrics)
    return times, values


