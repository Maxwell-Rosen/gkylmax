# Resolution profile comparisons

All plotting scripts use `plotting.mplstyle`: LaTeX text, a serif font, 12 pt
body and legend text, 18 pt axes titles, 22 pt axes labels, and 14 pt tick
labels. Rendering requires a working LaTeX installation and `dvipng` for PNG
output. Figure sizes allow room for the larger labels.

Run in the postgkyl environment, from this directory:

```bash
python z-scan/compare.py
python vpar-scan/compare.py
python mu-scan/compare.py
```

Each script discovers all numeric resolution folders beside itself, sorts them
numerically, and includes `../1x-beams` relative to this directory. Paths are
resolved from the scripts, so they also work from another working directory.
Shared plotting code lives in `compare_profiles.py`.

The default frame is 75. The four mapped profiles (density, parallel flow,
parallel temperature, perpendicular temperature) use the original split
linear/log layout. Output is `res-scan-z.pdf`, `res-scan-vpar.pdf`, or
`res-scan-mu.pdf` in the corresponding scan directory.

Missing moment or coordinate-map files are reported and skipped. At least two
available runs are required, and their frame times must match. Rerun after
unfinished simulations produce output to include them automatically. Legend
resolutions come from directory names; the separate baseline is labeled
`Original (1x-beams)`.

Options include `--frame 75`, `--output /path/to/figure.pdf`,
`--baseline /path/to/run`, and `--no-baseline`.

## Convergence diagnostics

Use the shared convergence script in the postgkyl environment:

```bash
python convergence.py --axis z
python convergence.py --axis vpar
python convergence.py --axis mu
```

The default axis is `z`. Each command writes seven figures as PDF and PNG,
CSV diagnostics, and `run_report.txt` to `<axis>-scan/convergence-plots/`.
The figures cover inventory and wall-flux histories, final-relaxation zooms,
convergence versus resolution, relative differences, and final ion profiles.
All axes use the same diagnostics and normalization.

Particle boundary-flux histories show `2 * bflux / source` for each wall,
using integrated ion particle rates. The both-wall convergence panel shows
`(lower bflux + upper bflux) / source`, equivalent to twice either wall flux
for symmetric boundaries. A dotted line at 1 marks source-loss balance.
The source and normalized ratios are exported alongside the raw fluxes in CSV.

By default, comparisons use frame 75, a trailing 15 µs averaging window,
and a 60 µs final-relaxation zoom. The highest included resolution is the
finite-resolution reference; bars show temporal ranges, not uncertainty.
The `../1x-beams` baseline is included, giving velocity scans their
`Nvpar=64` or `Nmu=32` reference with the current data. Input files must match
apart from the scanned resolution. Actual velocity cell counts come from
saved run statistics when available.

Options include `--no-baseline`, `--exclude 32`, `--frame 65`,
`--window-us 15`, `--zoom-us 60`, `--root /path/to/scan`, and
`--output /path/to/plots`. Paths default relative to the script, so it can
also be invoked from another working directory.

## Combined publication figures

Run `python publication_convergence.py` in the postgkyl environment to build
one overview and four supplementary figures for all three scans together.
The script reads the raw diagnostics using the definitions in `convergence.py`,
checks matching physics and fixed unscanned cell counts, and loads the shared
baseline once. Existing per-scan CSV files are not used as inputs.

Outputs are in `publication-figures/`, each as a vector PDF, editable SVG, and
300 dpi PNG:

1. `01_convergence_overview`: the same six ion quantities as the absolute
   diagnostics, organized by scan and by inventory/wall loss. Points show
   `100 * (Q_final - Q_reference_final) / abs(Q_reference_final)`, using the
   finest available grids for the z and vpar scans and `Nmu=32` for the mu scan
   as references. Finer mu runs remain plotted relative to `Nmu=32`.
   Only final values are shown, with no error bars. Both rows use independent
   linear y-axis scales with one centered shared y-label. The figure has no
   title and a compact one-row, six-column legend. Each quantity uses a distinct
   color, line style and marker, consistently across the three scan columns.
   `01_final_state_convergence.csv` records the values and relative differences.
2. `02_absolute_diagnostics`: the paper figure, showing final values only:
   ion number `Ni`, ion kinetic energy `Ei`, central ion density,
   source-normalized outward ion loss `Gamma_i / integral(S_i d^3x d^3v)`, mean ion
   wall heat flux `qi = (qi_lower + qi_upper)/2`, and summed ion wall power
   `Pi`. Units are SI as specified by the simulation author; ion number is
   a particle count, energy is in J, density in m^-3, boundary flux/source is dimensionless,
   heat flux in W m^-2, and power in W. The labels show display multipliers.
   There are no error bars, temporal means, electron contributions, daggers,
   panel text, titles, column subtitles or footnotes.
   `02_final_state_data.csv` contains the plotted values without display scaling,
   plus the separate summed particle loss and integrated source, both in s^-1.
3. `03_final_relaxation`: the last 60 µs of the same six absolute ion
   quantities as figure 2, in the same order and display units: `Ni` in
   10^19, `Ei` in kJ, central ion density in 10^19 m^-3, summed outward ion
   loss divided by the integrated source, mean ion wall heat flux in MW m^-2,
   and summed ion wall power in MW. No quantity is normalized by another
   simulation or temporal mean. A dashed line marks unity in the source-loss
   ratio row. Time is `t - t_end` in µs, with zero at the simulation end.
   Every completed simulation is included, with shared vertical scales per
   row, simple mathematical labels, and compact resolution legends matching
   figure 4's colors and line styles. There are no titles, subtitles, panel
   annotations, daggers, bottom explanatory text or averaging-window shading.
4. `04_final_profiles`: density in m^-3, flow in 10^6 m s^-1, and both
   temperatures in keV at the final comparison frame. Each profile is split
   at `z=0`, with a linear left half and a logarithmic right half. Nonpositive
   samples are masked on the log half, preserving gaps between valid samples
   and between DG cells. Signed values and red negative-temperature shading
   remain on the linear half. Grey dashed lines mark both mirror throats at
   `z=+/-0.98 m`, drawn above the simulation curves, with one shared legend
   entry `Z=+/-Z_m` in the `Nz` legend and no inline labels. Each column starts
   the same palette in increasing resolution order: blue, orange, green,
   purple, then brown. The production baseline `(400, 64, 32)` is purple in
   all three columns. Scan legends show completed resolutions. This figure
   has no title.
   The perpendicular-temperature log halves show the full positive data range,
   with 8% padding in log space, so the exhaust tails remain visible. Other log
   minima use the fifth percentile of positive right-half samples, rounded
   down to a decade: currently `1e-4` for flow and `1e-1` keV for parallel
   temperature. Linear temperature minima use the first percentile (or zero,
   if smaller), with 8% range padding, to keep isolated negative spikes from
   stretching the axes. Bounds pool each distinct run once and are shared
   across scans. Outliers remain in the data and are clipped by the viewport.
5. `05_particle_balance`: the `Gamma_i / integral(S_i d^3x d^3v)` row of
   figure 2 isolated into one row of three columns (`Nz`, `Nvpar`, `Nmu`).
   Every completed resolution is included on a shared ratio scale. A dashed
   horizontal line marks the expected value of 1, and open diamonds identify
   the shared baseline `(400, 64, 32)`. The right-hand axis shows the signed
   percent deviation `100 * (ratio - 1)`, so agreement with unity can be read
   directly. These are the same final values as figure 2, without temporal
   averaging or normalization by another simulation. The figure has no title.
   `05_particle_balance_data.csv` records every ratio, signed deviation,
   summed outward ion loss and integrated source.

All five figures use matching typography, mathematical labels, and clean
axes; figure 2 retains its master title. They have no error bars, daggers,
panel annotations, subtitles, or bottom explanatory text. Figure 1 uses the
same ion-only final quantities, summed particle loss, and mean ion wall heat
flux as figure 2, as do the time histories in figure 3; figure 4 contains
local ion profiles only.

The default data have 11 distinct completed simulations and 13 scan
memberships because the baseline `(Nz, Nvpar, Nmu) = (400, 64, 32)` appears
in all three scans. The `Nz=128` directory was deleted and is not included
in the final-state figures. Old CSVs for this resolution are not treated as current
results. `quality_checks.csv` reports negative samples in both temperature
components and reconstruction checks. A finite-grid reference is not a
continuum solution, and temporal
range bars are not uncertainty estimates.

`captions_and_methods.txt` provides figure captions, normalization details,
scan membership, and numerical differences for the finest pairs.
`all_diagnostics.csv` exports all 29 diagnostic quantities for every scan
member. `manifest.json` and `diagnostic_snapshot.npz` preserve the plotted
data, full time histories, profiles, paths, cell counts, and code changesets.

Use `python publication_convergence.py --reuse-data` to replot that snapshot
without rereading the simulations. Omit `--reuse-data` to incorporate new
output. Options include `--frame`, `--window-us`, `--zoom-us`, `--dpi`,
`--root`, `--baseline`, and `--output`. The data collection step checks the
same equilibrium filename while allowing its directory to differ between
machines, matching the scan-local convergence scripts.

Regenerate just the paper figure with `python publication_convergence.py --figures 2`.
Regenerate the overview and final profiles with
`python publication_convergence.py --figures 1 4`.
Regenerate only the final profiles with `python publication_convergence.py --figures 4`.
Regenerate only the final-relaxation histories with
`python publication_convergence.py --figures 3`.
Regenerate only the source-loss balance figure with
`python publication_convergence.py --figures 5`.
Add `--reuse-data` to render from the saved snapshot.

### Definitions checked against Gkeyll

The local Gkeyll source at `/home/mr1884/gkeyll` has changeset
`1dd609d0006c5bc2a0bc1bddb888d446ef2ab98c`, matching every simulation.
The integrated moment kernels, diagnostic reduction code, and saved-moment
Jacobian handling were inspected at that revision.

The ion inventory is

\[
N_i=\int n_i\,d^3x=\int f_i\,d^3x\,d^3v.
\]

It is component 0 of `ion_integrated_moms`. Gkeyll calculates cell integrals
and sums them over the domain and MPI ranks. `convergence.py` independently
checks this value by integrating the saved `ion_M0` with the geometric
Jacobian.

The plotted ion energy is kinetic energy:

\[
E_i=\frac{m_i}{2}\int M_{2i}\,d^3x
=\int\left(\frac{m_i v_\parallel^2}{2}+\mu B\right)f_i\,d^3x\,d^3v,
\qquad
M_{2i}=\int\left(v_\parallel^2+\frac{2\mu B}{m_i}\right)f_i\,d^3v.
\]

This includes bulk parallel-flow energy and both thermal components. In
continuum moment notation it is

\[
E_i=\int n_i\left(\frac{m_i u_\parallel^2}{2}
+\frac{T_\parallel}{2}+T_\perp\right)d^3x,
\]

where temperatures are in energy units (J). This formula does not assume a
Maxwellian distribution. The plotted values use the saved second moment
directly, avoiding extra numerical products of primitive variables.
`convergence.py` multiplies `ion_M2` by the configuration-space geometric
Jacobian, integrates the product of the two p=1 modal expansions exactly,
and multiplies by `m_i/2`.

The input selects `GKYL_F_MOMENT_HAMILTONIAN`. Its three integrated
components are particle number, parallel momentum, and the integral of
`m_i*v_parallel^2/2 + mu*B + e_i*phi` against the distribution. Therefore
the third component of `ion_integrated_moms` is
`E_i + integral(e_i*n_i*phi d^3x)`, with `e_i` denoting ion charge; it is
not kinetic energy alone. The script independently checks that identity.

The boundary ratio is

\[
\frac{\Gamma_i}{\int S_i\,d^3x\,d^3v},
\qquad \Gamma_i=\dot N_{i,\mathrm{lower}}+\dot N_{i,\mathrm{upper}}.
\]

Both wall losses are outward-positive. The denominator is component 0 of
`ion_source_integrated_moms`, calculated by applying the same integrated
moment operator to the phase-space source. The ratio is evaluated at the
final time, and unity means that boundary losses equal particle injection.
All spatial and velocity measures include the geometry and velocity-map
Jacobians used by the gyrokinetic discretization.

The exact source functions and line locations are recorded in
`publication-figures/captions_and_methods.txt` under GKEYLL SOURCE VERIFICATION.
