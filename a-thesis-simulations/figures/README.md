# Publication figures

Each study owns its entry script and `publication-figures/` outputs. Shared
`diagnostics.py`, `publication_style.py`, and `plotting.mplstyle` preserve the
diagnostic definitions, LaTeX typography, colors, and vector export settings of
the original resolution analysis. `publication_scans.py` handles the three new
studies. Simulation data stay in their original directories.

Run from this directory with the existing `pgkyl` Python environment, which
provides NumPy, Matplotlib, Postgkyl, and the required LaTeX tools:

```bash
python generate_all.py
python generate_all.py --reuse-data
python generate_all.py --only alpha-scan extra-time-dilation
```

The first command reads raw diagnostics. `--reuse-data` explicitly replots the
recorded NPZ snapshot, preserving the selected endpoint of an ongoing run.
Omit it to incorporate new output. `--root` on `generate_all.py` selects the
parent simulation directory. Scripts work independently of the current working
directory.

| Study | Script | Data |
| --- | --- | --- |
| Resolution | [publication_convergence.py](resolution-scan/publication_convergence.py) | `1x-beams-resolution-scan/{z,vpar,mu}-scan` and `1x-beams` |
| OAP alpha | [publication_figures.py](alpha-scan/publication_figures.py) | Completed `1x-beams-alpha-scan/alpha-scan/alpha-*` cases and `1x-beams` |
| Extra cycles / FDP duration | [publication_figures.py](extra-cycles-fdp-length/publication_figures.py) | `1x-beams`, `1x-beams-ecycles`, `1x-beams-anneal`, `1x-beams-fdp-continuation`, and the selected complete moment frame in `1x-beams-fdp-continuation-extra` |
| Extra time dilation | [publication_figures.py](extra-time-dilation/publication_figures.py) | Completed five-cycle runs in `1x-beams` and `1x-beams-extra-dilated` |

The original `1x-beams-resolution-scan/publication_convergence.py` remains a
compatibility entry point. Its default output now goes under `figures/`.
Previous output files in the simulation directory are retained as an archive.

The endpoint and history comparisons show five diagnostics. The ion heat flux `q_i` is the sum
of the energy fluxes integrated over both walls, displayed in **kW** and
exported as `qi_W` in plotted CSVs. It is not divided by wall area. The separate
power diagnostic has been removed from the figures. Raw snapshots retain the
reader's original keys (`Pi_lo` and `Pi_hi` for the integrated wall fluxes) so
existing snapshots remain usable with `--reuse-data`.

## Figures and data

The resolution script preserves the existing five figures and CLI options,
including the 60 μs relaxation window. Each new study writes the following six
figures. The alpha and cycle/FDP studies also each write a study-specific figure 7:

1. Endpoint ion number, kinetic energy, central density, source-loss ratio,
   and total wall heat flux.
2. The last 15 μs of these quantities.
3. The last 60 μs.
4. The full last continuous FDP interval, using each run's true time extent.
   For the alpha scan and extra cycles / FDP duration studies, this is a single
   column of five panels overlaying all runs over the last 90 μs on a shared
   time axis. This covers the full final continuous FDP interval for every
   alpha case.
5. Endpoint density, parallel flow, parallel temperature, and perpendicular
   temperature, with linear negative-z and logarithmic positive-z axes.
6. A separate source-loss ratio plot with deviation from unity.
7. **Alpha scan only:** a single row with the density difference
   `n_OAP - n_FDP` versus z in units of 10¹⁷ m⁻³ on the left and
   `integral |n_OAP - n_FDP| dz` versus alpha on the right.
   The figure has no title or panel subtitles; all five alpha values share
   one legend row. The alpha figures include the completed `2.5e-3` case
   alongside `0.1e-3`, `0.3e-3`, `0.5e-3`, and the `1e-3` baseline.
   The standard schedule selects frames **45 and 75**, separated
   by **90 μs**. Frames come from each run's phase schedule.

The cycle/FDP study's **figure 7**, `07_annealing_stages`, shows ion number,
kinetic energy, and source-loss ratio through the three annealing OAP stages
and final FDP. Each stage has its own time axis and the rows share vertical
scales. Dashed reference lines show the baseline endpoint. The accompanying
schedule CSV records nominal and saved stage boundaries.

The alpha study's figure 7 uses `ion_M0` and Postgkyl native modal subtraction
for the plotted density difference. The additional dimensionless relative
profile in the CSV uses weak division by `n_FDP`. FDP density must be positive
throughout every cell; no density floor is imposed.

For the absolute integral, the exact mean of `|delta n|` is computed in each
p=1 cell, splitting at any zero crossing inside it. These modal means are
weakly multiplied by `dz/dz_comp`, obtained with Postgkyl's local DG derivative
of `geo_corn_mc2nu_pos_deflated`, then integrated with Postgkyl's native DG
kernel. The installed `integrate(op="abs")` shortcut takes the absolute cell
mean and would miss cancellation within a cell, so it is not used here.
The measure is physical `dz` over the full domain, with units of m⁻².
The absolute integral accumulates differences without sign cancellation.
The summary CSV retains signed integrals and adds absolute integrals and their
percentages of the FDP line density. Profile CSV rows preserve both one-sided
DG traces and include the weak relative difference.

Every figure is exported as vector PDF, vector SVG, and 300 dpi PNG. Plotted
values are also exported as CSV. Each study includes `manifest.json`,
`diagnostic_snapshot.npz`, and `captions_and_methods.txt`; the new studies also
write `sampling.csv` and `window_statistics.csv`. Read the captions alongside
the figures when preparing a manuscript. The manifest records paths, phase
parameters, solver revisions, reconstruction checks, source hashes, and the
size and modification time of every consumed diagnostic.

For example:

```bash
python alpha-scan/publication_figures.py --figures 2 4 5
python alpha-scan/publication_figures.py --figures 7
python extra-cycles-fdp-length/publication_figures.py --reuse-data
python extra-cycles-fdp-length/publication_figures.py --refresh-runs annealing
python resolution-scan/publication_convergence.py --figures 3 --zoom-us 15 --output /tmp/resolution-last15
```

For the resolution script, `--root` means the resolution-scan data directory;
for the three new scripts, it means the parent simulation directory. All four
scripts support `--output`, `--dpi`, `--figures`, and `--reuse-data`.
Older alpha snapshots without the relative/absolute OAP comparison must be regenerated once
without `--reuse-data` before plotting figure 7 from the snapshot.

`--refresh-runs UID ...` reads only the named runs and retains the others from
the existing output snapshot. It is mutually exclusive with `--reuse-data`.
The annealing addition preserves the four previously recorded comparison runs,
including continuation frame 13: newer raw continuation diagnostics contain
backward time jumps from restart output and fail the reader's validation.
Per-run collection times and source inventories are retained in the manifest.
Use `--reuse-data` to reproduce this complete five-run set. Older cycle/FDP
snapshots need `--refresh-runs annealing` once to add the annealing data.

## Time windows and interpretation

The standard schedule has a **75 μs final extra FDP**, immediately preceded by
a **15 μs FDP with identical settings**. Figure 4 includes the combined **90 μs**
continuous interval. Its dotted line marks the configured extra-phase boundary.
The extra cycles / FDP duration study restricts every run in figure 4 to the
last up to 90 μs, retaining the filename `04_last_full_dynamics`; the other studies
show the entire final continuous interval.
The recorded annealing run has four regular cycles followed by OAP stages with alpha
0.001, 0.01, and 0.1 lasting 100, 10, and 1 ms. Its final FDP is **75 μs**
(frames 55–80), ending at **0.511135 s**. Its figure-4 curve begins at -75 μs;
the preceding annealing stage is not counted as FDP. The `256/` archive has a
different spatial resolution and is excluded from this comparison. The selected
annealing run uses the same 400 × 64 × 32 grid as the baseline.
The revised `1x-beams-anneal/sim.c` keeps four regular OAP/FDP pairs. The
first ramp stage supplies the fifth 100 ms OAP at alpha=0.001, matching the
baseline's five main OAPs; alpha=0.01 and 0.1 add 10 ms and 1 ms of annealing.
The **90 μs** final FDP (frames 55–85) now matches the baseline's continuous
final FDP. The endpoint is **0.511150 s**, exactly 11 ms after the baseline's
0.500150 s. The final frame spacing remains 3 μs, and all 12 allocated phases
are initialized. Rebuild the executable before running or resuming this input.
The saved figures and `--reuse-data` snapshot still describe the recorded run;
refresh the annealing data only after outputs for the revised schedule are available.
The last 60 μs is a separate view, not a claim about the phase duration. For the
extra-dilation run these durations are 25 ms, 5 ms, and 30 ms, respectively.
Continuations reset their clock; their full-interval panels show the available
continuation segment without stitching it to the parent run.

All times are **saved simulation-clock times**. The CFL/time-dilation cap is
not treated as a global conversion to physical time. Histories are aligned to
each run's own selected endpoint. Endpoint comparisons therefore compare the
outcomes of the stated schedules, not a common elapsed physical time.

The extra-dilation comparison changes both `cfl_factor_times_omega_max`
(0.1 to 0.001) and the FDP durations. It does not isolate the effect of the cap.
The separate `1x-beams-extra-dilated-fdp-continuation` is recorded as excluded;
at collection it contains only frame 0.

Energy and central density use saved spatial moments: approximately every
3 μs for the standard runs, 10 μs for the 300 μs continuation, 100 μs for the
ongoing continuation, and 1 ms for the extra-dilation run. Dots mark actual
samples. If a bracketing gap exceeds half the displayed window, only native
samples appear; no short transient is invented by interpolation. Statistics
still use piecewise-linear interpolation, so consult `sampling.csv` before
interpreting a short-window mean or range as evidence of relaxation. These
ranges describe temporal variation, not statistical uncertainty.

Incomplete alpha cases are excluded from endpoint comparisons. The ongoing
long continuation is explicitly labeled and uses an open diamond in the
endpoint figure. Its current saved state has negative parallel-temperature
samples, which remain visible on the linear profile panels. Nonpositive values
are masked on logarithmic panels. The plots do not certify a physically valid
steady state. Solver revisions differ between some runs and are listed in the
manifest and methods.

## Checks

The shared reader rejects nonfinite data, backward time, unsupported bases,
missing frames, and temporal extrapolation. It reconstructs particle number and
Hamiltonian energy independently. New studies also check the reviewed input
differences, saved endpoint times, phase-boundary times, fixed resolution,
geometry, continuation parent paths, and whether inputs changed while reading.

```bash
python -m unittest discover -s . -p 'test_*.py'
```

`validation.txt` records the validation of the supplied figures, including
comparison of the refactored resolution diagnostics with the original reader.
