# Beam alpha scan

This scan varies `alpha_oap`, the factor multiplying collisionless terms in the
OAP phases. `core/sim.c` comes from the `1x-beams` base case. All other physics,
the 400 × 64 × 32 ion grid, phase durations, and 75-frame schedule are unchanged;
`alpha_fdp` remains 1.

| Case under `alpha-scan/` | Alpha | Existing run and data |
| --- | --- | --- |
| `alpha-0.3e-3` | 0.0003 | Moved from `1x-beams-higher-alpha`; frames 0–75 |
| `alpha-1e-3` | 0.001 | Moved from `1x-beams`; base case, frames 0–75 |
| `alpha-2.5e-3` | 0.0025 | New input; ready to run |
| `alpha-5e-3` | 0.005 | New input; ready to run |
| `alpha-10e-3` | 0.01 | Moved from `1x-beams-lower-alpha`; geometry only |

Entire original directories were moved, including diagnostics, checkpoints,
logs, binaries, and the older `256/` data where present. Compatibility symlinks
at all three original paths keep existing analysis and continuation references
working without duplicating the data. The historical higher/lower names do not
describe the numerical ordering of `alpha_oap`.

The earlier `10e-3` job (`slurm-6671.out`) segfaulted during initialization and
wrote no ion or field checkpoints. Its existing files are preserved; it still
needs a successful run. Building the input does not establish that the earlier
runtime failure has been resolved.

## Prepare, build, and submit

Run from this directory (or invoke the script by its full path):

```bash
# Show the submissions without changing files or submitting jobs.
bash copy-run.sh --submit --dry-run

# Build and submit the unfinished 2.5e-3, 5e-3, and 10e-3 cases on Stellar AI.
bash copy-run.sh --submit

# Build only, or submit a selected case.
bash copy-run.sh --build
bash copy-run.sh --submit 5e-3

# Prepare another alpha using the same base input.
bash copy-run.sh --prepare 7.5e-3

# Use the Perlmutter modules and job script when running on Perlmutter.
bash copy-run.sh --machine perlmutter --submit
```

With no options, `copy-run.sh` only prepares missing folders. It never refreshes
an existing case from `core/`, overwrites its source, or copies output into a
different alpha case. An existing folder must contain the matching `alpha_oap`.
Use the listed spelling for existing alpha values.

Build and submit modes skip cases with nonempty final ion and field checkpoints
at frame 75. This completion check assumes the unchanged 75-frame schedule.
Unfinished cases are rebuilt against the current installed Gkeyll library and
shared headers. The job scripts restart from the latest saved ion frame, or start
from the beginning when none exists. The helper does not detect already queued
jobs; check the queue before submitting the same unfinished case again.

Each case has its own `sim.c`, `Makefile`, and Stellar AI and Perlmutter job
scripts. Run a job script from its case directory. The existing Makefile's
installation/header paths and the geometry path in `sim.c` are retained.

## Compare saved profiles

```bash
python3 compare_profiles.py
python3 compare_profiles.py --frame 65 --output alpha-profiles-65.pdf
```

The comparison discovers numeric `alpha-*` folders, plots ion density, parallel
flow, and parallel/perpendicular temperatures, and labels each curve with its
actual alpha. It uses each run's own geometry mapping, converts temperatures
using the ion mass in the input, and checks that saved times agree. Missing
frames are reported and skipped, so the existing `0.3e-3` and `1e-3` outputs can
be compared immediately; new cases join automatically once their output exists.
The default output is `alpha-profiles.pdf`. Requires NumPy, Matplotlib, and
Postgkyl in the active Python environment.
