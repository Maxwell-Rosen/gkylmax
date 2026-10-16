# Resolution profile comparisons

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
