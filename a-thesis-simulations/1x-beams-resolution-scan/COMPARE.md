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

The default frame is 65. The four mapped profiles (density, parallel flow,
parallel temperature, perpendicular temperature) use the original split
linear/log layout. Output is `res-scan-z.pdf`, `res-scan-vpar.pdf`, or
`res-scan-mu.pdf` in the corresponding scan directory.

Missing moment or coordinate-map files are reported and skipped. At least two
available runs are required, and their frame times must match. Rerun after
unfinished simulations produce output to include them automatically. Legend
resolutions come from directory names; the separate baseline is labeled
`Original (1x-beams)`.

Options include `--frame 50`, `--output /path/to/figure.pdf`,
`--baseline /path/to/run`, and `--no-baseline`.
