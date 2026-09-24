# New mapping comparison, 2026-09-22

The new 256-cell run does **not establish improved convergence overall**.
It agrees better with the original 400-cell reference for particle loss and
mean wall powers, but substantially changes stored energy and density.
Agreement with this finite-resolution reference is not a measurement of
continuum error; the new result could expose a mapping-dependent error shared
by the original scan. A refinement scan using the new mapping is needed to
resolve that question.

All six runs are compared at frame 65, t = 0.50012 s. The table uses means
over the same final 15 microseconds, [0.500105, 0.50012] s. Entries are signed
percentage differences from the **original-mapping 400-cell** run.

| Quantity | Original 256 | Original 320 | New mapping 256 |
| --- | ---: | ---: | ---: |
| Integrated ion number | -0.00652% | -0.00125% | +1.68136% |
| Total particle kinetic energy | +0.03309% | +0.02351% | +12.25248% |
| Central ion density | +0.00048% | +0.00078% | +0.75734% |
| Particle loss, both walls | -1.72718% | -0.74695% | -0.19275% |
| Ion wall power, both walls | -8.64887% | -3.22954% | +2.93637% |
| Total wall power, both walls (model) | -8.03881% | -3.01074% | +2.66058% |

The new mean particle loss is closer to the reference than either original
256 or 320. Its mean wall powers are also closer, but their instantaneous
final values overshoot the reference by +6.14% (ion) and +5.62% (total),
larger absolute differences than original 320 (-3.16% and -2.94%, respectively).
In the new run the total wall-power
max-minus-min range is 5.47% of its mean over the averaging window, versus
0.84% for original 256, 0.12% for original 320, and 0.24% for the reference.
Thus the apparent wall-power improvement depends appreciably on averaging
time. Range bars are temporal variation, not statistical error bars.

The new stored-energy difference is persistent: +12.2525% in the mean and
+12.2532% at the final time, while its temporal range is only 0.00286% of
the mean. The final profiles show a higher central parallel temperature,
consistent with this energy shift. Neither original 256, original 320,
original 400, nor new 256 has negative sampled final parallel temperatures;
the check uses two interior samples per cell, not a positivity proof.

The grid is redistributed, not uniformly refined. Both runs use
`GKYL_PMAP_CONSTANT_DB_NUMERIC`, with `map_strength = 1.0`. The new input
changes both `maximum_slope_at_min_B` and `maximum_slope_at_max_B` from 2 to
4, and `gaussian_std` from 0.5 to 0.25. Other `sim.c` content agrees after
normalizing Nz and removing comments/whitespace.

| Physical cell width | Original 256 | Original 400 | New mapping 256 |
| --- | ---: | ---: | ---: |
| At center | 0.033994 m | 0.021751 m | 0.046434 m |
| At outer wall | 0.039063 m | 0.025000 m | 0.078125 m |
| Smallest cell | 0.004073 m | 0.002619 m | 0.003822 m |

The new grid has slightly smaller minimum cells than original 256 but coarser
central and outer cells. No single effective Nz describes it. This also means
that ranking the grids by one global scalar can hide losses in local accuracy.

Start with [relative differences](relative_differences.pdf),
[resolution scan](convergence_vs_resolution.pdf),
[final profiles](final_profiles.pdf), and [grid spacing](grid_spacing.pdf).
Every plot is also available as PNG. The original seven figure types include
all six runs; the grid-spacing comparison is an additional eighth figure.
Scalar data and temporal ranges are in `convergence_summary.csv`, with the
new run identified by `mapping=nunif`.

Validation: all diagnostic code changesets are `02623e26d7a5`, comparison
times agree, and wall area factors agree. Independent ion-number and
Hamiltonian reconstructions agree at approximately 1e-15 relative precision,
including the new mapping. All original scalar results exactly match the
original folder's CSV, and the five original frame CSVs are byte-identical.
See `run_report.txt` and `input_checksums.txt` for provenance. This checks the
diagnostic reconstruction, not the simulation's discretization error.

Reproduce from the parent `z-scan` directory with the postgkyl environment:

```bash
python plot_convergence.py --nonuniform-run ../../../a-thesis-exploration/1x-beam-nunif-extra
```

This assessment and the input checksums describe the data analyzed on the
date above; rerunning the plotting script regenerates figures/CSVs and the
run report but does not update this interpretation or those checksums.
