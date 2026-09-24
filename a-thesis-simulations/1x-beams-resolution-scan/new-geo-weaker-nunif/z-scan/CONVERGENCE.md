# Resolution convergence plots

Run from the postgkyl environment after the simulations finish:

```bash
python plot_convergence.py
```

The script reads numeric resolution directories beside itself and the **400-cell
baseline at `../../1x-beams`**, as in `compare.py`. It checks that `sim.c` differs
only in `Nz`, uses the resolution in the output metadata, and compares the same
frame and time in every run. Directories without frame 65 or required diagnostics
are reported and skipped. It needs at least two available resolutions. It does
not launch simulations or change simulation output.

Dependencies: the installed `postgkyl`, NumPy (2.0 or newer), and Matplotlib.
Only small moment and geometry files are loaded, not distribution functions.

Useful options:

```bash
python plot_convergence.py --no-baseline
python plot_convergence.py --baseline /path/to/reference/run
python plot_convergence.py --window-us 30 --output convergence-30us
python plot_convergence.py --frame 50 --zoom-us 15 --window-us 9
```

To include the additional 256-cell mapping:

```bash
python plot_convergence.py --nonuniform-run ../../../a-thesis-exploration/1x-beam-nunif-extra
```

This writes to `convergence-plots-with-nunif/` by default. The original scan
and its 400-cell reference are retained. Both mappings already are nonuniform;
the extra run is labeled **new mapping**, plotted as a separate point at its
actual cell count, and never connected into the original refinement sequence.
Its `sim.c` may differ only in `Nz`, `maximum_slope_at_min_B`,
`maximum_slope_at_max_B`, and `gaussian_std`. Other differences are rejected.
The CSV `mapping` column distinguishes the two 256-cell runs, and the extra
frame-history filename ends in `_nunif.csv`. An eighth figure, `grid_spacing`,
shows physical cell widths from the saved coordinate maps. See
`convergence-plots-with-nunif/ASSESSMENT.md` for the comparison produced on
2026-09-22. Passing `--output` still selects a different destination.

The defaults match these input files: final frame **65**, final time
**0.50012 s**, a **60 µs** zoom covering the last uninterrupted full-dynamics
relaxation, **15 µs** trailing averages, and **940 eV** Boltzmann electrons.
For another comparison frame, choose windows that stay within its relevant
phase. The full-history x-axis is the OAP/FDP **simulation clock**, not an
unmodified physical evolution time. Interpret convergence primarily in the
final full-dynamics interval. The finite-resolution reference is the largest
available `Nz`; it is not an exact solution.

## Outputs

`convergence-plots/` contains seven figures, each as PDF and PNG:

- `inventories_history`: integrated ion number, ion-plus-electron number, ion
  and electron kinetic energies, total kinetic energy, and central density.
- `inventories_final_relaxation`: the same quantities over the final 60 µs.
- `wall_fluxes_history`: outward particle, ion energy, and total heat flux
  densities at each wall separately.
- `wall_fluxes_final_relaxation`: the same wall quantities over the final 60 µs.
- `convergence_vs_resolution`: final values and trailing means versus `Nz` for
  integrated ion number, total kinetic energy, central density, both-wall
  particle loss rate, ion wall power, and total wall power. Bars show the
  **temporal min/max**, not statistical uncertainty.
- `relative_differences`: signed percentage differences from the finest
  available resolution, for both final values and trailing means.
- `final_profiles`: mapped physical-z profiles of ion density, flow, and
  parallel/perpendicular temperatures. Temperatures use the actual mass and
  charge in the diagnostics. Negative temperatures remain visible.

`convergence_summary.csv` includes final values, means, ranges, and relative
changes for **all** scalar diagnostics, including each wall and each species.
`frame_diagnostics_Nz*.csv` contains time histories sampled at the saved frame
times. All CSV values use unscaled diagnostic units (wall heat flux in W/m²,
not the MW/m² used in figures). `run_report.txt` records paths, code changesets,
comparison time/window, skipped runs, geometry factors, and reconstruction checks.
Repeated runs overwrite these generated files in the selected output folder.

## Definitions and normalization

These are 1D simulations of a flux tube with no specified finite transverse
extent. **Integrated numbers, energies, loss rates, and powers retain the
simulation's per-unit transverse-coordinate normalization.** They are not
whole-device particle counts, joules, or watts. Plot labels say “simulation
normalization” to preserve that distinction. Relative resolution comparisons
are unaffected by a common transverse normalization. Wall **flux densities**
are converted to physical m⁻² s⁻¹ and W/m² using the geometry.

- Ion integrated Hamiltonian diagnostics contain `[N_i, parallel momentum,
  H_i]`, where `H_i = integral (K_i + q_i n_i phi) J dz_comp`.
  Component 2 is already mass-weighted energy: do not multiply it by mass/2.
- Electron integrated diagnostics contain `[N_e, M1, M2_parallel, M2_perp]`.
  Their kinetic energy is `E_e = m_e/2 * (M2_parallel + M2_perp)`.
- Ion kinetic energy is reconstructed as
  `E_i = m_i/2 * integral J M2_i dz_comp`. The script integrates products of
  the p=1 modal coefficients exactly, including the nonuniform geometry.
  Total stored particle energy is `E_i + E_e`, including bulk-flow energy;
  total particle number is `N_i + N_e` (ion inventory is shown separately).
  The saved `H_i` is also exported as `Hi` in both CSV files.
- `field_energy` is **not** added to particle energy: for this Boltzmann-field
  configuration the code computes a phi-squared diagnostic with unit weight,
  not a physical electrostatic energy in joules. The script does not claim
  a conserved total energy for the driven, open, OAP/FDP model.
- Central density is the ion density evaluated at **physical z=0** using the
  saved coordinate map. If this is a DG interface, both traces are averaged.
- Boundary diagnostics already use **outward-positive** particle/energy loss
  at both ends. No lower-wall sign flip or extra mass factor is applied.
- Surface area per unit transverse measure is
  `A = J * |grad(z_comp)|`, read from `geo_surf0_lenr` at each wall. The surface
  file stores lower faces in its first column and an upper-boundary copy in
  its second; these columns are not p=1 volume-modal coefficients. The physical
  wall-normal particle flux is `Gamma_i = boundary_M0 / A`.
- **Total heat at the wall is a model estimate**, since kinetic-electron wall
  fluxes were not saved. With a grounded wall (`phi_wall=0`), singly charged
  ions, ambipolar losses (`Gamma_e=Gamma_i`), and escaping Maxwellian Boltzmann
  electrons, the ion arrival energy flux is `q_i = boundary_H_i / A`, the
  electron arrival energy flux is `q_e = 2 T_e Gamma_i`, and
  `q_total = q_i + q_e`. `T_e` is in joules in this expression. Ion Hamiltonian
  flux already includes acceleration through the sheath potential; do not
  add `e phi_s Gamma_i` again. This is total incident particle energy, including
  directed motion, not just the conductive third central moment. It excludes
  radiation, recombination, secondary emission, and other wall physics absent
  from this model.

The script independently reconstructs integrated ion number and Hamiltonian
energy from the spatial diagnostics and records their maximum relative
mismatch. The checks on the existing 128-cell output agree to floating-point
precision. Moment definitions and surface normalization were checked against
local Gkeyll sources, with diagnostic build changeset `02623e26d7a5`.

Time averages are trapezoidal integrals over identical time bounds. Integrated
and boundary diagnostics retain their native high time cadence for plotting
and averaging; central density and reconstructed ion/total kinetic energy use
saved frames and linear time interpolation. The default 15 µs mean therefore
uses only six spatial frames for those latter quantities. Use the reported
temporal ranges to assess settling, and avoid interpreting relative differences
alone as a measured convergence order.
