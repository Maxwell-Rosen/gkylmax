# Theoretical Potential Assessment (Rosen-Dougherty)

## Overview

Modular framework for computing theoretical confinement potential using the Rosen-Dougherty formulation for Boltzmann electron simulations only. Successfully validates against reference implementation (`calc_pastukhov_potential.py`).

## Status

✅ **WORKING** — Root finding algorithm successfully finds theoretical potential values matching reference implementation exactly.

## Files

### Core Modules

**`rosen_dougherty_potential.py`**
- Shared module implementing Rosen-Dougherty confinement time calculations
- Uses same postgkyl loading framework as `calc_pastukhov_potential.py` reference
- Key functions:
  - `load_magnetic_field_ratio(run_path, sim_prefix)` → Loads B_throat/B_0 ratio
  - `compute_electron_collision_freq(density, te_ev)` → Collision frequency via Coulomb formula
  - `compute_particle_confinement_time(...)` → tau_pi = ∫M0 dx / ∫source_M0 dx
  - `rosen_dougherty_confinement_time(ephi_over_te, R, nu_e, ...)` → RD confinement time
  - `compute_theoretical_potential(...)` → Solves for e*phi/Te using root finding

**`assess-theoretical-potential-boltzmann.py`**
- Assessment script for Boltzmann simulations
- Runs as part of `master.py` pipeline
- Computes and reports:
  - Mirror ratio `R = B_throat / B_0`
  - Electron collision frequency
  - Theoretical e*phi/Te 
- Handles validation against reference implementation

### Integration

**`master.py`** updated to include:
```python
PlotTask(script="assess-theoretical-potential-boltzmann.py", enabled=True)
```

## Physics

### Formulation

The script solves for the normalized potential $P = e\phi/T_e$ using:

$$\tau_{pi} = \frac{1}{\nu_e} \cdot \frac{\ln\left(\frac{w+1}{w-1}\right) - c}{2Z_{pfl} \cdot (1 - \text{erf}(a))}$$

Where:
- $w = \sqrt{1 + \frac{2P}{RZ_{pfl}}}$ (Rosen-Dougherty parameter)
- $a = \sqrt{P + \ln(w)}$
- $\tau_{pi}$ = particle confinement time from simulation integration
- $\nu_e$ = electron-electron collision frequency (Coulomb, from Najmabadi 1984)
- $c = 1.117$ (Dougherty coefficient)
- $R = B_{throat}/ B_0$ (mirror ratio)

### Data Sources (Boltzmann)

- **Magnetic field**: `{sim_name}-bmag.gkyl` (loads R at z=0 and z=0.98)
- **Electron density**: `{sim_name}-ion_M0_{frame}.gkyl` (quasineutrality: n_e = n_i)
- **Confinement time**: Computed from ion moments integrals
  - $\int_{-0.98}^{0.98} M0 \, dx$ (trapped particles)
  - $\int_{-0.98}^{0.98} \text{source\_M0} \, dx$ (source injection)

## Usage

Run full pipeline:
```bash
python master.py
```

This executes all four tasks:
1. `potential.py` - Normalized potential plot
2. `electric-field.py` - Physical E_z field (V/m)
3. `electric-field-normalized.py` - Normalized E_z = -∇(eφ/Te)
4. `assess-theoretical-potential-boltzmann.py` - Theoretical assessment

## Validation Results

### Test Case: stellar-lorentzian1x-orbit-average-beams (frame 65)

**Reference Script Output** (`calc_pastukhov_potential.py`):
```
Dougherty e*phi/Te = 7.256
```

**Module Output** (`assess-theoretical-potential-boltzmann.py`):
```
[assess-theoretical-potential] Theoretical e*phi/Te (Dougherty coeff=1.117): 7.2556
```

✅ **Match**: Within numerical precision (7.256 vs 7.2556)

### Console Output Example

```
[assess-theoretical-potential] Computing Rosen-Dougherty theoretical potential...
[assess-theoretical-potential] Boltzmann run: /path/to/stellar-lorentzian1x-orbit-average-beams
[assess-theoretical-potential] Frame: 65, Te(z=0): 940.0 eV
[assess-theoretical-potential] Detected sim_prefix: gk_lorentzian_mirror
[assess-theoretical-potential] ==========================================
[assess-theoretical-potential] Theoretical Results (Rosen-Dougherty):
[assess-theoretical-potential]   Mirror ratio R = B_throat/B_0: 31.3921
[assess-theoretical-potential]   Electron collision freq [s^-1]: 5.050e+04
[assess-theoretical-potential]   Theoretical e*phi/Te (Dougherty coeff=1.117): 7.2556
[assess-theoretical-potential] ==========================================
```

## Implementation Details

### Root Finding Algorithm

1. **Primary method**: `scipy.optimize.ridder` with initially bracketed region `[zpfl=1.0, 20.0]`
2. **Adaptive bracketing**: If initial bracket fails, expands to `[1.0, 100.0]`, `[1.0, 500.0]`, `[1.0, 1000.0]`
3. **Fallback method**: `scipy.optimize.brent` if all ridder attempts fail

### Bug Fix History

**Issue**: Initial implementation had incorrect operator precedence in formula simplification
- Error: `(1-erf(...))` term was placed in numerator instead of denominator
- Symptom: Root equation always positive, ridder method couldn't find bracket with sign change
- Fix: Corrected formula to match reference exactly:
  ```python
  # WRONG: loss_rosen = (1/nu_e) * 1 / (2*Zpfl / (log(...) - coeff) / (1 - erf(...)))
  # RIGHT: loss_rosen = (1/nu_e) * (log(...) - coeff) / (2*Zpfl * (1 - erf(...)))
  ```
- Validation: Result now matches reference implementation exactly (7.256)

### Postgkyl API Pattern

```python
bmag_data = pg.GData(file_path, mapc2p_name=map_path)
pg.GInterpModal(bmag_data).interpolate(comp, overwrite=True)
_, value = pg.data.select(bmag_data, z0=z_value)
```

## Files and Sizes

### Generated Output
- All three PDF plots generated successfully (~700 KB total):
  - `potential-overlay-kinetic-vs-boltzmann.pdf` (295 KB)
  - `electric-field-vm-overlay-kinetic-vs-boltzmann.pdf` (168 KB)
  - `electric-field-normalized-overlay-kinetic-vs-boltzmann.pdf` (229 KB)

## Notes

- **Modularity**: RD calculations packaged separately for reuse in other analysis scripts
- **Quasineutrality**: For Boltzmann electrons, electron density n_e is inferred from ion density n_i
- **Physical validity**: Results are physical only when parameters are within model validity range
- **Extensibility**: Coefficients (zpfl, coeff) parameterized for different collision assumptions
