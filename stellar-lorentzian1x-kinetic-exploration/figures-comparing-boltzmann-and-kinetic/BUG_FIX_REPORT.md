# Root Finding Bug Fix Report

## Problem

The initial `rosen_dougherty_potential.py` module could not find roots in the Rosen-Dougherty confinement time equation, even though the reference implementation (`calc_pastukhov_potential.py`) worked perfectly.

### Symptoms
- Root finding failed with: `ValueError: f(a) and f(b) must have different signs`
- Root equation was **always positive** across all test points
- Expected result: e*phi/Te ≈ 7.256
- Actual behavior: Bracketing failed at [1.0, 20], [1.0, 100], [1.0, 500], [1.0, 1000]

### Root Cause: Formula Operator Precedence Bug

The Rosen-Dougherty confinement time formula in the reference code:
```python
Loss_Rosen = 1/nuElc * 1/(2*ZpFl / (log(...) - coeff) * (1-erf(...)))
```

In Python, `/` and `*` have equal precedence and evaluate **left-to-right**.

**INCORRECT interpretation (my first implementation):**
```python
loss_rosen = (1.0 / nu_e) * (
    1.0 / (
        2.0 * Zpfl / (log(...) - coeff) / (1.0 - erf(...))
    )
)
```

With left-to-right evaluation: `2*Zpfl / A / B = (2*Zpfl/A) / B = 2*Zpfl / (A*B)`

So the final formula became:
$$\text{loss\_rosen} = \frac{1}{\nu_e} \cdot \frac{\ln(...) - c \cdot (1-\text{erf}(...))}{2 \cdot Z_{pfl}}$$

The $(1-\text{erf}(...))$ term ended up in the **NUMERATOR** instead of denominator!

**CORRECT formula:**
```python
loss_rosen = (1.0 / nu_e) * (
    (np.log((w_term + 1.0) / (w_term - 1.0)) - coeff)
    / (2.0 * Zpfl * (1.0 - erf(a_term)))
)
```

Correctly placing $(1-\text{erf}(...))$ in the **DENOMINATOR**:
$$\text{loss\_rosen} = \frac{1}{\nu_e} \cdot \frac{\ln(...) - c}{2 \cdot Z_{pfl} \cdot (1-\text{erf}(...))}$$

## Impact

### Before Fix
- Root equation $\tau_{pi} - \text{RD\_time}(P)$ was always slightly positive
- No sign change in any bracket → root finding impossible
- Module was non-functional

### After Fix
- Root equation exhibits proper sign change from positive to negative
- Root finding finds root in first bracket [1.0, 20.0]
- Result matches reference implementation exactly: **7.256**

## Debug Process

1. **Added extensive debug output** to test root equation at multiple points:
   ```
   [DEBUG] root_eq(   1.000) = +1.075873e-01  ← positive
   [DEBUG] root_eq(   5.000) = +9.635083e-02  ← positive
   [DEBUG] root_eq(   7.256) = -4.369903e-05  ← SIGN CHANGE!
   [DEBUG] root_eq(  10.000) = -1.516227e+00  ← negative
   [DEBUG] root_eq(  20.000) = -2.847710e+04  ← large negative
   ```

2. **Identified the issue**: Function value at 7.256 changes from positive to negative
   - This is exactly where the root should be!
   - The problem was in the formula, not the algorithm

3. **Verified the fix**: All three implementations now agree:
   - Reference: `7.256`
   - Module (fixed): `7.256`
   - Within numerical precision: ✅

## Code Changes

**File**: `rosen_dougherty_potential.py`
**Function**: `rosen_dougherty_confinement_time()`
**Lines**: 178-194

```diff
- loss_rosen = (1.0 / nu_e) * (
-     1.0
-     / (
-         2.0 * Zpfl
-         / (np.log((w_term + 1.0) / (w_term - 1.0)) - coeff)
-         / (1.0 - erf(a_term))
-     )
- )

+ loss_rosen = (1.0 / nu_e) * (
+     (np.log((w_term + 1.0) / (w_term - 1.0)) - coeff)
+     / (2.0 * Zpfl * (1.0 - erf(a_term)))
+ )
```

## Validation

### Test Case
- Simulation: `stellar-lorentzian1x-orbit-average-beams`
- Frame: 65
- Temperature: 940 eV
- Mirror ratio: 31.39

### Results

| Implementation | Result | Status |
|---|---|---|
| Reference (`calc_pastukhov_potential.py`) | 7.256 | ✅ |
| Module (fixed) | 7.256 | ✅ Match |
| Pipeline (`master.py`) | 7.2556 | ✅ Match |

### Full Pipeline Output
```
[assess-theoretical-potential] ==========================================
[assess-theoretical-potential] Theoretical Results (Rosen-Dougherty):
[assess-theoretical-potential]   Mirror ratio R = B_throat/B_0: 31.3921
[assess-theoretical-potential]   Electron collision freq [s^-1]: 5.050e+04
[assess-theoretical-potential]   Theoretical e*phi/Te (Dougherty coeff=1.117): 7.2556
[assess-theoretical-potential] ==========================================
```

## Lessons Learned

1. **Operator precedence matters**: Python evaluates `/` and `*` left-to-right, not as a unit
2. **Always test simple cases**: Debug output at known/expected points reveals the issue quickly
3. **Validate against reference**: The exact match with `calc_pastukhov_potential.py` confirmed the fix
4. **Explicit is better than implicit**: Writing out the full formula explicitly with clear numerator/denominator made the bug obvious

## Files Updated

- ✅ `rosen_dougherty_potential.py` — Formula corrected
- ✅ `THEORETICAL_POTENTIAL_README.md` — Documentation updated with bug history
- ✅ `master.py` — Already integrated (unchanged)
- ✅ `assess-theoretical-potential-boltzmann.py` — Works correctly (unchanged)
