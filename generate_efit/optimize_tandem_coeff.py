import numpy as np
import scipy.optimize as sco
import matplotlib.pyplot as plt

L_center = 4.0
L_plugs = 2.0
L_to_wall = 0.5

Z_INNER = L_center / 2.0
Z_OUTER = L_center / 2.0 + L_plugs
Z_END_MIN = L_center / 2.0 + L_plugs / 2.0
Z_WALL = Z_OUTER + L_to_wall

B_END_TARGET = 0.5



def B_axis(z, mcB_inner, mcB_outer, gamma_inner, gamma_outer):
  # Axis field from psi(R, z) = 0.5 * R^2 * B(z).
  inner = mcB_inner * (
      1.0 / (np.pi * gamma_inner * (1.0 + ((z - Z_INNER) / gamma_inner) ** 2))
      + 1.0 / (np.pi * gamma_inner * (1.0 + ((z + Z_INNER) / gamma_inner) ** 2))
  )
  outer = mcB_outer * (
      1.0 / (np.pi * gamma_outer * (1.0 + ((z - Z_OUTER) / gamma_outer) ** 2))
      + 1.0 / (np.pi * gamma_outer * (1.0 + ((z + Z_OUTER) / gamma_outer) ** 2))
  )
  return inner + outer

R_values     = [3, 5, 10, 15, 22, 32, 40, 50]


def residuals(params, desired_R, ref_params):
  mcB_inner, mcB_outer, gamma_inner, gamma_outer = params

  # Evaluate at positive-z characteristic points (symmetry gives matching -z values).
  b_inner = B_axis(Z_INNER, mcB_inner, mcB_outer, gamma_inner, gamma_outer)
  b_outer = B_axis(Z_OUTER, mcB_inner, mcB_outer, gamma_inner, gamma_outer)
  b_end = B_axis(Z_END_MIN, mcB_inner, mcB_outer, gamma_inner, gamma_outer)

  b_peak_target = desired_R * B_END_TARGET

  eq_equal_peaks = b_inner - b_outer
  eq_peak_level = 0.5 * (b_inner + b_outer) - b_peak_target
  eq_end = b_end - B_END_TARGET
  eq_ratio = (0.5 * (b_inner + b_outer) / b_end) - desired_R

  # Mild regularization to keep widths from drifting to extremes.
  reg = 1e-3 * (np.array([gamma_inner, gamma_outer]) - np.array([ref_params[2], ref_params[3]]))

  return np.array([eq_equal_peaks, eq_peak_level, eq_end, eq_ratio, reg[0], reg[1]])


def optimize_for_ratio(desired_R, initial_guess):
  bounds_lower = np.array([1e-6, 1e-6, 1e-3, 1e-3])
  bounds_upper = np.array([50.0, 50.0, 5.0, 5.0])

  result = sco.least_squares(
      residuals,
      x0=np.array(initial_guess, dtype=float),
      bounds=(bounds_lower, bounds_upper),
      args=(desired_R, initial_guess),
      xtol=1e-12,
      ftol=1e-12,
      gtol=1e-12,
      max_nfev=20000,
  )

  if not result.success:
    raise RuntimeError(f"Optimization failed for R={desired_R}: {result.message}")

  return result.x


def summarize_solution(desired_R, params):
  mcB_inner, mcB_outer, gamma_inner, gamma_outer = params
  b_inner = B_axis(Z_INNER, mcB_inner, mcB_outer, gamma_inner, gamma_outer)
  b_outer = B_axis(Z_OUTER, mcB_inner, mcB_outer, gamma_inner, gamma_outer)
  b_end = B_axis(Z_END_MIN, mcB_inner, mcB_outer, gamma_inner, gamma_outer)
  b_max = max(b_inner, b_outer)
  r_actual = b_max / b_end

  return {
      "R_target": desired_R,
      "mcB_inner": mcB_inner,
      "mcB_outer": mcB_outer,
      "gamma_inner": gamma_inner,
      "gamma_outer": gamma_outer,
      "B_inner_peak": b_inner,
      "B_outer_peak": b_outer,
      "B_end": b_end,
      "B_max": b_max,
      "R_actual": r_actual,
      "peak_mismatch": abs(b_inner - b_outer),
      "end_error": abs(b_end - B_END_TARGET),
      "ratio_error": abs(r_actual - desired_R),
  }


if __name__ == "__main__":
  # Seed from a moderate-R hand-tuned profile and continue from previous solution.
  guess = [2.7, 2.7, 0.35, 0.35]
  results = []

  print("--- Tandem Mirror Coefficient Optimization ---")
  print(
      "Target constraints: equal peak B, B_end = "
      f"{B_END_TARGET:.3f} T at z=+-{Z_END_MIN:.3f}, mirror ratio R = Bmax/B_end"
  )
  print(
      f"{'R':>6s}  {'mcB_i':>10s}  {'mcB_o':>10s}  {'g_i':>9s}  {'g_o':>9s}  "
      f"{'Bpk_i':>10s}  {'Bpk_o':>10s}  {'Bend':>10s}  {'R_act':>10s}"
  )
  print("-" * 106)

  for desired_R in R_values:
    params = optimize_for_ratio(desired_R, guess)
    summary = summarize_solution(desired_R, params)
    results.append(summary)

    print(
        f"{desired_R:6.1f}  {summary['mcB_inner']:10.6f}  {summary['mcB_outer']:10.6f}  "
        f"{summary['gamma_inner']:9.6f}  {summary['gamma_outer']:9.6f}  "
        f"{summary['B_inner_peak']:10.6f}  {summary['B_outer_peak']:10.6f}  "
        f"{summary['B_end']:10.6f}  {summary['R_actual']:10.6f}"
    )

    guess = params

  max_peak_mismatch = max(r["peak_mismatch"] for r in results)
  max_end_error = max(r["end_error"] for r in results)
  max_ratio_error = max(r["ratio_error"] for r in results)

  print("\n--- Constraint Check (max abs error across all R) ---")
  print(f"peak mismatch |B_inner - B_outer| : {max_peak_mismatch:.3e}")
  print(f"end-cell error |B_end - 0.5|      : {max_end_error:.3e}")
  print(f"ratio error |R_actual - R_target| : {max_ratio_error:.3e}")

  # Plot all optimized profiles.
  z = np.linspace(-Z_WALL, Z_WALL, 2001)
  plt.figure(figsize=(11, 7))
  colors = plt.cm.Set2(np.linspace(0.0, 1.0, len(results)))
  for i, res in enumerate(results):
    b = B_axis(z, res["mcB_inner"], res["mcB_outer"], res["gamma_inner"], res["gamma_outer"])
    plt.plot(z, b, lw=2.0, color=colors[i], label=f"R={res['R_target']:.0f}")

  plt.axhline(B_END_TARGET, color="k", ls="--", lw=1.0, label=f"B_end target = {B_END_TARGET:.2f} T")
  plt.axvline(Z_INNER, color="gray", ls=":", lw=1.0)
  plt.axvline(-Z_INNER, color="gray", ls=":", lw=1.0)
  plt.axvline(Z_OUTER, color="gray", ls=":", lw=1.0)
  plt.axvline(-Z_OUTER, color="gray", ls=":", lw=1.0)
  plt.xlabel("z")
  plt.ylabel("B(z) [T]")
  plt.title("Optimized Tandem Mirror B(z) Profiles")
  plt.grid(alpha=0.2)
  plt.legend(ncol=2, fontsize=9)
  plt.tight_layout()
  plt.show()

  # Copy-paste friendly arrays.
  print("\n--- Arrays for copy-paste ---")
  print("R_values = [" + ", ".join(str(int(r["R_target"])) for r in results) + "]")
  print("mcB_inner_values = [" + ", ".join(f"{r['mcB_inner']:.6f}" for r in results) + "]")
  print("mcB_outer_values = [" + ", ".join(f"{r['mcB_outer']:.6f}" for r in results) + "]")
  print("gamma_inner_values = [" + ", ".join(f"{r['gamma_inner']:.6f}" for r in results) + "]")
  print("gamma_outer_values = [" + ", ".join(f"{r['gamma_outer']:.6f}" for r in results) + "]")
