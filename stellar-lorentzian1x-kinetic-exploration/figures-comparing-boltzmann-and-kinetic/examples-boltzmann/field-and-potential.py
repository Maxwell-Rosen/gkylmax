"""Plot potential and electric-field overlays for Boltzmann vs kinetic runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import postgkyl as pg


plt.rcParams.update(
	{
		"text.usetex": True,
		"font.family": "serif",
		"font.size": 12,
		"axes.titlesize": 14,
		"axes.labelsize": 14,
		"legend.fontsize": 11,
		"xtick.labelsize": 12,
		"ytick.labelsize": 12,
	}
)


KINETIC_COLOR = "tab:blue"
BOLTZMANN_COLOR = "tab:red"

ELECTRON_MASS = 9.1093837015e-31
ELEMENTARY_CHARGE = 1.602176634e-19

Z_CENTER = 0.0
Z_MIRROR_THROAT = 0.98
Z_SHEATH = 2.5


@dataclass(frozen=True)
class RunSpec:
	label: str
	data_path: Path
	last_frame: int
	color: str


def _find_field_file(run: RunSpec) -> Path:
	pattern = f"*-field_{run.last_frame}.gkyl"
	candidates = sorted(run.data_path.glob(pattern))
	if not candidates:
		raise FileNotFoundError(
			f"No field file found for {run.label} with pattern {pattern} in {run.data_path}"
		)
	return candidates[0]


def _sim_prefix_from_field(field_file: Path) -> str:
	marker = "-field_"
	stem = field_file.stem
	if marker not in stem:
		raise ValueError(f"Could not infer simulation prefix from {field_file.name}")
	return stem.split(marker, 1)[0]


def _find_map_file(run: RunSpec, sim_prefix: str) -> Path | None:
	candidates = [
		run.data_path / f"{sim_prefix}-mc2nu_pos.gkyl",
		run.data_path / f"{sim_prefix}-mc2nu_pos_deflated.gkyl",
		run.data_path / f"{sim_prefix}-mapc2p.gkyl",
	]
	for path in candidates:
		if path.exists():
			return path
	return None


def _candidate_map_components(map_path: Path) -> tuple[int, ...]:
	name = map_path.name
	if "mc2nu_pos_deflated" in name:
		return (0,)
	if "mc2nu_pos" in name:
		return (2, 0, 1)
	if "mapc2p" in name:
		return (1, 0, 2)
	return (0, 1, 2)


def _load_z_from_map(map_path: Path) -> np.ndarray | None:
	map_data = pg.GData(str(map_path))
	map_interp = pg.GInterpModal(map_data)

	best: np.ndarray | None = None
	best_span = -np.inf

	for comp in _candidate_map_components(map_path):
		try:
			_, z_raw = map_interp.interpolate(comp=comp)
		except Exception:
			continue

		z_try = _to_1d(z_raw)
		if z_try.size < 4:
			continue

		z_span = float(np.max(z_try) - np.min(z_try))
		if z_span <= 1e-8:
			continue

		unique_count = np.unique(np.round(z_try, 12)).size
		if unique_count < max(8, int(0.2 * z_try.size)):
			continue

		if z_span > best_span:
			best = z_try
			best_span = z_span

	return best


def _to_1d(arr: Any) -> np.ndarray:
	out = np.asarray(arr)
	out = np.squeeze(out)
	if out.ndim == 0:
		out = out.reshape(1)
	return out


def _load_profile(run: RunSpec) -> tuple[np.ndarray, np.ndarray]:
	field_file = _find_field_file(run)
	sim_prefix = _sim_prefix_from_field(field_file)
	map_file = _find_map_file(run, sim_prefix)

	field_data = pg.GData(str(field_file))
	field_interp = pg.GInterpModal(field_data)
	x_coords, phi_raw = field_interp.interpolate(comp=0)

	phi = _to_1d(phi_raw)

	z_from_map: np.ndarray | None = None
	if map_file is not None:
		z_from_map = _load_z_from_map(map_file)

	if z_from_map is not None:
		z = z_from_map
	else:
		# Fallback for runs without mapc2p-like files.
		z = _to_1d(x_coords[0])

	size = min(z.size, phi.size)
	z = z[:size]
	phi = phi[:size]

	sort_idx = np.argsort(z)
	z_sorted = z[sort_idx]
	phi_sorted = phi[sort_idx]

	# Enforce unique coordinates before gradient calculations.
	z_unique, inverse = np.unique(np.round(z_sorted, 12), return_inverse=True)
	phi_unique = np.zeros_like(z_unique, dtype=float)
	counts = np.zeros_like(z_unique, dtype=float)
	for i, idx in enumerate(inverse):
		phi_unique[idx] += float(phi_sorted[i])
		counts[idx] += 1.0
	phi_unique /= np.maximum(counts, 1.0)

	return z_unique, phi_unique


def _load_kinetic_te0_ev(run: RunSpec, sim_prefix: str, map_file: Path | None) -> float:
	"""Load kinetic electron temperature at z=0 from BiMaxwellian moments."""
	moments_file = run.data_path / f"{sim_prefix}-elc_BiMaxwellianMoments_{run.last_frame}.gkyl"
	if not moments_file.exists():
		raise FileNotFoundError(
			f"Missing kinetic moments file for Te(z=0): {moments_file}"
		)

	# Keep moments loading identical to validated analysis scripts:
	# GData(path) + GInterpModal(..., 1, "ms") with no mapc2p_name.
	moments_data = pg.GData(str(moments_file))
	interp = pg.GInterpModal(moments_data, 1, "ms")

	x_coords, tpar_raw = interp.interpolate(2)
	_, tperp_raw = interp.interpolate(3)

	tpar = _to_1d(tpar_raw)
	tperp = _to_1d(tperp_raw)

	z_from_map: np.ndarray | None = None
	if map_file is not None:
		z_from_map = _load_z_from_map(map_file)

	if z_from_map is not None:
		z = z_from_map
	else:
		z = _to_1d(x_coords[0])

	size = min(z.size, tpar.size, tperp.size)
	z = z[:size]
	tpar = tpar[:size]
	tperp = tperp[:size]

	sort_idx = np.argsort(z)
	z_sorted = z[sort_idx]
	tpar_sorted = tpar[sort_idx]
	tperp_sorted = tperp[sort_idx]

	z_unique, inverse = np.unique(np.round(z_sorted, 12), return_inverse=True)
	tpar_unique = np.zeros_like(z_unique, dtype=float)
	tperp_unique = np.zeros_like(z_unique, dtype=float)
	counts = np.zeros_like(z_unique, dtype=float)
	for i, idx in enumerate(inverse):
		tpar_unique[idx] += float(tpar_sorted[i])
		tperp_unique[idx] += float(tperp_sorted[i])
		counts[idx] += 1.0
	tpar_unique /= np.maximum(counts, 1.0)
	tperp_unique /= np.maximum(counts, 1.0)

	tpar_j = float(np.interp(Z_CENTER, z_unique, tpar_unique)) * ELECTRON_MASS
	tperp_j = float(np.interp(Z_CENTER, z_unique, tperp_unique)) * ELECTRON_MASS
	te0_j = (tpar_j + 2.0 * tperp_j) / 3.0
	if te0_j <= 0.0:
		raise ValueError(f"Non-positive kinetic Te(z=0): {te0_j}")

	return te0_j / ELEMENTARY_CHARGE


def _delta_labels(z: np.ndarray, phi_norm: np.ndarray, prefix: str) -> tuple[str, float, float]:
	phi_center = float(np.interp(Z_CENTER, z, phi_norm))
	phi_throat = float(np.interp(Z_MIRROR_THROAT, z, phi_norm))
	phi_sheath = float(np.interp(Z_SHEATH, z, phi_norm))

	d_center_to_throat = phi_center - phi_throat
	d_throat_to_sheath = phi_throat - phi_sheath

	label = "\n".join(
		[
			f"{prefix}:",
			f"$\\Delta (e\\phi/T_e)_{{0\\to0.98}} = {d_center_to_throat:.3f}$",
			f"$\\Delta (e\\phi/T_e)_{{0.98\\to2.5}} = {d_throat_to_sheath:.3f}$",
		]
	)
	return label, d_center_to_throat, d_throat_to_sheath


def _plot_potential(
	out_dir: Path,
	z_kin: np.ndarray,
	phi_kin_norm: np.ndarray,
	z_boltz: np.ndarray,
	phi_boltz_norm: np.ndarray,
) -> None:
	kin_label, _, _ = _delta_labels(z_kin, phi_kin_norm, "Kinetic")
	boltz_label, _, _ = _delta_labels(z_boltz, phi_boltz_norm, "Boltzmann")

	fig, ax = plt.subplots(figsize=(7.0, 4.4))
	ax.plot(z_kin, phi_kin_norm, color=KINETIC_COLOR, lw=2.0, label="Kinetic electrons")
	ax.plot(z_boltz, phi_boltz_norm, color=BOLTZMANN_COLOR, lw=2.0, label="Boltzmann electrons")

	for marker_x, marker_label in (
		(Z_CENTER, "center"),
		(Z_MIRROR_THROAT, "mirror throat"),
		(Z_SHEATH, "sheath"),
	):
		ax.axvline(marker_x, color="0.35", lw=1.0, ls="--")
		ax.text(marker_x, 0.98, marker_label, transform=ax.get_xaxis_transform(),
				ha="center", va="top", fontsize=9)

	ax.set_xlabel(r"$z$ [m]")
	ax.set_ylabel(r"$e\phi / T_e$")
	ax.set_title("Potential profile at last frame (normalized)")
	ax.grid(alpha=0.3)
	ax.legend(loc="best")

	annotation = kin_label + "\n\n" + boltz_label
	ax.text(
		0.02,
		0.02,
		annotation,
		transform=ax.transAxes,
		ha="left",
		va="bottom",
		fontsize=10.5,
		linespacing=1.25,
		bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.95},
	)

	fig.tight_layout()
	fig.savefig(out_dir / "potential-overlay-kinetic-vs-boltzmann.pdf")
	plt.close(fig)


def _plot_electric_field(
	out_dir: Path,
	z_kin: np.ndarray,
	phi_kin_norm: np.ndarray,
	z_boltz: np.ndarray,
	phi_boltz_norm: np.ndarray,
) -> None:
	e_kin = -np.gradient(phi_kin_norm, z_kin)
	e_boltz = -np.gradient(phi_boltz_norm, z_boltz)

	fig, ax = plt.subplots(figsize=(7.0, 4.4))
	ax.plot(z_kin, e_kin, color=KINETIC_COLOR, lw=2.0, label="Kinetic electrons")
	ax.plot(z_boltz, e_boltz, color=BOLTZMANN_COLOR, lw=2.0, label="Boltzmann electrons")

	for marker_x in (Z_CENTER, Z_MIRROR_THROAT, Z_SHEATH):
		ax.axvline(marker_x, color="0.35", lw=1.0, ls="--")

	ax.set_xlabel(r"$z$ [m]")
	ax.set_ylabel(r"$-\nabla (e\phi/T_e)$ [m$^{-1}$]")
	ax.set_title("Electric field proxy from normalized potential")
	ax.grid(alpha=0.3)
	ax.legend(loc="best")

	fig.tight_layout()
	fig.savefig(out_dir / "electric-field-overlay-kinetic-vs-boltzmann.pdf")
	plt.close(fig)


def run(context: dict[str, Any]) -> None:
	out_dir = Path(context["figures_dir"]) / "field-and-potential"
	out_dir.mkdir(parents=True, exist_ok=True)

	kinetic_cfg = context["kinetic"]
	boltzmann_cfg = context["boltzmann"]

	kinetic = RunSpec(
		label=kinetic_cfg["label"],
		data_path=Path(kinetic_cfg["data_path"]),
		last_frame=int(kinetic_cfg["last_frame"]),
		color=KINETIC_COLOR,
	)
	boltzmann = RunSpec(
		label=boltzmann_cfg["label"],
		data_path=Path(boltzmann_cfg["data_path"]),
		last_frame=int(boltzmann_cfg["last_frame"]),
		color=BOLTZMANN_COLOR,
	)

	z_kin, phi_kin = _load_profile(kinetic)
	z_boltz, phi_boltz = _load_profile(boltzmann)

	kin_field_file = _find_field_file(kinetic)
	kin_prefix = _sim_prefix_from_field(kin_field_file)
	kin_map_file = _find_map_file(kinetic, kin_prefix)

	kin_te0_ev = _load_kinetic_te0_ev(kinetic, kin_prefix, kin_map_file)
	boltz_te0_ev = float(context["constants"]["boltzmann_te0_ev"])

	phi_kin_norm = phi_kin / kin_te0_ev
	phi_boltz_norm = phi_boltz / boltz_te0_ev

	_plot_potential(out_dir, z_kin, phi_kin_norm, z_boltz, phi_boltz_norm)
	_plot_electric_field(out_dir, z_kin, phi_kin_norm, z_boltz, phi_boltz_norm)

	kin_label, kin_d1, kin_d2 = _delta_labels(z_kin, phi_kin_norm, "Kinetic")
	boltz_label, boltz_d1, boltz_d2 = _delta_labels(z_boltz, phi_boltz_norm, "Boltzmann")
	print("[field-and-potential] Saved figures to", out_dir)
	print(
		"[field-and-potential] Te used [eV]: "
		f"kinetic(z=0)={kin_te0_ev:.3f}, boltzmann(master)={boltz_te0_ev:.3f}"
	)
	print("[field-and-potential]", kin_label)
	print("[field-and-potential]", boltz_label)
	print(
		"[field-and-potential] Delta summary "
		f"(center->throat, throat->sheath) [dimensionless]: "
		f"kinetic=({kin_d1:.3f}, {kin_d2:.3f}), "
		f"boltzmann=({boltz_d1:.3f}, {boltz_d2:.3f})"
	)


def main(context: dict[str, Any] | None = None) -> None:
	if context is None:
		raise ValueError(
			"This subfile is intended to be called by master.py with context."
		)
	run(context)

