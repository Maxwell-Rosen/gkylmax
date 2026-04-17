"""Distribution-function comparison figure (modular, configurable).

Layout defaults to 3 rows:
1) Boltzmann ions
2) Kinetic ions
3) Boltzmann electrons

Columns are configurable z-slices. By default, includes z=0.8 m in addition to
the original locations used in the legacy plotting script.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import numpy as np
import postgkyl as pg


plt.rcParams.update(
	{
		"text.usetex": True,
		"font.family": "serif",
		"font.size": 12,
		"axes.titlesize": 16,
		"axes.labelsize": 20,
		"legend.fontsize": 11,
		"xtick.labelsize": 12,
		"ytick.labelsize": 12,
	}
)


@dataclass(frozen=True)
class SpeciesConfig:
	name: str
	mass: float
	charge: float
	temp_joule: float
	label: str


@dataclass(frozen=True)
class RunSpeciesConfig:
	run_key: str
	species: SpeciesConfig
	row_label: str


@dataclass
class DistributionData:
	z: np.ndarray
	vpar_edges: np.ndarray
	mu_edges: np.ndarray
	f: np.ndarray
	bmag: np.ndarray
	phi: np.ndarray
	jb: np.ndarray


E_CHARGE = 1.602176634e-19
MASS_PROTON = 1.67262192369e-27
MASS_ELECTRON = 9.1093837015e-31

MASS_ION = 2.014 * MASS_PROTON
CHARGE_ION = E_CHARGE
CHARGE_ELECTRON = -E_CHARGE

TE0_J = 940.0 * E_CHARGE
N0 = 3.0e19
MU0 = 4.0 * np.pi * 1e-7
B0_REF = 0.53
BETA = 0.4
TAU = (B0_REF**2) * BETA / (2.0 * MU0 * N0 * TE0_J) - 1.0
TI0_J = TAU * TE0_J


ION_CFG = SpeciesConfig(
	name="ion",
	mass=MASS_ION,
	charge=CHARGE_ION,
	temp_joule=TI0_J,
	label="Ion",
)
ELC_CFG = SpeciesConfig(
	name="elc",
	mass=MASS_ELECTRON,
	charge=CHARGE_ELECTRON,
	temp_joule=TE0_J,
	label="Electron",
)


# Edit this list to change/add rows.
ROW_SPECS: list[RunSpeciesConfig] = [
	RunSpeciesConfig("boltzmann", ION_CFG, r"Boltzmann ions"),
	RunSpeciesConfig("kinetic", ION_CFG, r"Kinetic ions"),
	RunSpeciesConfig("kinetic", ELC_CFG, r"Kinetic electrons"),
]

# Edit this list to change/add z locations.
PLOT_Z_LOCATIONS = [0.0, 0.9, 0.98, 2.5]


def _to_1d(arr: Any) -> np.ndarray:
	out = np.asarray(arr)
	out = np.squeeze(out)
	if out.ndim == 0:
		out = out.reshape(1)
	return out


def _sim_prefix(run_path: Path, frame: int) -> str:
	files = sorted(run_path.glob(f"*-field_{frame}.gkyl"))
	if not files:
		raise FileNotFoundError(f"No field file found in {run_path} for frame {frame}")
	return files[0].stem.split("-field_", 1)[0]


def _load_z_nodes(run_path: Path, prefix: str) -> np.ndarray:
	map_path = run_path / f"{prefix}-mapc2p.gkyl"
	if not map_path.exists():
		raise FileNotFoundError(f"Missing map file: {map_path}")
	data = pg.GData(str(map_path))
	interp = pg.GInterpModal(data, 1, "ms")
	_, z_nodes = interp.interpolate(1)
	return _to_1d(z_nodes)


def _load_distribution_data(run_path: Path, prefix: str, frame: int, species: SpeciesConfig) -> DistributionData:
	f_path = run_path / f"{prefix}-{species.name}_{frame}.gkyl"
	map_vel_path = run_path / f"{prefix}-{species.name}_mapc2p_vel.gkyl"
	jv_path = run_path / f"{prefix}-{species.name}_jacobvel.gkyl"
	bmag_path = run_path / f"{prefix}-bmag.gkyl"
	phi_path = run_path / f"{prefix}-field_{frame}.gkyl"
	jb_path = run_path / f"{prefix}-jacobtot.gkyl"

	for required in (f_path, map_vel_path, jv_path, bmag_path, phi_path, jb_path):
		if not required.exists():
			raise FileNotFoundError(f"Missing required file: {required}")

	f_data = pg.GData(str(f_path), mapc2p_vel_name=str(map_vel_path))
	jv_data = pg.GData(str(jv_path))
	f_data._values = f_data.get_values() / jv_data.get_values()
	f_interp = pg.GInterpModal(f_data, 1, "gkhyb")

	bmag_interp = pg.GInterpModal(pg.GData(str(bmag_path)), 1, "ms")
	phi_interp = pg.GInterpModal(pg.GData(str(phi_path)), 1, "ms")
	jb_interp = pg.GInterpModal(pg.GData(str(jb_path)), 1, "ms")

	x_grid, f_raw = f_interp.interpolate()
	_, bmag_raw = bmag_interp.interpolate(0)
	_, phi_raw = phi_interp.interpolate(0)
	_, jb_raw = jb_interp.interpolate(0)

	z = _load_z_nodes(run_path, prefix)
	f = np.abs(np.squeeze(f_raw))
	f[f < 1e-16] = 1e-16

	vpar_edges = _to_1d(x_grid[1])
	mu_edges = _to_1d(x_grid[2])

	return DistributionData(
		z=z,
		vpar_edges=vpar_edges,
		mu_edges=mu_edges,
		f=f,
		bmag=_to_1d(bmag_raw),
		phi=_to_1d(phi_raw),
		jb=_to_1d(jb_raw),
	)


def _z_indices(z_nodes: np.ndarray, z_locations: list[float]) -> list[int]:
	return [int(np.argmin(np.abs(z_nodes - z0))) for z0 in z_locations]


def _normalize_axes(data: DistributionData, species: SpeciesConfig) -> tuple[np.ndarray, np.ndarray, float]:
	vth = np.sqrt(species.temp_joule / species.mass)
	vpar_norm = data.vpar_edges / vth

	b_center = float(np.interp(0.0, data.z, data.bmag))
	mu_ref = 0.5 * species.mass * vth**2 / b_center
	mu_norm = data.mu_edges / mu_ref
	return vpar_norm, mu_norm, mu_ref


def _loss_cone_curve(
	data: DistributionData,
	species: SpeciesConfig,
	z_eval: float,
	z_throat: float = 0.98,
) -> tuple[np.ndarray, np.ndarray]:
	b_eval = float(np.interp(z_eval, data.z, data.bmag))
	b_throat = float(np.interp(z_throat, data.z, data.bmag))
	phi_eval = float(np.interp(z_eval, data.z, data.phi))
	phi_throat = float(np.interp(z_throat, data.z, data.phi))

	delta_b = b_throat - b_eval
	if np.isclose(delta_b, 0.0):
		# At/near the throat the classical loss-cone denominator vanishes.
		vpar_norm, _, _ = _normalize_axes(data, species)
		return vpar_norm, np.full_like(vpar_norm, np.nan, dtype=float)

	vpar = data.vpar_edges
	mu_loss = (
		0.5 * species.mass * vpar**2
		+ species.charge * (phi_eval - phi_throat)
	) / delta_b

	vpar_norm, _, mu_ref = _normalize_axes(data, species)
	mu_loss_norm = mu_loss / mu_ref
	return vpar_norm, mu_loss_norm


def _panel_label(index: int) -> str:
	return rf"$\mathbf{{({chr(ord('a') + index)})}}$"


def _add_zoom_inset(
	ax_main,
	x_nodal: np.ndarray,
	y_nodal: np.ndarray,
	panel: np.ndarray,
	norm,
	*,
	xlim: tuple[float, float],
	ylim: tuple[float, float],
	xticks: list[float],
	yticks: list[float],
) -> None:
	"""Add legacy-style inset used for throat/exhaust panels."""
	ax_ins = inset_axes(ax_main, width="60%", height="42%", loc="upper center", borderpad=1.0)
	ax_ins.pcolormesh(
		x_nodal,
		y_nodal,
		panel,
		cmap="inferno",
		norm=norm,
		rasterized=True,
	)
	ax_ins.set_xticks(xticks)
	ax_ins.set_yticks(yticks)
	ax_ins.set_xlim(*xlim)
	ax_ins.set_ylim(*ylim)
	ax_ins.tick_params(labelsize=8)
	for spine in ax_ins.spines.values():
		spine.set_edgecolor("0.6")
		spine.set_linewidth(1.1)
	mark_inset(ax_main, ax_ins, loc1=3, loc2=4, fc="none", ec="white", lw=0.8)


def _add_exhaust_inset(
	ax_main,
	x_nodal: np.ndarray,
	y_nodal: np.ndarray,
	panel: np.ndarray,
	norm,
	*,
	xlim: tuple[float, float],
	ylim: tuple[float, float],
	xticks: list[float],
	yticks: list[float],
) -> None:
	"""Add legacy-style inset for far-exhaust panels."""
	ax_ins = inset_axes(ax_main, width="60%", height="42%", loc="upper center", borderpad=1.0)
	ax_ins.pcolormesh(
		x_nodal,
		y_nodal,
		panel,
		cmap="inferno",
		norm=norm,
		rasterized=True,
	)
	ax_ins.set_xticks(xticks)
	ax_ins.set_yticks(yticks)
	ax_ins.set_xlim(*xlim)
	ax_ins.set_ylim(*ylim)
	ax_ins.tick_params(labelsize=8)
	for spine in ax_ins.spines.values():
		spine.set_edgecolor("0.6")
		spine.set_linewidth(1.1)
	mark_inset(ax_main, ax_ins, loc1=3, loc2=4, fc="none", ec="white", lw=0.8)


def run(context: dict[str, Any]) -> None:
	out_dir = Path(context["figures_dir"]) / "distributions"
	out_dir.mkdir(parents=True, exist_ok=True)

	run_cfg = {
		"kinetic": context["kinetic"],
		"boltzmann": context["boltzmann"],
	}

	loaded: dict[tuple[str, str], DistributionData] = {}
	for row in ROW_SPECS:
		key = (row.run_key, row.species.name)
		if key in loaded:
			continue
		cfg = run_cfg[row.run_key]
		run_path = Path(cfg["data_path"])
		frame = int(cfg["last_frame"])
		prefix = _sim_prefix(run_path, frame)
		loaded[key] = _load_distribution_data(run_path, prefix, frame, row.species)

	n_rows = len(ROW_SPECS)
	n_cols = len(PLOT_Z_LOCATIONS)
	plot_z_locations = PLOT_Z_LOCATIONS
	fig, ax = plt.subplots(n_rows, n_cols, figsize=(3.2 * n_cols, 2.5 * n_rows), squeeze=False)

	# Separate inset extents for ions (top rows) and electrons (bottom row).
	ion_throat_inset = {
		"xlim": (-0.5, 2.5),
		"ylim": (0.0, 0.5),
		"xticks": [-0.5, 2.5],
		"yticks": [0.0, 0.5],
	}
	ion_exhaust_inset = {
		"xlim": (0.0, 5.0),
		"ylim": (0.0, 1.0),
		"xticks": [0.0, 5.0],
		"yticks": [0.0, 1.0],
	}
	elc_throat_inset = {
		"xlim": (-0.5, 2.5),
		"ylim": (0.0, 0.5),
		"xticks": [-0.5, 2.5],
		"yticks": [0.0, 0.5],
	}
	elc_exhaust_inset = {
		"xlim": (-1.0, 5.0),
		"ylim": (0.0, 0.5),
		"xticks": [-1.0, 5.0],
		"yticks": [0.0, 0.5],
	}

	panel_idx = 0
	for i, row in enumerate(ROW_SPECS):
		is_electron_row = row.species.name == "elc"
		throat_inset = elc_throat_inset if is_electron_row else ion_throat_inset
		exhaust_inset = elc_exhaust_inset if is_electron_row else ion_exhaust_inset
		data = loaded[(row.run_key, row.species.name)]
		z_idx = _z_indices(data.z, plot_z_locations)
		vpar_norm, mu_norm, _ = _normalize_axes(data, row.species)
		x_nodal = np.outer(vpar_norm, np.ones_like(mu_norm))
		y_nodal = np.outer(np.ones_like(vpar_norm), mu_norm)

		# Match legacy style: first column linear, remaining columns log-scale.
		selected = [data.f[j, :, :] / data.jb[j] for j in z_idx]
		vmax = float(np.max(selected[0]))
		vmin_log = max(1e-9, float(np.min(selected[0][selected[0] > 0])))

		for j, zi in enumerate(z_idx):
			panel = selected[j]
			z_here = float(data.z[zi])
			vpar_loss, mu_loss = None, None
			if z_here < 0.97:
				vpar_loss, mu_loss = _loss_cone_curve(data, row.species, z_eval=z_here)

			is_linear_panel = j == 0 and i != 2
			norm = Normalize(vmin=0.0, vmax=vmax) if is_linear_panel else LogNorm(vmin=vmin_log, vmax=vmax)
			pcm = ax[i, j].pcolormesh(
				x_nodal,
				y_nodal,
				panel,
				cmap="inferno",
				norm=norm,
				rasterized=True,
			)
			if vpar_loss is not None and mu_loss is not None:
				ax[i, j].plot(vpar_loss, mu_loss, color="white", linestyle="--", linewidth=1.2)

			# Mirror throat / exhaust insets to match legacy distribution plot style.
			if np.isclose(z_here, 0.98, atol=0.03):
				_add_zoom_inset(
					ax[i, j],
					x_nodal,
					y_nodal,
					panel,
					norm,
					**throat_inset,
				)
			if np.isclose(z_here, 2.5, atol=0.08):
				_add_exhaust_inset(
					ax[i, j],
					x_nodal,
					y_nodal,
					panel,
					norm,
					**exhaust_inset,
				)

			ax[i, j].set_xlim(-8.0, 8.0)
			ax[i, j].set_ylim(0.0, 8.0)
			ax[i, j].set_title(rf"$z={plot_z_locations[j]:.2f}\,\mathrm{{m}}$")
			ax[i, j].text(-7.5, 6.7, _panel_label(panel_idx), color="white", fontsize=14)
			panel_idx += 1

			# Keep left-column linear colorbars and rightmost-column colorbars.
			if is_linear_panel or j == n_cols - 1:
				cbar = plt.colorbar(pcm, ax=ax[i, j], orientation="vertical", pad=0.02)
				cbar.ax.tick_params(labelsize=9)
				if j == n_cols - 1:
					cbar.set_label(r"$f$", fontsize=14)

		ax[i, 0].set_ylabel(row.row_label + "\n" + r"$\mu / \mu_t$")

	for j in range(n_cols):
		ax[-1, j].set_xlabel(r"$v_{||}/v_t$")

	fig.tight_layout()
	out_file = out_dir / "distribution-functions-comparison.pdf"
	fig.savefig(out_file)
	plt.close(fig)
	print("[distributions] Saved figure to", out_file)
	print("[distributions] z-locations:", plot_z_locations)


def main(context: dict[str, Any] | None = None) -> None:
	if context is None:
		raise ValueError("This subfile is intended to be called by master.py with context.")
	run(context)

