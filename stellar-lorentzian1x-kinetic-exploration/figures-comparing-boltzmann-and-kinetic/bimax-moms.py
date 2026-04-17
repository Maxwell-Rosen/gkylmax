"""Plot BiMaxwellian moments comparing Boltzmann-ion and kinetic species.

Overlayed traces:
- Kinetic ions      : solid red
- Kinetic electrons : solid blue
- Boltzmann ions    : dashed red (drawn last so it stays visible)
"""

from __future__ import annotations

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


E_CHARGE = 1.602176634e-19
MASS_ION = 2.014 * 1.67262192369e-27
MASS_ELC = 9.1093837015e-31
TE0_EV = 940.0
CS_ION = np.sqrt((TE0_EV * E_CHARGE) / MASS_ION)

COLOR_ION = "tab:red"
COLOR_ION_BOLTZ = "orangered"
COLOR_ELC = "tab:blue"


def _to_1d(arr: Any) -> np.ndarray:
	out = np.asarray(arr)
	out = np.squeeze(out)
	if out.ndim == 0:
		out = out.reshape(1)
	return out


def _find_prefix(data_path: Path, frame: int) -> str:
	candidates = sorted(data_path.glob(f"*-field_{frame}.gkyl"))
	if not candidates:
		raise FileNotFoundError(
			f"No field file found for frame {frame} in {data_path}"
		)
	return candidates[0].stem.split("-field_", 1)[0]


def _find_map_file(data_path: Path, prefix: str) -> Path:
	candidates = [
		data_path / f"{prefix}-mapc2p.gkyl",
		data_path / f"{prefix}-mc2nu_pos.gkyl",
		data_path / f"{prefix}-mc2nu_pos_deflated.gkyl",
	]
	for c in candidates:
		if c.exists():
			return c
	raise FileNotFoundError(f"No map file found for {prefix} in {data_path}")


def _load_z_nodes(map_file: Path) -> np.ndarray:
	data = pg.GData(str(map_file))
	interp = pg.GInterpModal(data, 1, "ms")
	_, nodes = interp.interpolate(1)
	return _to_1d(nodes)


def _load_moments(data_path: Path, prefix: str, species: str, frame: int) -> dict[str, np.ndarray]:
	file_path = data_path / f"{prefix}-{species}_BiMaxwellianMoments_{frame}.gkyl"
	if not file_path.exists():
		raise FileNotFoundError(f"Missing moments file: {file_path}")

	data = pg.GData(str(file_path))
	interp = pg.GInterpModal(data, 1, "ms")

	_, dens = interp.interpolate(0)
	_, upar = interp.interpolate(1)
	_, tpar_div_m = interp.interpolate(2)
	_, tperp_div_m = interp.interpolate(3)

	mass = MASS_ION if species == "ion" else MASS_ELC
	return {
		"dens": _to_1d(dens),
		# Match legacy plotting convention: normalize all u_parallel by ion sound speed.
		"upar": _to_1d(upar) / CS_ION,
		"tpar": _to_1d(tpar_div_m) * mass / E_CHARGE / 1e3,
		"tperp": _to_1d(tperp_div_m) * mass / E_CHARGE / 1e3,
	}


def _split_midplane(z: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
	left = z < 0.0
	return z[left], y[left], z[~left], y[~left]


def _plot_row(ax, z: np.ndarray, moms: dict[str, np.ndarray], *, color: str, ls: str, lw: float, label: str | None = None) -> None:
	style = dict(color=color, linestyle=ls, linewidth=lw)
	y_keys = ("dens", "upar", "tpar", "tperp")
	panel_coords = ((0, 0), (0, 2), (1, 0), (1, 2))

	for idx, (key, (r, c)) in enumerate(zip(y_keys, panel_coords)):
		zl, yl, zr, yr = _split_midplane(z, moms[key])
		if key == "upar":
			# Right panel is log-scale; use magnitude to keep values plottable.
			yr = np.abs(yr)
		if key == "tpar":
			yr = np.abs(yr)
		kw = style.copy()
		if idx == 0 and label is not None:
			kw["label"] = label
		ax[r, c].plot(zl, yl, **kw)
		ax[r, c + 1].plot(zr, yr, **kw)


def _configure_axes(ax, x_min: float, x_max: float) -> None:
	# density
	ax[0, 0].set_ylabel(r"$n$ (m$^{-3}$)")
	ax[0, 0].set_xlim(x_min, 0.0)
	ax[0, 0].set_ylim(0.0, 3.2e19)
	ax[0, 1].set_yscale("log")
	ax[0, 1].yaxis.tick_right()
	ax[0, 1].yaxis.set_label_position("right")
	ax[0, 1].set_xlim(0.0, x_max)
	ax[0, 1].set_ylim(1e13, 4e19)

	# parallel velocity
	ax[0, 2].set_ylabel(r"$u_{||} / c_s$")
	ax[0, 2].set_xlim(x_min, 0.0)
	ax[0, 2].set_ylim(-7.5, 7.5)
	ax[0, 3].set_yscale("log")
	ax[0, 3].yaxis.tick_right()
	ax[0, 3].yaxis.set_label_position("right")
	ax[0, 3].set_xlim(0.0, x_max)
	ax[0, 3].set_ylim(1e-3, 10.0)

	# T_parallel
	ax[1, 0].set_ylabel(r"$T_{||}$ (keV)")
	ax[1, 0].set_xlim(x_min, 0.0)
	ax[1, 0].set_ylim(0.0, 15.0)
	ax[1, 1].set_yscale("log")
	ax[1, 1].yaxis.tick_right()
	ax[1, 1].yaxis.set_label_position("right")
	ax[1, 1].set_xlim(0.0, x_max)
	ax[1, 1].set_ylim(1e-3, 20.0)

	# T_perp
	ax[1, 2].set_ylabel(r"$T_{\perp}$ (keV)")
	ax[1, 2].set_xlim(x_min, 0.0)
	ax[1, 2].set_ylim(0.0, 30.0)
	ax[1, 3].set_yscale("log")
	ax[1, 3].yaxis.tick_right()
	ax[1, 3].yaxis.set_label_position("right")
	ax[1, 3].set_xlim(0.0, x_max)
	ax[1, 3].set_ylim(1e-3, 100.0)

	ax[1, 1].set_xlabel(r"$z$ (m)")
	ax[1, 1].xaxis.set_label_coords(0.0, -0.15)
	ax[1, 3].set_xlabel(r"$z$ (m)")
	ax[1, 3].xaxis.set_label_coords(0.0, -0.15)


def _add_mirror_throat_lines(ax) -> None:
	kw = dict(color="grey", linestyle="--", linewidth=1)
	ax[0, 0].plot([-0.98, -0.98], [0, 4e19], **kw)
	ax[0, 1].plot([0.98, 0.98], [1e13, 1e20], **kw)
	ax[0, 2].plot([-0.98, -0.98], [-10, 10], **kw)
	ax[0, 3].plot([0.98, 0.98], [1e-3, 10], **kw)
	ax[1, 0].plot([-0.98, -0.98], [0, 15000], **kw)
	ax[1, 1].plot([0.98, 0.98], [1e-3, 15000], **kw)
	ax[1, 2].plot([-0.98, -0.98], [0, 35000], **kw)
	ax[1, 3].plot([0.98, 0.98], [1e-3, 3.5e6], **kw)


def _merge_subplot_pairs(ax) -> None:
	"""Replicate moments_poa.py pair-merging geometry for cols (0,1) and (2,3)."""
	for i in range(2):
		for j in range(4):
			pos = ax[i, j].get_position()
			ax[i, j].set_position([pos.x0, pos.y0, pos.width * 1.5, pos.height])

	for i in range(2):
		for j in range(4):
			pos = ax[i, j].get_position()
			if j == 0:
				ax[i, j].set_position([pos.x0 + 0.02, pos.y0, pos.width, pos.height])
			elif j == 1:
				left_edge = ax[i, 0].get_position().x0 + ax[i, 0].get_position().width
				ax[i, j].set_position([left_edge, pos.y0, pos.width, pos.height])
			elif j == 3:
				left_edge = ax[i, 2].get_position().x0 + ax[i, 2].get_position().width
				ax[i, j].set_position([left_edge, pos.y0, pos.width, pos.height])


def run(context: dict[str, Any]) -> None:
	out_dir = Path(context["figures_dir"]) / "bimax-moms"
	out_dir.mkdir(parents=True, exist_ok=True)

	kin_cfg = context["kinetic"]
	boltz_cfg = context["boltzmann"]

	kin_path = Path(kin_cfg["data_path"])
	boltz_path = Path(boltz_cfg["data_path"])
	kin_frame = int(kin_cfg["last_frame"])
	boltz_frame = int(boltz_cfg["last_frame"])

	kin_prefix = _find_prefix(kin_path, kin_frame)
	boltz_prefix = _find_prefix(boltz_path, boltz_frame)

	z_kin = _load_z_nodes(_find_map_file(kin_path, kin_prefix))
	z_boltz = _load_z_nodes(_find_map_file(boltz_path, boltz_prefix))

	kin_ion = _load_moments(kin_path, kin_prefix, "ion", kin_frame)
	kin_elc = _load_moments(kin_path, kin_prefix, "elc", kin_frame)
	boltz_ion = _load_moments(boltz_path, boltz_prefix, "ion", boltz_frame)

	fig, ax = plt.subplots(2, 4, figsize=(10, 6))

	# Draw solids first.
	_plot_row(
		ax,
		z_kin,
		kin_ion,
		color=COLOR_ION,
		ls="-",
		lw=1.5,
		label="Ion (Kin $e^-$)",
	)
	_plot_row(
		ax,
		z_kin,
		kin_elc,
		color=COLOR_ELC,
		ls="-",
		lw=1.5,
		label="Elc (Kin $e^-$)",
	)
	# Draw dashed Boltzmann-ion last so it stays on top.
	_plot_row(
		ax,
		z_boltz,
		boltz_ion,
		color=COLOR_ION_BOLTZ,
		ls="--",
		lw=1.6,
		label="Ion (Boltzmann $e^-$)",
	)

	x_min = float(min(np.min(z_kin), np.min(z_boltz)))
	x_max = float(max(np.max(z_kin), np.max(z_boltz)))
	_configure_axes(ax, x_min, x_max)
	_add_mirror_throat_lines(ax)

	handles, labels = ax[0, 0].get_legend_handles_labels()
	ax[0, 3].legend(handles, labels, loc="lower right", framealpha=1.0)

	fig.tight_layout()
	_merge_subplot_pairs(ax)
	out_file = out_dir / "bimax-moms-kinetic-vs-boltzmann.pdf"
	fig.savefig(out_file)
	plt.close(fig)
	print("[bimax-moms] Saved figure to", out_file)


def main(context: dict[str, Any] | None = None) -> None:
	if context is None:
		raise ValueError("This subfile is intended to be called by master.py with context.")
	run(context)

