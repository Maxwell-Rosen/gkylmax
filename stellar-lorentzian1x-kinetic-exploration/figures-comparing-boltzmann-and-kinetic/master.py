

"""Master orchestration for Boltzmann vs kinetic postprocessing.

This file owns the shared, high-level configuration and dispatches work to
sub-plot scripts (for example: ``potential.py`` and ``electric-field.py``).

Expected subfile interface:
- Preferred: ``run(context: dict) -> None``
- Fallback: ``main(context: dict) -> None``
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
from types import ModuleType
from typing import Any


@dataclass(frozen=True)
class SimulationInput:
	label: str
	data_path: str
	last_frame: int


@dataclass(frozen=True)
class PlotTask:
	script: str
	enabled: bool = True


@dataclass(frozen=True)
class MasterConfig:
	kinetic: SimulationInput
	boltzmann: SimulationInput
	figures_dir: str
	tasks: tuple[PlotTask, ...]


def build_master_config() -> MasterConfig:
	"""Central place for critical shared settings used by all plot scripts."""
	kinetic = SimulationInput(
		label="kinetic-elc",
		data_path=(
			"/scratch/gpfs/mr1884/scratch/gkylmax/"
			"stellar-lorentzian1x-kinetic-exploration/initial-test"
		),
		last_frame=140,
	)

	boltzmann = SimulationInput(
		label="boltzmann-elc",
		data_path=(
			"/scratch/gpfs/mr1884/scratch/gkylmax/"
			"stellar-lorentzian1x-orbit-average-beams"
		),
		last_frame=65,
	)

	tasks = (
		PlotTask(script="potential.py", enabled=True),
		PlotTask(script="electric-field.py", enabled=True),
		PlotTask(script="electric-field-normalized.py", enabled=True),
		PlotTask(script="assess-theoretical-potential-boltzmann.py", enabled=True),
	)

	return MasterConfig(
		kinetic=kinetic,
		boltzmann=boltzmann,
		figures_dir="figures",
		tasks=tasks,
	)


def validate_paths(config: MasterConfig) -> None:
	"""Fail early if key simulation paths are not present."""
	for sim in (config.kinetic, config.boltzmann):
		path = Path(sim.data_path)
		if not path.exists():
			raise FileNotFoundError(f"Missing input path for {sim.label}: {path}")


def load_plot_module(script_path: Path) -> ModuleType:
	"""Load a plotting script by file path (supports hyphenated file names)."""
	module_name = script_path.stem.replace("-", "_")
	spec = spec_from_file_location(module_name, script_path)
	if spec is None or spec.loader is None:
		raise ImportError(f"Failed to create import spec for: {script_path}")

	module = module_from_spec(spec)
	# Register module before execution so decorators (for example dataclass)
	# can resolve cls.__module__ through sys.modules during import-time checks.
	sys.modules[module_name] = module
	try:
		spec.loader.exec_module(module)
	except Exception:
		sys.modules.pop(module_name, None)
		raise
	return module


def build_context(config: MasterConfig, task: PlotTask) -> dict[str, Any]:
	"""Shared payload passed into each plotting subfile."""
	out_dir = Path(config.figures_dir)
	out_dir.mkdir(parents=True, exist_ok=True)

	return {
		"task": asdict(task),
		"kinetic": asdict(config.kinetic),
		"boltzmann": asdict(config.boltzmann),
		"constants": {
			# Keep key normalization values centralized in the master config context.
			"boltzmann_te0_ev": 940.0,
		},
		"figures_dir": str(out_dir.resolve()),
		"repo_root": str(Path.cwd().resolve()),
	}


def run_task(config: MasterConfig, task: PlotTask) -> None:
	"""Execute one plotting task using the standard subfile interface."""
	if not task.enabled:
		print(f"[skip] {task.script} (disabled)")
		return

	script_path = Path(__file__).with_name(task.script)
	if not script_path.exists():
		raise FileNotFoundError(f"Plot task script not found: {script_path}")

	module = load_plot_module(script_path)
	context = build_context(config, task)

	if hasattr(module, "run"):
		print(f"[run] {task.script} via run(context)")
		module.run(context)  # type: ignore[attr-defined]
		return

	if hasattr(module, "main"):
		print(f"[run] {task.script} via main(context)")
		module.main(context)  # type: ignore[attr-defined]
		return

	raise AttributeError(
		f"{task.script} must define run(context) or main(context)."
	)


def main() -> None:
	config = build_master_config()
	validate_paths(config)

	print("Starting postprocessing tasks...")
	for task in config.tasks:
		run_task(config, task)
	print("All enabled tasks finished.")


if __name__ == "__main__":
	main()