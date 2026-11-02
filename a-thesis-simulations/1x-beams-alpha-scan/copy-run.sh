#!/usr/bin/env bash
# Prepare independent alpha cases without replacing existing inputs or output.
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
mode=prepare
machine=stellar-ai
dry_run=false
alphas=()

usage() {
  cat <<'EOF'
Usage: bash copy-run.sh [--prepare|--build|--submit] [--dry-run]
                        [--machine stellar-ai|perlmutter] [ALPHA ...]

Defaults: prepare 0.3e-3 1e-3 2.5e-3 5e-3 10e-3, using Stellar AI.
  --prepare  Create missing cases from core/; preserve existing cases (default).
  --build    Prepare and rebuild unfinished cases with the cluster modules.
  --submit   Prepare, rebuild, and submit unfinished cases with sbatch.
  --dry-run  Report actions without creating, building, or submitting anything.

Examples:
  bash copy-run.sh --submit --dry-run
  bash copy-run.sh --submit
  bash copy-run.sh --submit 2.5e-3 5e-3 10e-3
  bash copy-run.sh --prepare 7.5e-3

Existing ion and field checkpoints at frame 75 mark a completed case.
Job scripts resume unfinished cases from their latest ion checkpoint.
EOF
}

while (( $# )); do
  case "$1" in
    --prepare|--build|--submit) mode=${1#--}; shift ;;
    --dry-run) dry_run=true; shift ;;
    --machine)
      [[ $# -ge 2 ]] || { echo 'Missing --machine value' >&2; exit 1; }
      machine=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    -*) echo "Unknown option: $1" >&2; exit 1 ;;
    *) alphas+=("$1"); shift ;;
  esac
done
case "$machine" in
  stellar-ai|perlmutter) ;;
  *) echo "Unknown machine: $machine" >&2; exit 1 ;;
esac
if (( ${#alphas[@]} == 0 )); then
  alphas=(0.3e-3 1e-3 2.5e-3 5e-3 10e-3)
fi

python3 - "$script_dir" "$dry_run" "${alphas[@]}" <<'PY'
from decimal import Decimal
from pathlib import Path
import re
import shutil
import sys
import tempfile

root = Path(sys.argv[1])
dry_run = sys.argv[2] == 'true'
alphas = sys.argv[3:]
files = ('sim.c', 'Makefile', 'jobscript-gkyl-stellar-ai', 'jobscript-gkyl-perlmutter')
pattern = r'(\bdouble alpha_oap\s*=\s*)([^;]+)(;)'
known = {}
for folder in (root / 'alpha-scan').glob('alpha-*'):
    try:
        known[Decimal(folder.name.removeprefix('alpha-'))] = folder.name
    except ArithmeticError:
        pass
# Validate the entire request before making changes.
for alpha in alphas:
    if not re.fullmatch(r'(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', alpha):
        sys.exit(f'Invalid alpha: {alpha}')
    value = Decimal(alpha)
    if not value.is_finite() or value <= 0:
        sys.exit(f'Alpha must be finite and positive: {alpha}')
    name = f'alpha-{alpha}'
    if value in known and known[value] != name:
        sys.exit(f'{alpha} already exists as {known[value]}; use that spelling.')
    known[value] = name
    folder = root / 'alpha-scan' / name
    source = folder if folder.exists() else root / 'core'
    for filename in files:
        if not (source / filename).is_file():
            sys.exit(f'Missing input: {source / filename}')
    matches = re.findall(pattern, (source / 'sim.c').read_text())
    if len(matches) != 1:
        sys.exit(f'Expected one alpha_oap assignment in {source / "sim.c"}')
    if folder.exists() and Decimal(matches[0][1].strip()) != value:
        sys.exit(f'{folder}: alpha_oap does not match its folder; refusing to overwrite it.')

for alpha in alphas:
    folder = root / 'alpha-scan' / f'alpha-{alpha}'
    if folder.exists():
        print(f'KEEP {folder.relative_to(root)} (existing inputs and data)')
        continue
    print(f'CREATE {folder.relative_to(root)} from core/')
    if dry_run:
        continue
    folder.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.prepare-', dir=folder.parent) as temporary:
        staging = Path(temporary) / folder.name
        staging.mkdir()
        for filename in files:
            shutil.copy2(root / 'core' / filename, staging / filename)
        source = staging / 'sim.c'
        source.write_text(re.sub(pattern, lambda m: m[1] + alpha + m[3], source.read_text()))
        for job in staging.glob('jobscript-*'):
            text, count = re.subn(r'^#SBATCH -J .*$', f'#SBATCH -J poa-a-{alpha}',
                                 job.read_text(), flags=re.MULTILINE)
            if count != 1:
                sys.exit(f'Expected one job name in {job}')
            job.write_text(text)
        staging.rename(folder)
PY

[[ $mode != prepare ]] || exit 0
modules_loaded=false
for alpha in "${alphas[@]}"; do
  folder="$script_dir/alpha-scan/alpha-$alpha"
  if [[ -s "$folder/zzim-ion_75.gkyl" && -s "$folder/zzim-field_75.gkyl" ]]; then
    echo "SKIP alpha=$alpha (final ion and field checkpoints already exist)"
    continue
  fi
  echo "BUILD alpha=$alpha ($machine)"
  if [[ $mode == submit ]]; then
    echo "SUBMIT alpha=$alpha with jobscript-gkyl-$machine (resume if checkpoints exist)"
  fi
  [[ $dry_run == false ]] || continue
  if [[ $modules_loaded == false ]]; then
    : "${PS1:=}"
    if [[ $machine == stellar-ai ]]; then
      module purge
      module load nvhpc/26.5
      module load cudatoolkit/13.3
      module load openmpi/nvhpc-26.5/5.0.10
      module load openblas/0.3.x
      export LD_LIBRARY_PATH="${NVHPC_ROOT:?}/comm_libs/nccl/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    else
      module load PrgEnv-gnu/8.6.0
      module load craype-accel-nvidia80
      module load cray-mpich/9.0.1
      module load cudatoolkit/13.0
      module load nccl/2.29.2-cu13
      module load cray-libsci/25.09.0
    fi
    modules_loaded=true
  fi
  (
    cd "$folder"
    # Rebuild against current shared headers/libraries, including migrated runs.
    make -B sim
    if [[ $mode == submit ]]; then
      sbatch "jobscript-gkyl-$machine"
    fi
  )
done
