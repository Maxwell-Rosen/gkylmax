#!/usr/bin/env bash

set -euo pipefail

print_help() {
	cat <<'EOF'
Usage:
	./make-initial-condition.sh RESZ RESVPAR RESMU

Description:
	Copies required initial condition files into:
		initial-condition-RESZ-RESVPAR-RESMU

Required files in the current directory:
	- gk_lorentzian_mirror-ion_0.gkyl
	- gk_lorentzian_mirror-ion_jacobvel.gkyl
	- gk_lorentzian_mirror-jacobtot_inv.gkyl

Examples:
	./make-initial-condition.sh 96 48 32
	./make-initial-condition.sh nz96 nv64 nmu24
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
	print_help
	exit 0
fi

if [[ "$#" -lt 3 ]]; then
	echo "Error: not enough resolution arguments." >&2
	echo "Need 3 values: RESZ RESVPAR RESMU" >&2
	echo >&2
	print_help >&2
	exit 1
fi

if [[ "$#" -gt 3 ]]; then
	echo "Error: too many arguments." >&2
	echo "Expected exactly 3 values: RESZ RESVPAR RESMU" >&2
	echo >&2
	print_help >&2
	exit 1
fi

resz="$1"
resvpar="$2"
resmu="$3"

target_dir="initial-condition-${resz}-${resvpar}-${resmu}"

files=(
	"gk_lorentzian_mirror-ion_0.gkyl"
	"gk_lorentzian_mirror-ion_jacobvel.gkyl"
	"gk_lorentzian_mirror-jacobtot_inv.gkyl"
)

for f in "${files[@]}"; do
	if [[ ! -f "$f" ]]; then
		echo "Error: required file not found: $f" >&2
		exit 1
	fi
done

mkdir -p "$target_dir"
cp -v "${files[@]}" "$target_dir/"

echo "Copied files to: $target_dir"
