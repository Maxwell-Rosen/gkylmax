#!/usr/bin/env python3
"""Regenerate the four publication studies using the active Python environment."""
import argparse
from pathlib import Path
import subprocess
import sys

FIGURES = Path(__file__).resolve().parent
STUDIES = ('resolution-scan', 'alpha-scan', 'extra-cycles-fdp-length', 'extra-time-dilation')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=FIGURES.parent,
                        help='Parent directory containing the simulation directories')
    parser.add_argument('--only', nargs='+', choices=STUDIES, default=STUDIES)
    parser.add_argument('--reuse-data', action='store_true')
    parser.add_argument('--dpi', type=int, default=300)
    args = parser.parse_args()
    for study in args.only:
        script = 'publication_convergence.py' if study == 'resolution-scan' else 'publication_figures.py'
        command = [sys.executable, str(FIGURES / study / script), '--dpi', str(args.dpi)]
        if study == 'resolution-scan':
            command += ['--root', str(args.root / '1x-beams-resolution-scan'),
                        '--baseline', str(args.root / '1x-beams')]
        else:
            command += ['--root', str(args.root)]
        if args.reuse_data:
            command += ['--reuse-data']
        subprocess.run(command, check=True)


if __name__ == '__main__':
    main()
