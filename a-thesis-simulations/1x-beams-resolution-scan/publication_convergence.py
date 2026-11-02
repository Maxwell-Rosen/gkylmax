#!/usr/bin/env python3
"""Compatibility entry point; canonical analysis now lives under figures/."""
from pathlib import Path
import runpy

TARGET = (Path(__file__).resolve().parents[1] / 'figures' / 'resolution-scan'
          / 'publication_convergence.py')

if __name__ == '__main__':
    runpy.run_path(str(TARGET), run_name='__main__')
else:
    # Preserve imports of the original script's public analysis functions.
    globals().update({key: value for key, value in runpy.run_path(str(TARGET)).items()
                      if not key.startswith('__')})
