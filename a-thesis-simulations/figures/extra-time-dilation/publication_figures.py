#!/usr/bin/env python3
"""Generate the publication figures for this study; see ../README.md."""
import sys
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent))
from publication_scans import main

if __name__ == '__main__':
    main('extra-time-dilation', SCRIPT_DIR)
