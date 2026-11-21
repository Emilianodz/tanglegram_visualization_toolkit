#!/usr/bin/env python3
"""Entry point: VAMPhyRE vs GTDB-Tk tanglegram with rogue-branch analysis.

All logic lives in the `tanglegram_utils` package; see `--help` and README.md."""
import sys

from tanglegram_utils.cli import main

if __name__ == "__main__":
    sys.exit(main())
