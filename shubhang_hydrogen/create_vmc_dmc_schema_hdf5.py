#!/usr/bin/env python3
"""Entrypoint for constructing VMC/DMC schema HDF5 files.

This forwards execution to the maintained generator script in:
  .../database_work/configurations/create_vmc_dmc_schema_hdf5.py
"""

from __future__ import annotations

import runpy
from pathlib import Path


def main() -> None:
    here = Path(__file__).resolve()
    target = (
        here.parents[2]
        / "configurations"
        / "create_vmc_dmc_schema_hdf5.py"
    )
    if not target.exists():
        raise FileNotFoundError(f"Could not find generator script at: {target}")
    runpy.run_path(str(target), run_name="__main__")


if __name__ == "__main__":
    main()
