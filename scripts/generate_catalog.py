"""Generate database.csv from HDF5 root attributes in s3://phy240060/QMCHAMM/."""

from __future__ import annotations

import os
from pathlib import Path

import h5py
import pandas as pd
import s3fs
from dotenv import load_dotenv

S3_PREFIX = "phy240060/QMCHAMM/"
OUTPUT_FILE = Path(__file__).resolve().parent.parent / "catalog.csv"


def load_env() -> dict[str, str]:
    repo_root = Path(__file__).resolve().parent.parent
    load_dotenv(repo_root / ".env")
    return {
        "endpoint_url": os.environ["S3_ENDPOINT_URL"],
        "key": os.environ["AWS_ACCESS_KEY_ID"],
        "secret": os.environ["AWS_SECRET_ACCESS_KEY"],
    }


def make_fs(cfg: dict[str, str]) -> s3fs.S3FileSystem:
    return s3fs.S3FileSystem(
        key=cfg["key"],
        secret=cfg["secret"],
        client_kwargs={"endpoint_url": cfg["endpoint_url"]},
    )


def collect_attrs(fs: s3fs.S3FileSystem, path: str) -> dict[str, object]:
    with fs.open(path, "rb") as f:
        with h5py.File(f, "r") as h5:
            return {k: h5.attrs[k] for k in h5.attrs.keys()}


def main() -> None:
    cfg = load_env()
    fs = make_fs(cfg)

    h5_files = sorted(p for p in fs.ls(S3_PREFIX) if p.endswith(".h5"))
    if not h5_files:
        print(f"No .h5 files found under s3://{S3_PREFIX}")
        return

    records = []
    for path in h5_files:
        uri = f"s3://{path}"
        try:
            attrs = collect_attrs(fs, path)
            records.append({"uri": uri, **attrs})
        except Exception as exc:
            print(f"Error reading {uri}: {exc}")

    if not records:
        print("No attributes collected.")
        return

    df = pd.DataFrame(records)
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"Wrote {len(df)} rows to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()