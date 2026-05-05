"""List root-level attributes of every HDF5 file in s3://phy240060/QMCHAMM/."""

from __future__ import annotations

import os
from pathlib import Path

import h5py
import pandas as pd
import s3fs
from dotenv import load_dotenv
from rich.console import Console
from rich.table import Table

S3_PREFIX = "phy240060/QMCHAMM/"


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


def safe_str(val: object) -> str:
    try:
        if pd.isna(val):
            return ""
    except (TypeError, ValueError):
        pass
    return str(val)


def build_dataframe(fs: s3fs.S3FileSystem, h5_files: list[str]) -> pd.DataFrame:
    records = []
    errors = []
    for path in h5_files:
        name = path.rsplit("/", 1)[-1]
        try:
            attrs = collect_attrs(fs, path)
            records.append({"filename": name, **attrs})
        except Exception as exc:
            errors.append((name, str(exc)))

    df = pd.DataFrame(records).set_index("filename") if records else pd.DataFrame()
    return df, errors


def main() -> None:
    console = Console()
    cfg = load_env()
    fs = make_fs(cfg)

    h5_files = sorted(p for p in fs.ls(S3_PREFIX) if p.endswith(".h5"))
    if not h5_files:
        console.print(f"[yellow]No .h5 files found under s3://{S3_PREFIX}[/yellow]")
        return

    df, errors = build_dataframe(fs, h5_files)

    for name, exc in errors:
        console.print(f"[red]Error reading {name}: {exc}[/red]")

    if df.empty:
        console.print("[yellow]No attributes collected.[/yellow]")
        return

    if "system" not in df.columns:
        console.print("[yellow]No 'system' attribute found in any file.[/yellow]")
        return

    for system, group in df.groupby("system"):
        table = Table(title=f"system = {system}", show_lines=True)
        table.add_column("File", style="cyan", no_wrap=True)
        for col in df.columns:
            if col != "system":
                table.add_column(col, style="white")

        for filename, row in group.iterrows():
            other_vals = [safe_str(row[col]) if col in row.index else "" for col in df.columns if col != "system"]
            table.add_row(filename, *other_vals)

        console.print(table)


if __name__ == "__main__":
    main()