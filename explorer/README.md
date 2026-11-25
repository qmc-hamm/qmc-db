# QMC-HAMM database explorer

A local, read-only web explorer for the QMC-HAMM HDF5 database. It mirrors the
OSN bucket (`https://uri.osn.mghpcc.org/phy240060`, prefix `QMCHAMM/`) into a
local cache with unsigned HTTPS GETs (no credentials), builds an index, and
serves a browser UI plus an extxyz export.

## Quickstart

```
cd qmc-db
python explorer/harvest.py  # first run about 32 min, 7.3 GB into explorer/qmc-db-cache/h5 (git-ignored)
python explorer/serve.py    # prints the URL
```

Open http://127.0.0.1:8771 (`--port` to change).

`harvest.py` flags: `--cache DIR` (or `QMCDB_CACHE`), `--data DIR` (default
`explorer/data/`, git-ignored), `--workers N` (default 16), `--limit N` (requires `--data <scratch dir>`; refused with the default `explorer/data/` so a partial run cannot overwrite the real index),
`--local-dir DIR` (repeatable; default `../Li7P3S11_QMC_DFT_dataset`),
`--no-download`. `serve.py` flags: `--port`, `--data`, `--cache`.

Outputs in `explorer/data/`: `index.json.gz`, `manifest.json`,
`harvest_report.json`, `eos_fit.json`, `eos_residuals.json`.

## Refresh

Re-run `harvest.py`. It lists the bucket, compares ETags with `manifest.json`,
and downloads only new or changed objects; unchanged files are not re-parsed unless `NORMALIZE_VERSION` (in `normalize.py`, stored in the index as `normalize_version`) differs from the stored one; then every cached file is re-parsed (no download).
With `--no-download` it parses what is already cached and counts changed bucket
objects it did not fetch as `stale` (their previous rows are kept; see
`harvest_report.json`). Restart `serve.py` to load the new index.

## Flags

| Flag | Meaning |
|------|---------|
| `P_estimated` | Pressure was estimated from the EOS fit, not recorded in the file. Shown with hollow markers and `(est.)` in the table. `P_used` is the pressure shown/exported. |
| `P_extrapolated` | The estimate lies outside the density-of-training-data mask of the EOS fit (the (rs, T) bin has too few recorded-pressure training points, or is outside the training range). Shown as `(extrap.)` in the table. |
| `rs_suspect` | Hydrogen rows only: computed `rs` is missing or outside 0.8-3.0, or differs from the recorded `rs` by more than 1%. `None` for non-hydrogen tabs. |
| `fsc_applied` | QMC files only: an `observables/fsc_kinetic_energy` is present or `qmc_quality` is 10. |
| `config_id` | Shared by a DMC/VMC file and the DFT file(s) of the same configuration: `min(uuid)` of the connected component of the geometry-verified provenance links. |
| `linked_dft_uuids` | UUIDs of the DFT files a QMC file depends on. `linked_qmc_uuids` is the reverse (any method other than DFT, including RQMC, counts as QMC). |
| `geom_mismatch` | A link exists but the positions of the two files differ. Mismatched links are excluded from `linked_dft_uuids` / `linked_qmc_uuids` and listed in `geom_mismatch`. |

Files with no T or no P stay in the index and are counted in the no-T/P bar.

## Export (extxyz)

`POST /api/export.xyz` with form field `payload` = JSON `{"uuids": [...], "include_linked_dft": bool, "include_suspect": bool}`.
One frame per file. `rs_suspect` rows (including linked DFT rows added by expansion) are dropped unless `include_suspect` is true (default false). The `X-Export-Summary` response header carries
`planned_frames`, `unknown_uuids`, `missing_links`, `suspect_skipped`. The page shows the planned counts in a confirm dialog before submitting. Other routes: `/`, `/static/*`,
`/data/<file>`, `/api/file/<uuid>`, `/h5/<uuid>.h5`.

Energy rule (verbatim from the constraints):

> extxyz `energy`, any `total_energy` key, and ASE `get_potential_energy()` are **only** the FSC-corrected `observables/total_energy`. The uncorrected value appears only as the info key `energy_uncorrected`. It is never substituted. A missing corrected value gives `energy=NaN`.

Frame contents:

- Atoms: symbols (from `structure/xyz`, else the `species` attribute), `structure/positions`, `structure/lattice_vectors`, `structure/pbc`.
- Calculator: `energy` (eV-schema value from the index, NaN if missing), `forces` when `observables/forces` has shape (N,3), and `stress` (`ricky_legacy` only, see Known data issues) when `observables/stress` is (3,3), (1,3,3) or 6 values. The schema stores stress in GPa; the export converts to ASE convention (eV/A^3, ASE sign, so P = -trace/3). Checked on Ricky DFT files: `-trace(ASE stress)/3 * 160.21766208` equals the file's `observables/pressure` (136.92 and 143.50 GPa for the two fixtures, recorded P 150).
- `arrays`: `forces_error` when `observables/forces_error` is (N,3).
- `info`: `uuid`, `source`, `method`, `T`, `P` (= `P_used`), `P_estimated`, `P_extrapolated`, `rs` (= computed rs), `config_id`, `qmc_quality`, `phase`, `creation_path`, `fsc_applied`, `energy_uncorrected`, `energy_error`, `linked_dft_uuids` and `linked_qmc_uuids` (comma-joined; the latter on DFT frames), `geom_mismatch` (comma-joined, when non-empty), `rs_suspect` (True, when set), `forces_shape_mismatch` (only when forces were omitted for shape). Keys with no value are absent.

Pairing DMC and DFT frames for delta-learning (export with `include_linked_dft`):

```python
import ase.io
frames = ase.io.read("qmchamm_H_N.xyz", ":")
by_cfg = {}
for a in frames:
    by_cfg.setdefault(a.info["config_id"], {}).setdefault(a.info["method"], []).append(a)
pairs = [(g["DMC"][0], g["DFT"][0]) for g in by_cfg.values() if "DMC" in g and "DFT" in g]
dE = [q.get_potential_energy() - d.get_potential_energy() for q, d in pairs]
```

Check the `geom_mismatch` info key and drop NaN energies before training.

## Current data (harvest 2026-11-19, 17:16-17:48 ET)

88,570 rows = 88,370 bucket files + 200 local Li7P3S11 files. Parse errors 0, download errors 0.

- By tab: hydrogen 88,278; graphene 84; defects 8; li7p3s11 200.
- By source: ricky_legacy 54,123; ilnur_leonardo 14,149; kevin_ly_atomic_h 11,982; aurora_sg 8,024; unknown 292 (graphene, defects and Li7P3S11 carry no source attribute).
- By method: DFT 52,360; DMC 25,957; VMC 8,249; RQMC 2,000 (kevin_ly); CASCI 4.
- Links: 55,581 found, 55,423 unioned, 0 missing, 158 `geom_mismatch`; 38,801 `config_id` components.
- `rs_suspect`: 3,394 hydrogen rows (aurora_sg DMC 1,640, VMC 1,640, DFT 86; ilnur_leonardo DFT 28).
- `P_estimated`: aurora_sg 2,556 (all extrapolated); ilnur_leonardo 2,691 (369 extrapolated).
- No T and no P: all 11,982 kevin_ly_atomic_h rows.

## Known data issues

These are properties of the source data, not of the explorer. Nothing is corrected or dropped.

- Li7P3S11: in every `*_QMC.h5` (100 files), `structure/positions` equals the DFT fractional coordinates times 0.529177. This is a units bug in the source data; those files are all `geom_mismatch`.
- 58 Ilnur CEIMC DMC/VMC files are `geom_mismatch` against their linked DFT file.
- Broken Aurora shells: `rs_suspect` flags 1,640 DMC, 1,640 VMC and 86 DFT aurora_sg rows, where computed and recorded rs disagree. Ilnur DFT has 28 more (3,394 in total). The export skips them by default.
- ALL 1,838 Aurora DMC rows with an energy have `total_energy` exactly equal to their VMC parent's. In the raw files the VMC kinetic and total energies equal the DMC FSC-corrected values. This looks like a source bug in the shubhang writer, to be reported upstream.
- 1,690 Aurora DMC rows have E null (1,640 suspect plus 50 others), and their VMC rows are null too. Their export has `energy=NaN`.
- 23 Ilnur DFT files have a zero cell and E=0 (all are `rs_suspect`).
- Stress units are not uniform across sources: Ricky stress matches the schema (GPa, ASE sign). Ilnur `observables/stress` is labelled GPa but its values are about 1e5 times smaller than the recorded P implies, and Aurora stress has a different sign and a zero third row (QE kbar parser). Export rule: stress is exported only for `ricky_legacy` rows (the only source whose stress matches the recorded pressure). For every other source with stress it is omitted, the frame gets info `stress_omitted = "unverified units"`, and `stress_omitted` is counted in the server-side export summary.
- kevin_ly_atomic_h (11,982 files) has no T and no P, so it is absent from T-P plots and counted in the no-T/P bar.
- Defect files carry no per-atom symbols (`structure/xyz` absent, `species` not a single element), so the export skips them and counts them in `symbols_unknown`.

## EOS status

P is estimated only for files with no recorded pressure, using a fit of P(rs, T)
to the recorded-pressure files (`eos_fit.json`): a=4, b=2, cross-validated RMS 0.04
(relative). Training range: P 50-600 GPa, rs 1.2-1.84, T 200-4000 K.

Median |relative residual| by source and creation path (`eos_residuals.json`):

| Source | Path | n | Median |resid| |
|--------|------|---|------|
| aurora_sg | M18 | 1,036 | 4.9% |
| ilnur_leonardo | CEIMC | 2,680 | 2.1% |
| ilnur_leonardo | DFT | 818 | 1.2% |
| ilnur_leonardo | M29 | 323 | 2.5% |
| ricky_legacy | n2p2 blyp | 5,415 | 2.7% |
| ricky_legacy | npt | 10,923 | 2.7% |
| ricky_legacy | nvt | 3,313 | 4.6% |
| ricky_legacy | pbe | 146 | 15% |
| ricky_legacy | vdw-df | 460 | 3.7% |

Backup-EOS trigger (median residual above 5% for a source or region): met.

- ricky pbe is at 15% (n=146).
- Pathology: at rs = 1.29 the fit gives P = 381, 367, 348 GPa at T = 1000, 2000, 3000 K, so dP/dT < 0, which is unphysical. Near rs of about 1.29 the training data only reach T <= 2000 K. All aurora_sg estimated pressures are extrapolated and fall in this region.

The planned backup, a tabulated ab initio EOS (Militzer & Hubbard 2013 or
Chabrier-Debras 2021), has not been built. Treat `P_estimated` values with
`P_extrapolated` set, and those near rs 1.29 at T > 2000 K, as unreliable.

## Tests

The `.h5` fixtures in `explorer/tests/fixtures/` are not in git (git-ignored); get them from the maintainer.

```
cd qmc-db
~/.venv/py312/bin/python -m pytest explorer/tests -q && node --test explorer/tests/fmt.test.mjs
```
