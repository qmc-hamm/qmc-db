# Legacy QMC-HAMM database (Ricky) ingestion

Converts the published QMC-HAMM dataset from <https://qmc-hamm.hub.yt/data.html>
(hosted on `girder.hub.yt`) into the project HDF5 schema, so it sits alongside
the Aurora/SG hydrogen data in one indexable, uploadable database.

## What was ingested

| | source frames | HDF5 files written |
|---|---|---|
| DMC (`b_qmc`) | 17,528 | 17,528 |
| DFT (`a_dft`) | 36,595 | 36,595 |
| **total** | | **54,123** (2.4 GB) |

Covering 65 (P,T) conditions: T = 600–2200 K, P = 50–210 GPa, rs = 1.40–1.84.

The DFT set is the subset of `a_dft` frames that actually back a DMC
configuration, plus the standalone BOPIMC hcp trajectory
(`f15b_bopimc-hcp-n96-rs1.49-t800.traj`). DFT frames with no QMC child were not
converted, which is what keeps this at gigabytes rather than terabytes.

## How to tell this data apart

Every file carries these root attributes:

```
source_database = 'qmc-hamm.hub.yt'
source_dataset  = 'ricky_legacy'
author          = 'Hongwei Niu (Ricky); QMC-HAMM Team'
```

with the full origin (Girder item id, item name, config index, trajectory file)
recorded under `/provenance`.

## Layout

```
aurora_backup/ricky_legacy_database/P<P>T<T>/<uuid>.h5
aurora_backup/ricky_legacy_database/build_database_ledger.txt
```

## Conventions worth knowing

* **Finite-size corrections.** The published `energy` is already corrected.
  `observables/total_energy` holds it as-is; the uncorrected value is written to
  `total_energy_uncorrected = energy - (dv + dt)`, with `dv`/`dt` kept as
  `fsc_potential_energy` / `fsc_kinetic_energy`. `qmc_quality = 10`.
* **DMC observables.** The legacy `dmc_mean` trajectories publish only the total
  energy. Kinetic/potential/structure-factor fields exist but are `NaN`/empty so
  the whole database exposes one uniform set of names.
* **Provenance links.** Each DMC file's `dependent_uuids` points at the DFT
  file(s) for the same geometry. Links are resolved through the metadata
  `conf_uuid` (base-name matching collides between `f54-ipc` and `f54-ipq*`) and
  then confirmed by comparing lattice and positions under PBC.
* **BOPIMC stress sign.** That trajectory stores stress with the opposite sign
  convention to ASE, which yielded pressures near −156 GPa. Files where the
  virial pressure came out negative at rs < 2.0 have stress and pressure negated
  and carry `provenance/stress_sign_convention_flipped = True`.
* **Classical MD** (`f54-hs1`) is included, distinguished by the `ensemble`,
  `quantum_ions` and `config_dft` root attributes.
* **Phase.** All legacy DMC configs and the BOPIMC DFT set are labeled
  `phase="hcp"`. Ionic S(k) on sampled DMC frames shows solid-like Bragg peaks
  (S_max often 20–40, comparable to BOPIMC), especially at T ≲ 1200 K.

## Scripts

| script | purpose |
|---|---|
| `ricky_common.py` | shared constants, HDF5 writers, extxyz rendering, unit/stress conversion, item-name parsing |
| `build_ricky_database.py` | the converter: plans DFT↔QMC links, then writes both phases in parallel |
| `verify_ricky_transfer.py` | four checks — coverage vs source, schema, fidelity against the original trajectories, provenance link integrity |
| `fix_bopimc_stress_sign.py` | standalone patch for the BOPIMC sign issue (now also handled inline by the converter) |
| `build_database_index.py` | indexes **all** sources into `aurora_backup/database_index.csv` |

## Rebuilding / re-verifying

```bash
PY=/scratch/sgoswam3/ricky_qmchamm/venv/bin/python

# convert (idempotent: the ledger is consulted, existing files are skipped)
$PY build_ricky_database.py --workers 32

# check it
$PY verify_ricky_transfer.py

# refresh the whole-database index, with content hashes for the upload registry
$PY build_database_index.py --hash --workers 32
```

Then refresh the coverage plots:

```bash
cd ../../../HP_dataset_ilnur_coverage
$PY update_coverage_with_sg_database.py
```

## Upload

See `../sync_database_to_osn.py` and `../submit_database_sync.sh`. The registry
at `aurora_backup/upload_registry.csv` plus a live bucket listing ensure nothing
is sent twice.
