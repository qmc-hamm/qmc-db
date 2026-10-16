# Ilnur hydrogen data — ingest brief for the Leonardo agent

**Read this file first, then `README.md` in this same repository.**
This brief is the operating manual for converting Ilnur’s raw high-pressure
hydrogen calculations on Leonardo into the QMC-HAMM HDF5 database and
uploading the new files to Open Storage Network (OSN).

**Written:** 2026-10-16  
**Owner:** Shubhang Goswami (`sgoswam3`)  
**Repo you are in:** https://github.com/qmc-hamm/qmc-db  
**Companion schema prose:** `README.md` (same directory)  
**Reference writers already in this repo:**

- `shubhang_hydrogen/create_hdf5_schema.py` — Aurora/SG QMCPACK run directories
- `ricky_legacy/ricky_common.py` and `ricky_legacy/build_ricky_database.py` — legacy Girder trajectories
- `ricky_legacy/verify_ricky_transfer.py` — the verification pattern to copy
- `sync_database_to_osn.py` — bucket sync (needs the Illinois index; see §8)

Do not invent a third HDF5 layout. Copy the group names, dataset names, units,
and finite-size-correction convention from those writers.

---

## 0. Mission

You are on Leonardo, logged in as the user who has Ilnur’s raw hydrogen data.
Your job is:

1. Find the raw data tree (the human will point you at a folder, or you search
   under the usual Leonardo project/scratch areas for Ilnur / `isaitov` /
   `HYDROGEN` / `HP_dataset`).
2. Inventory every calculation in that tree before writing any HDF5.
3. Convert **each energy/force evaluation of each configuration** into one
   schema-compliant `{uuid}.h5`.
4. Validate against the known 4,717-frame catalog in §3 (a lower bound, not a cap).
5. Upload **only the new files** to `s3://phy240060/QMCHAMM/`.
6. Leave a ledger, a local index, and a short report of what was ingested,
   what was skipped, and what failed.

The community database already contains Aurora/SG and Ricky-legacy hydrogen.
Ilnur’s frames are **not** in that HDF5 archive yet. Coverage plots on Illinois
only overlay an older XYZ export. You are doing the real ingest.

---

## 1. What the database is

A community archive of **single-configuration** energy and force evaluations
for high-pressure hydrogen. One HDF5 file = one evaluation of one geometry.
Files are uploaded flat to:

```
s3://phy240060/QMCHAMM/<uuid>.h5
endpoint: https://uri.osn.mghpcc.org
```

Design rules:

1. **Self-contained.** Structure, observables, parameters, code inputs, and
   upstream UUIDs live inside the file. After upload, a filesystem path on
   Leonardo is not a substitute for provenance.
2. **Indexable.** Root attributes are the query surface (T, P, rs, method,
   system, uuid, source tags, phase).
3. **Deduplicable.** UUID filenames plus a sha256 ledger. Never upload the
   same content twice. Never mint a new UUID for a frame you already wrote.
4. **Uniform.** Every DMC file exposes the same observable names. Missing
   numbers are NaN. Missing vectors are empty arrays of the right rank.
   Do not omit the dataset.

### Already in the archive (Illinois, August 2026)

| Source | `source_dataset` | Files | What it is |
|---|---|---|---|
| Aurora / Shubhang QMC | `aurora_sg` | 5,463 | VMC + DMC + backing DFT from Nexus/QMCPACK |
| Legacy QMC-HAMM (“Ricky”) | `ricky_legacy` | 54,123 | 17,528 DMC + 36,595 DFT |
| **Total indexed** | | **59,586** | about 3.24 GB |

Those files must not be regenerated or re-uploaded from Leonardo.
Ilnur is a **new** `source_dataset`.

Illinois on-disk copies (you generally cannot see these from Leonardo):

```
/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/
├── database_index.csv
├── database_index_summary.json
├── ricky_legacy_database/P<P>T<T>/*.h5
├── HP_dataset_ilnur_coverage/          # XYZ-overlay plots, NOT the HDF5 DB
└── database_work/qmchamm_database/     # this git repo
```

---

## 2. Scientific decisions already made — do not reopen them

These were confirmed while building the Illinois archive. Apply them.

### 2.1 Ilnur’s `CEIMC` / `DFT` / `M29` prefixes are creation-path tags

Every frame in the known catalog has a prefix ending in `-dmc`.
The token before the phase (`CEIMC`, `DFT`, `M29`) records **how the ionic
configuration was generated**, not the electronic-structure method of the
stored energy.

| Prefix token | What it means | What to write in HDF5 |
|---|---|---|
| `CEIMC-...-dmc` | geometry from a CEIMC path; energy is DMC | `method="DMC"`, `calculation_type="qmc"` |
| `DFT-...-dmc` | geometry from a DFT / AIMD-like path; energy is DMC | same |
| `M29-...-dmc` | geometry from the M29 model path; energy is DMC | same |

Put the creation path in:

- `method_kws` (for example `ceimc_path`, `dft_path`, or `m29_path`, plus any
  ansatz tags you can prove from the inputs, such as `slater_jastrow`)
- `starting_configuration_model_name` (the path tag, or a model name if the
  raw files name one)
- `provenance/config_creation_method` as a UTF-8 string dataset

**Do not** set root `method` to `DFT` just because the folder or prefix says
`DFT`. That mistake was made once in the coverage tables and then corrected:
those frames are DMC results. Coverage-plot colors still use the old labels
`CEIMC` and `DFT` only as plot series. The HDF5 `method` field must be `DMC`.

### 2.2 The only DFT-only hydrogen set already in the database is Ricky BOPIMC

`f15b_bopimc` (2,729 frames) is QE SCF, `method=DFT`, empty `dependent_uuids`,
all near `rs≈1.49`. It is tagged `QMCHAMM_DFT` in plots and is excluded from
default coverage figures. It is unrelated to Ilnur. Do not mix it in.

If Ilnur’s tree also contains genuine DFT SCF or VMC evaluations that are
**not** the DMC energy (separate output files, separate totals), write those
as their own HDF5 files (`method=DFT` or `method=VMC`) and point the DMC
file’s `dependent_uuids` at them. Do not invent a DFT file when you only have
a DMC energy.

### 2.3 One file per evaluation

Never pack a trajectory, a twist grid, or a (P, T) series into one HDF5.
NSCF, `pw2qmcpack`, orbital optimization, and pseudopotential generation are
not database rows. Embed their inputs under `/code` of the energy evaluation
they support.

### 2.4 Finite-size correction

**The stored `observables/total_energy` is the finite-size-corrected value**
when a correction exists.

```
total_energy_uncorrected = total_energy - (fsc_potential_energy + fsc_kinetic_energy)
```

Same convention as Aurora XYZ headers (`energy`, `fsc_dv_ev`, `fsc_dt_ev`)
and the legacy Girder trajectories. Do not invert it.

If you cannot find a correction:

- store the published total as `total_energy`
- set `fsc_potential_energy` and `fsc_kinetic_energy` to NaN
- set `observables.attrs["finite_size_corrected"] = False`
- set `qmc_quality` to `7` for DMC (see §4.6), not `10`

If the raw energy is uncorrected and you do have Δv and Δt in eV, the value
you write to `total_energy` is `raw + Δv + Δt`. Record which file the deltas
came from in provenance. Do not silently assume the XYZ `energy` field is
corrected or uncorrected — check the header keys and any `*energycorrections*`
file and write the choice into provenance.

### 2.5 Pressure on the root attribute is the nominal thermodynamic pressure

Root `pressure` (GPa) is the target / label pressure (the `P250` in the
prefix), not necessarily the virial pressure of this electronic calculation.

`observables/pressure` is the virial pressure from **this** evaluation’s
stress, in GPa. They often differ. If the prefix has no `P` (the DFT-fluid
`RS` series), **omit** root `pressure` or store NaN. Do not back out a fake
nominal pressure from rs.

### 2.6 Phase labels

Known ionic labels in the XYZ catalog:

| Label in the prefix | Frames in the old XYZ export | Natoms |
|---|---|---|
| `PC48` | 1,573 | 96 |
| `FLUID` | 1,145 | 96 |
| `C2C` | 886 | 96 |
| `CMCA12` | 762 | 96 |
| `CSIV` | 351 | 128 |

Store:

- `phase` — use `liquid` when the source label is `FLUID` / `Fluid` / `LIQUID`;
  otherwise keep `C2C`, `CMCA12`, `CSIV`, `PC48`
- `phase_label_raw` — the original token, unchanged

Do not run the Ricky S(k) classifier on Ilnur data unless the human asks.
The Ricky classifier was built for a different legacy set (hcp / liquid /
ambiguous). Ilnur already has crystallographic labels.

### 2.7 `rs`

Wigner–Seitz radius in **bohr**, one electron per proton:

```python
BOHR_PER_ANG = 1.0 / 0.529177210903
r_ang = (3.0 * V_ang3 / (4.0 * math.pi * natoms)) ** (1.0 / 3.0)
rs = r_ang * BOHR_PER_ANG
```

Compute `rs` from the cell you store. If the prefix also carries an `RS`
label (DFT fluid: 1.35, 1.37, 1.39, 1.41, 1.43, 1.45), keep that label in
provenance as `rs_label`. The root `rs` attribute is the value from the cell,
which should agree with the label to ~1e-3 or better. If it does not, stop
and record the mismatch; do not overwrite the cell-derived value with the label.

### 2.8 Author string

Do not invent a legal name. The historical Delta path owner is the account
`isaitov`. Use `author = "Ilnur"` unless a README or input file in the raw
tree gives a fuller name. Record the account and the machine in `/code`.

### 2.9 Source tags (use these exact strings)

| Attribute | Value |
|---|---|
| `source_database` | `leonardo` |
| `source_dataset` | `ilnur_leonardo` |
| `system` | `Hydrogen` |
| `formula` | `H` |
| `species` | `["H"]` |

Also mirror `source_database` and `source_dataset` as datasets under
`/provenance`.

---

## 3. Known catalog (lower bound, from an older XYZ export)

On Delta, two extended-XYZ files were once used **only for coverage plots**.
They were never converted to HDF5:

```
/work/hdd/bcqo/isaitov/HYDROGEN/HP_dataset_01/dataset_smfsc_ceimc_dft.xyz
/work/hdd/bcqo/isaitov/HYDROGEN/HP_dataset_01/dataset_smfsc_m29.xyz
```

Those paths are Delta (`bcqo`), not Leonardo. The same files, or the run
directories that produced them, may have been copied to Leonardo under a
different prefix. Search by filename as well as by directory.

Illinois snapshot of that export (4,717 frames):

```
.../aurora_backup/HP_dataset_ilnur_coverage/frames_TP_rs_ilnur_only.csv
```

You may not have that CSV on Leonardo. The counts below are the acceptance
check. If your inventory is **smaller**, you have not found the full set —
keep looking and say what is missing. If it is **larger**, ingest the extra
runs too. The human asked for all data in Ilnur’s folder, not only this export.

| | Count |
|---|---|
| Frames | 4,717 |
| `CEIMC` path | 2,679 (file `dataset_smfsc_ceimc_dft.xyz`) |
| `DFT` path | 1,715 (same file) |
| `M29` path | 323 (file `dataset_smfsc_m29.xyz`) |
| Atoms | 4,366 frames with 96 atoms; 351 `CSIV` frames with 128 atoms |
| T (K) | 200, 290, 400, 430, 480, 600, 675, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000 |
| Nominal P (GPa), when present | 250, 275, 300, 325, 350, 375, 400, 450, 475, 500, 550, 600 |
| Frames with **no** nominal P | 897, all `DFT-FLUIDT*RS*` |
| Cell-derived rs | about 1.200 to 1.450 bohr |

Prefix grammar (every known frame matched this):

```
^(CEIMC|DFT|M29)-(.+)-config-(\d+)-dmc$
```

The middle token is `<PHASE>T<temperature>P<pressure>` or, for DFT fluid,
`FLUIDT<temperature>RS<rs>`.

Examples:

```
CEIMC-C2CT200P250-config-100-dmc
DFT-FLUIDT600RS1.35-config-1-dmc
M29-FLUIDT2000P450-config-12-dmc
CEIMC-CSIVT200P600-config-51-dmc
```

Extended XYZ frame layout used by the coverage extractor:

```
<natoms>
<key="value" ... Lattice="ax ay az bx by bz cx cy cz" prefix=CEIMC-...-dmc ...>
species x y z
...
```

Header keys the Aurora writer already understands, if they appear, are:
`Lattice`, `pressure`, `temperature`, `config`, `energy`,
`electron_kinetic_energy`, `potential_energy`, `fsc_dv_ev`, `fsc_dt_ev`,
`rs`, `QMC-quality`, `prefix`, `uuid`, `modelname`. Forces, when present,
are three extra columns on each atom line, in eV/Å. Lattice in these XYZ
files is Ångström.

**The XYZ export may contain energies and cells but not the QMCPACK inputs,
pseudopotentials, or per-twist traces.** Prefer the raw run directory when
both exist. Use the XYZ as a checklist and as a fallback for geometry and
scalars when the run directory is missing. If a frame exists only as XYZ,
still write an HDF5, set `code/location` to the XYZ path, and lower
`qmc_quality` if forces or FSC evidence are absent (§4.6).

Configs per thermodynamic condition in the export: median about 46, maximum 100.

---

## 4. HDF5 schema you must match

Authoritative short prose is `README.md`. The production writers are stricter
than the oldest README tables. Follow this section when they disagree with a
casual reading of the README (the README still shows a generic `energy`
dataset name; production files use `total_energy`).

### 4.1 Groups (all required)

```
/
├── .attrs
├── parameters/
├── code/
├── structure/
├── observables/
└── provenance/
```

### 4.2 Root attributes

| Attribute | Type | Rule |
|---|---|---|
| `uuid` | hex string, no dashes | `uuid.uuid4().hex`. Filename is `{uuid}.h5`. |
| `system` | string | `Hydrogen` |
| `formula` | string | `H` |
| `natoms` | int | |
| `species` | bytes/UTF-8 array | `H` |
| `calculation_type` | string | `qmc` for VMC and DMC; `scf` for DFT |
| `method` | string | `DMC`, `VMC`, or `DFT` |
| `method_kws` | string array | functional / ansatz / creation-path tags |
| `temperature` | float, K | from prefix or input |
| `pressure` | float, GPa | nominal label; omit/NaN if the series is rs-only |
| `rs` | float, bohr | from the stored cell |
| `creation_date` | string | date of the **calculation** if known, else `unknown`. Not “today”. |
| `author` | string | `Ilnur` unless the tree says otherwise |
| `name_in_system` | string | the full prefix, e.g. `CEIMC-C2CT200P250-config-100-dmc` |
| `config_number` | int | the integer in `-config-N-` |
| `qmc_quality` | int | VMC/DMC only. **Absent on DFT.** |
| `phase` | string | §2.6 |
| `phase_label_raw` | string | original token |
| `ensemble` | string | only if inputs say `npt` or `nvt` |
| `quantum_ions` | bool | only if inputs prove it. Do not guess from the word CEIMC. |
| `starting_configuration_model_name` | string | creation path or model name |
| `source_database` | string | `leonardo` |
| `source_dataset` | string | `ilnur_leonardo` |

### 4.3 Units

| Quantity | Unit in the file |
|---|---|
| Positions, lattice | Å |
| Energies | eV |
| Forces | eV/Å |
| Stress | GPa, shape (3, 3) |
| Pressure | GPa |
| Temperature | K |
| `rs` | bohr |
| Plane-wave cutoffs in `/parameters` | eV (Ry × 13.605693122994) |

Constants:

```python
BOHR_TO_ANG = 0.529177210903
RY_TO_EV = 13.605693122994
HARTREE_TO_EV = 27.211386245988
EV_PER_A3_TO_GPA = 160.21766208
```

QMCPACK XML lattices are often in bohr. Convert before writing.
ASE stress is often Voigt `(xx, yy, zz, yz, xz, xy)` in eV/Å³ with the
convention stress = −pressure. Virial pressure:

```python
P_GPa = -trace(stress_matrix_GPa) / 3.0
```

If virial pressures come out largely **negative** for dense hydrogen
(`rs < 2`), the sign is probably flipped. Flip stress and
`observables/pressure`, and set
`provenance/stress_sign_convention_flipped = True`. Check a few points
against the nominal P(rs) before flipping a whole dataset.

### 4.4 `/structure`

| Dataset | Shape | Notes |
|---|---|---|
| `lattice_vectors` | (3, 3) | Å, rows are lattice vectors |
| `positions` | (N, 3) | Cartesian Å |
| `fractional_positions` | (N, 3) | |
| `pbc` | (3,) bool | hydrogen cells here are periodic |
| `xyz` | scalar UTF-8 | full extended-XYZ text of **this frame only** |

String datasets: `h5py.string_dtype(encoding="utf-8")`.
String **lists** (`dependent_uuids`, `method_kws`, `species` when stored as
datasets): always a **1-D array**, never a scalar string. This is a real
bug that has broken indexing before.

### 4.5 `/observables`

Errors use the postfix `_error`.

DFT, when you truly have an SCF evaluation:

| Dataset | Shape |
|---|---|
| `total_energy` | (1,) eV |
| `forces` | (1, N, 3) eV/Å |
| `stress` | (3, 3) GPa |
| `pressure` | (1,) virial GPa |

DMC and VMC — write **all** of these; NaN or empty if unknown:

| Dataset | Notes |
|---|---|
| `total_energy`, `total_energy_error` | corrected total when FSC exists |
| `kinetic_energy`, `kinetic_energy_error` | |
| `potential_energy`, `potential_energy_error` | |
| `forces`, `forces_error` | (1, N, 3) |
| `structure_factor`, `structure_factor_error`, `structure_factor_ks` | empty `(0, 3)` if you do not have S(k) |
| `fsc_potential_energy` | Δv in eV (DMC) |
| `fsc_kinetic_energy` | Δt in eV (DMC) |
| `total_energy_uncorrected` | |
| `kinetic_energy_uncorrected`, `potential_energy_uncorrected` | NaN allowed |

Group attribute:

```python
observables.attrs["finite_size_corrected"] = True  # or False
```

Twist-averaged scalars are what belong in these datasets. If the raw file is
a per-twist table (QMCPACK `qmc_scalardata_*.h5` or a scalar.dat per twist),
average over twists the same way the run’s own analysis did, and store the
twist count and weights in `/parameters/other_input` or `/code`. Do not store
one arbitrary twist as if it were the average.

### 4.6 `qmc_quality` (integer, VMC/DMC only)

Match `create_hdf5_schema._compute_qmc_quality`:

| Value | Meaning |
|---|---|
| `1` | No usable total energy, or forces entirely missing |
| `7` | DMC energy exists but there is no finite-size-correction evidence |
| `10` | Energy, forces, and FSC evidence are all present (or the source header explicitly sets `QMC-quality`) |

An explicit `QMC-quality` in a header overrides the heuristic.
DFT files must not contain `qmc_quality`. Validators treat that as failure.

### 4.7 `/provenance`

| Dataset | Role |
|---|---|
| `dependent_uuids` | 1-D UTF-8 array. UUIDs of upstream calculations (DFT behind this DMC). Empty array if none. |
| `uuid_in_system` | original prefix, or an upstream uuid if one exists |
| `source_files` | names or relative paths of the inputs you read |
| `source_database`, `source_dataset` | `leonardo`, `ilnur_leonardo` |
| `config_creation_method` | `CEIMC`, `DFT`, or `M29` |
| `frame_index` | config number |

Link DFT → DMC only after a PBC-aware geometry match (lattice + wrapped
positions, tight tolerance). **Basename matching is not sufficient.** That
collision happened in the Ricky transfer (`f54-ipc` vs `f54-ipq`). Prefer a
stable id (`conf_uuid`, prefix, geometry hash).

### 4.8 `/parameters` and `/code`

Read values from the actual input. Do **not** copy the Aurora defaults
(`kpoint_grid = [6,6,6]`, `time_step = 0.005`, `n_steps = 50`). Those are
specific to Shubhang’s Nexus runs.

- DFT: `kpoint_grid`, `ecutwfc`, `ecutrho`, `sigma`, `smearing`, compact `other_input` JSON
- QMC: `time_step`, `n_steps` or `n_blocks`, `target_walkers`, `kpoint_grid` if known
- `/code` attributes: `name` (likely `QMCPACK`), `version`, `compiler`, `parallelization`
- `/code` datasets: `input_file`, `stdout` (truncate huge logs if needed but record that you truncated and keep the energy table), `system` = `leonardo`, `location` = the raw directory
- Pseudopotential filenames, and the UPF/xml bytes if they are small enough to embed

Keep `/parameters/other_input` small. Large blobs go under `/code` once.

### 4.9 Minimal writer shape

```python
import uuid
import json
import h5py
import numpy as np

uid = uuid.uuid4().hex
dt = h5py.string_dtype("utf-8")
with h5py.File(f"{uid}.h5", "w") as h5:
    h5.attrs["uuid"] = uid
    h5.attrs["system"] = "Hydrogen"
    h5.attrs["formula"] = "H"
    h5.attrs["natoms"] = int(natoms)
    h5.attrs["species"] = np.array(["H"], dtype="S")
    h5.attrs["calculation_type"] = "qmc"
    h5.attrs["method"] = "DMC"
    h5.attrs["method_kws"] = np.array(["ceimc_path", "slater_jastrow"], dtype="S")
    h5.attrs["temperature"] = float(T_K)
    h5.attrs["pressure"] = float(P_GPa)          # skip this line if P is unknown
    h5.attrs["rs"] = float(rs_bohr)
    h5.attrs["creation_date"] = creation_date    # from the run, else "unknown"
    h5.attrs["author"] = "Ilnur"
    h5.attrs["name_in_system"] = prefix          # full prefix
    h5.attrs["config_number"] = int(config_n)
    h5.attrs["phase"] = phase                    # FLUID -> liquid
    h5.attrs["phase_label_raw"] = phase_raw
    h5.attrs["starting_configuration_model_name"] = "CEIMC"
    h5.attrs["source_database"] = "leonardo"
    h5.attrs["source_dataset"] = "ilnur_leonardo"
    h5.attrs["qmc_quality"] = 10                 # or 7 or 1; never on DFT

    # create parameters/, code/, structure/, observables/, provenance/
    gpv = h5.create_group("provenance")
    gpv.create_dataset(
        "dependent_uuids",
        data=np.asarray(dft_uuids, dtype=object),  # empty list is valid
        dtype=dt,
    )
```

---

## 5. How to find and read the raw tree

Do this before any conversion.

1. Ask the human which directory is Ilnur’s if it is ambiguous. Likely clues:
   account or group `isaitov`, directory names `HYDROGEN`, `HP_dataset`,
   `HP_dataset_01`, files `dataset_smfsc_ceimc_dft.xyz`,
   `dataset_smfsc_m29.xyz`, prefixes `CEIMC-`, `DFT-`, `M29-`.
2. `find` that root for: `*.xyz`, `*.in.xml`, `*.struct.xyz`, `qmc_scalardata*.h5`,
   `*energycorrections*`, `scalar.dat`, `*.wfs.xml`, `pw.x` outputs (`*.out`),
   `UPF` / pseudopotentials. Record counts. Do not descend into huge
   wavefunction or restart directories just to hash them.
3. Classify each leaf:
   - **DMC evaluation** of one config (the thing that becomes `method=DMC`)
   - **VMC** evaluation of that same config, if a separate total exists
   - **DFT SCF** that generated orbitals for that config, if outputs exist
   - **geometry-only** or MD trajectory frames with no QMC energy (do not
     upload these as DMC)
   - junk: restarts, band-structure NSCF, failed runs with no energy
4. Write an inventory TSV **before** creating HDF5:

```
inventory.tsv columns:
source_key, kind (dmc/vmc/dft/skip), prefix, T_K, P_GPa, phase_raw,
natoms, rs_or_blank, energy_or_blank, has_forces, has_fsc,
raw_dir, xyz_file, frame_index, notes
```

`source_key` must be stable across reruns, for example:

```
ilnur:<relative-path-of-run-or-xyz>:<frame-index>:<dmc|vmc|dft>
```

5. Compare inventory DMC rows to §3. Report:
   - how many of the 4,717 prefixes you found
   - prefixes in §3 that you did not find
   - extra runs the old XYZ export does not know about

### 5.1 Where the numbers usually live (Aurora layout, for analogy)

Shubhang’s runs, which you should treat as a **pattern**, not as Ilnur’s
directory names:

| Data | Typical file |
|---|---|
| Twist-averaged scalars | `qmc_scalardata_*.h5` |
| Finite-size corrections | `*_energycorrections.txt` (`Δv`, `Δt`) |
| Cell in bohr | `dmc/*.in.xml` |
| Positions | `dmc/*.struct.xyz` (sometimes missing the lattice header — then trust the XML cell) |
| Broken schema symptom | HDF5 with `qmc_quality=1`, NaN energies, `rs` NaN, identity lattice |

If Ilnur’s `*.struct.xyz` has no `Lattice=` and no energy in the header, do
not trust it for the cell. Read the QMCPACK XML. This exact failure produced
1,640 empty Aurora shells (100–210 GPa, 1600–4000 K). A re-ingest script
exists only on the Illinois working tree
(`shubhang_hydrogen/reingest_scalardata_h5.py`); it is **not** necessarily in
git. Do not run it on Ilnur data.

### 5.2 Disk

The existing archive is a few GB because it stores energies, forces, inputs,
and modest pseudopotentials — not orbitals or walker populations. Match that.
If a single config’s raw directory is many GB of `*.h5` wavefunctions, link
to it in `code/location` and do not copy it into the schema file.

Check quota before a full write. Output for ~5–10k files should be on the
order of **a few GB**, not tens of TB. If your estimate is far above that,
you are embedding something you should not.

---

## 6. Conversion procedure

Work in a directory you can afford to rewrite, for example
`$SCRATCH/ilnur_qmchamm/` or a project folder the human names. Suggested
layout (local only; the bucket is flat):

```
ilnur_database/
  P<nominalP>T<T>/<uuid>.h5     # use Pna when nominal P is absent, e.g. PnaT600/
  build_ledger.tsv
  inventory.tsv
  index.csv
  report.md
```

For rs-only fluid frames use `PnaT<T>` so they are not filed under a fake pressure.

### 6.1 Ledger (idempotency)

Append-only TSV:

```
source_key  uuid  relative_h5  sha256  nbytes  status  qmc_quality  utc_timestamp
```

On startup, load the ledger. If `source_key` already has `status=ok` and the
file still exists with the same sha256, skip it. If you must regenerate,
keep the **same uuid** (read it from the ledger) so the bucket object name
does not change. New UUIDs for the same frame cause duplicate uploads.

Use a file lock if you parallelize (`fcntl.flock`), one writer for the ledger.

### 6.2 Per frame

1. Parse prefix → creation method, phase, T, P or rs label, config number.
2. Read cell, positions, species. Convert units to Å.
3. Compute `rs` from the cell. Compare to `rs_label` if any.
4. Read total / kinetic / potential and errors. Apply §2.4.
5. Read forces if present. Shape `(1, N, 3)` eV/Å.
6. Embed a one-frame extended XYZ in `structure/xyz`.
7. Fill `/code` from the input XML and a truncated stdout.
8. Set `qmc_quality`.
9. Write to a temp name, fsync, then rename to `{uuid}.h5`.
10. Append the ledger row.

Parallelism: process-pool over configs is fine. Set
`OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` inside workers.
Do not run one process per twist.

### 6.3 DFT and VMC companions

Only when those evaluations exist as real outputs:

1. Write DFT first (`calculation_type=scf`, no `qmc_quality`).
2. Write VMC with `dependent_uuids = [dft_uuid]` if the orbitals came from that DFT.
3. Write DMC with `dependent_uuids` listing the DFT (and VMC if you want the
   chain). At minimum, DMC must list the DFT uuid when a DFT file was written
   for the same geometry.

Verify geometry before linking (§4.7). A sample of a few hundred links is
enough if the matcher is deterministic; then check that **zero** sampled
links fail.

### 6.4 What “done” means for a frame

A DMC file is acceptable to upload when:

- filename equals `{uuid}.h5` and `attrs["uuid"]` matches
- required groups exist
- `method == "DMC"` and `calculation_type == "qmc"`
- `natoms` matches the position array
- `rs` is finite and in a plausible window for this dataset (about 1.15–1.50
  for the known catalog; extra runs may sit slightly outside — flag, don’t delete)
- `total_energy` is finite and negative for hydrogen at these sizes
  (96 atoms is hundreds of eV in magnitude; if you see ~0 or positive totals,
  you have the wrong field or the wrong unit)
- forces either have a finite array with shape `(1, natoms, 3)` or
  `qmc_quality` is 1
- `dependent_uuids` exists (possibly empty)
- creation-path tag is in `method_kws` and provenance, and is **not** the
  root `method`

---

## 7. Validation before upload

Run this on every new file:

```python
import h5py
from pathlib import Path

REQUIRED_GROUPS = {"parameters", "code", "structure", "observables", "provenance"}
ROOT = {"uuid", "system", "formula", "natoms", "method", "calculation_type",
        "author", "source_dataset", "name_in_system"}

def check(path: Path) -> None:
    with h5py.File(path, "r") as f:
        assert ROOT <= set(f.attrs), path
        assert REQUIRED_GROUPS <= set(f), path
        assert path.name == f"{f.attrs['uuid']}.h5"
        assert str(f.attrs["source_dataset"]) == "ilnur_leonardo"
        method = str(f.attrs["method"]).upper()
        if method in ("VMC", "DMC"):
            assert "qmc_quality" in f.attrs
            assert int(f.attrs["qmc_quality"]) in (1, 7, 10)
            assert "total_energy" in f["observables"]
            assert "finite_size_corrected" in f["observables"].attrs
        if method == "DFT":
            assert "qmc_quality" not in f.attrs
        if method == "DMC":
            # creation path must not have been written as the method
            assert method == "DMC"
        deps = f["provenance"]["dependent_uuids"]
        assert deps.shape[0] >= 0  # 1-D, possibly empty
        assert f["structure"]["positions"].shape[0] == int(f.attrs["natoms"])
```

Then a fidelity check against the source, modeled on
`ricky_legacy/verify_ricky_transfer.py`:

1. **Coverage.** Every inventory row with `kind=dmc` and `status` not `skip`
   has one ledger `ok` row. (P, T, phase) histogram matches §3 for the
   overlapping prefixes.
2. **Schema.** The checker above, plus unique uuids.
3. **Fidelity.** Re-read the raw energy, cell, and forces and compare to the
   HDF5 within a tight tolerance (1e-6 relative on energies, 1e-5 Å on
   positions after wrapping). Do this for every frame if it is cheap, or for
   a stratified sample of at least 200 frames covering every creation path
   and every phase, and say which you did.
4. **Provenance.** Every uuid in `dependent_uuids` exists in the ledger and
   the geometries match under PBC.

Write the results into `report.md`. Do not upload if (1) or (2) fails.
If (3) fails on more than a few frames, stop and fix the reader.

Also confirm you did **not** set `method=DFT` on the 1,715 `DFT-*` path frames
or `method=CEIMC` on anything. `CEIMC` is not a legal `method` value.
Legal methods in this archive are `DFT`, `VMC`, and `DMC`.

---

## 8. Upload to OSN

Bucket:

```
s3://phy240060/QMCHAMM/
endpoint: https://uri.osn.mghpcc.org
```

Credentials (not in git, not in this file):

- RW key for upload, RO key for list/check
- https://coldfront.osn.mghpcc.org/allocation/744/
- parent project https://coldfront.osn.mghpcc.org/project/329/

Scripts in this repo:

- `upload_h5_to_aws.sh` — prompts for a source directory, dry-run defaults on
  (`DO_DRYRUN=1`). Set `DO_DRYRUN=0` only after the dry-run looks right.
- `osn_check_and_restore_h5.sh` — `list`, `check`, `restore-missing`, `restore-file`
- `sync_database_to_osn.py` — the Illinois indexer-aware syncer. It expects
  `database_index.csv` of the **existing** 59k files. On Leonardo, prefer
  `upload_h5_to_aws.sh` for the new Ilnur tree, after you have listed the
  bucket and subtracted uuids that are already there.

Upload rules:

1. Dry-run first.
2. Skip any object whose key `<uuid>.h5` is already in the bucket.
3. Skip any local file whose sha256 is already in your ledger as uploaded.
4. Upload only `ilnur_leonardo` files you just wrote. Do not upload Aurora
   or Ricky files if copies happen to be on Leonardo.
5. Record uploaded keys in `upload_registry.csv` (`uuid, sha256, nbytes, utc`).
6. Credentials stay in the process environment. Do not write them into the
   repo, the report, or the HDF5 files.
7. The bucket layout is flat. Local `P*T*/` folders are for humans.

`upload_h5_to_aws.sh` calls `python3 -m awscli`. If that module is missing,
install it in a user environment or use an equivalent `aws s3 cp` against
`--endpoint-url https://uri.osn.mghpcc.org`.

After a successful upload, the Illinois side still needs:

- a copy or rsync of `ilnur_database/` back under
  `/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/`
- `python3 ricky_legacy/build_database_index.py --hash` so
  `database_index.csv` learns the new files
- a refresh of `HP_dataset_ilnur_coverage/` so plots read HDF5 instead of
  only the XYZ overlay

You can do the upload from Leonardo. Leave the Illinois index rebuild as an
explicit item in `report.md` unless you have the Illinois filesystem mounted.

---

## 9. Coverage plots (context only)

On Illinois, `HP_dataset_ilnur_coverage/update_coverage_with_sg_database.py`
merges indexed HDF5 with the Ilnur XYZ snapshot. Default plots include Ilnur
because he is the only source covering ~250–600 GPa. Without him, T–P plots
stop near 351 GPa.

Plot series names vs HDF5 `method`:

| Plot `method` column | Meaning | HDF5 `method` |
|---|---|---|
| `CEIMC` | Ilnur DMC, CEIMC creation path | `DMC` |
| `DFT` | Ilnur DMC, DFT creation path | `DMC` |
| `M29` | Ilnur DMC, M29 creation path | `DMC` |
| `MLIAP` | Aurora/SG | `DMC` (and related) |
| `QMCHAMM_DMC` | Ricky legacy DMC | `DMC` |
| `QMCHAMM_DFT` | Ricky BOPIMC only | `DFT` |

When the HDF5 ingest exists, future plots should key off
`source_dataset=ilnur_leonardo` plus `provenance/config_creation_method`,
not off a fake `method=DFT`. Until the Illinois index is rebuilt, do not
delete `frames_TP_rs_ilnur_only.csv`.

Latest merged coverage **including** the XYZ overlay (not all of it is in HDF5):
29,898 frames, 200–4000 K, 50–600 GPa. Of those, 4,717 are Ilnur-overlay only.

---

## 10. Other pitfalls already paid for

- **Empty shells.** A schema-shaped HDF5 with NaN observables and an identity
  lattice is a failed ingest, not a successful one. `qmc_quality=1` plus
  finite raw `qmc_scalardata` means “reader looked at the wrong file”.
- **Don’t double-count DFT.** In the Ricky set, ~33.8k DFT files exist only
  as orbital generators for DMC geometries. Coverage counts the DMC frame
  once. Do the same here.
- **Projects quota on Illinois** was near 27/30 TB during the Ricky ingest.
  Keep schema files lean if you copy them back.
- **`qharv` is not available** in the default Illinois Python. Twist averaging
  there was done with pandas/h5py. Don’t depend on `qharv` unless Leonardo
  actually has it.
- **Git state of this repo on Illinois (2026-10-16).** `main` matches
  `origin/main` at `9d25f2c`. Local uncommitted work (Ricky phase-classifier
  edits, `sync_database_to_osn.py`, the scalardata re-ingest script) may not
  be in the commit you pulled. Don’t need those files for Ilnur.
- The sibling repo `configurations/` (https://github.com/qmc-hamm/configurations)
  is an older schema experiment. Production writers are in **this** repo
  (`qmc-db`).

---

## 11. Order of work

1. Locate the raw root. Write down the absolute path in `report.md`.
2. Build `inventory.tsv`. Diff it against §3. Show the human the diff if
   more than a few percent of the 4,717 prefixes are missing, and keep going
   on what you did find.
3. Implement the writer by adapting `shubhang_hydrogen/create_hdf5_schema.py`.
   Put the script in your work directory (or `ilnur_leonardo/` in a branch).
   Do not commit secrets.
4. Convert one frame from each creation path (`CEIMC`, `DFT` fluid, `DFT`
   with a pressure, `M29`, and one `CSIV` 128-atom frame). Open the files
   and check §6.4 by hand.
5. Run the full conversion with a ledger.
6. Run §7.
7. Dry-run the upload. Then upload.
8. Write `report.md` with counts by creation path, phase, `qmc_quality`,
   number uploaded, number skipped, failures, and the Illinois follow-ups
   (index rebuild, coverage refresh, rsync).

If a raw format does not match anything in §3 or §5, stop after the inventory
and describe the format with two redacted examples (paths, header lines,
dataset keys). Do not guess a unit conversion for an unfamiliar energy column.

---

## 12. Do / don’t

**Do**

- One evaluation → one `{uuid}.h5`.
- Set `method=DMC` for every `*-dmc` energy, and record CEIMC/DFT/M29 as the
  creation path.
- Keep DMC observable names uniform, with NaN placeholders.
- Treat a published QMC `energy` as FSC-corrected only after you have checked;
  record the decision.
- Tag files `source_dataset=ilnur_leonardo`.
- Keep a ledger and sha256s before upload.
- Geometry-check provenance links under PBC.
- Ingest runs that are **not** in the 4,717-frame export, if they are real
  QMC evaluations in the same folder.

**Don’t**

- Store a multi-frame XYZ as a single database object.
- Put `qmc_quality` on DFT files.
- Set `method` to `CEIMC`, `M29`, or `DFT` for a DMC energy.
- Trust basename-only DFT↔DMC links.
- Copy Aurora’s hardcoded k-point grid and timestep into Ilnur files.
- Embed orbitals, walker dumps, or restart files.
- Upload the existing Aurora/Ricky objects again.
- Write OSN keys into git or into the HDF5.
- Delete or rewrite the Illinois XYZ coverage CSV from Leonardo.
- Invent `quantum_ions`, ensemble, or stress signs you cannot justify from
  the inputs.
