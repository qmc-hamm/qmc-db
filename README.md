

```
/
├── .attrs
├── parameters/
├── code/
├── structure/
├── observables/
└── provenance/        
```

## Root attributes (.attrs)

This should be indexable high-level metadata for database integration. 

| Attribute Name     | Type             | Description                                       |
| ------------------ | ---------------- | ------------------------------------------------- |
| `system`           | string           | Human-readable system name (e.g. "Si_bulk_2x2x2") |
| `formula`          | string           | Chemical formula (e.g. "Si", "Fe2O3")             |
| `natoms`           | int              | Number of atoms                                   |
| `species`          | string[]         | Unique atomic species                             |
| `calculation_type` | string           | "scf", "relax", "md", etc.                        |
| `method`           | string           | Method: "VMC", "DMC", "DFT"              |
| `method_kws`           | string[]     | Exchange correlation functional, ansatz              |
| `temperature`      | float            | Temperature (K), if relevant                      |
| `pressure`         | float            | Pressure (GPa), if relevant                       |
| `creation_date`    | string (ISO8601) | File creation timestamp                           |
| `uuid`             | string           | Unique ID                                         |
| `author`           | string[]         | Author(s) of data                                 |

## Parameters

```
/parameters
    ├── kpoints
    ├── ecutwfc
    ├── basis
    └── other_input
```
Recommended parameters for plane-wave calculations

| Attribute        | Type   | Description              |
| ---------------- | ------ | ------------------------ |
| `ecutwfc`        | float  | Wavefunction cutoff (eV) |
| `ecutrho`        | float  | Density cutoff (eV)      |
| `smearing`       | string | Smearing type            |
| `sigma`          | float  | Smearing width           |


Recommended parameters for Gaussian basis calculations
| Attribute        | Type   | Description              |
| ---------------- | ------ | ------------------------ |
| `basis`        | float  | Wavefunction cutoff (eV) |


Recommended parameters for PBC calculations

| Attribute        | Type   | Description              |
| ---------------- | ------ | ------------------------ |
| `kpoint_grid`    | int[3] | Monkhorst-Pack grid      |

Recommended parameters for all calculations

| Attribute        | Type   | Description              |
| ---------------- | ------ | ------------------------ |
| `spin_polarized` | bool   | Spin calculation         |


Recommended parameters for DMC calculations

| Attribute        | Type   | Description              |
| ---------------- | ------ | ------------------------ |
| `time_step`      | float  | MD timestep (fs)         |
| `n_steps`        | int    | Number of steps          |



## Code (reproducibility metadata)

```
/code
    ├── .attrs
    ├── input_file
    ├── stdout
    └── git
```

.attrs

| Attribute         | Type   | Description             |
| ----------------- | ------ | ----------------------- |
| `name`            | string | Code name (e.g. "QMCPACK", "PYQMC") |
| `version`         | string | Code version            |
| `compiler`        | string | Compiler used           |
| `parallelization` | string | MPI/OpenMP info         |



| Dataset      | Type   | Description     |
| ------------ | ------ | --------------- |
| `input_file` | string | Full input file |
| `stdout`     | string | Output log      |
| `system` | string | HPC system calculation was done on        |
| `location` | string | Storage location on the HPC resource      |


## Structure

```
/structure
    ├── lattice_vectors
    ├── positions
    ├── fractional_positions
    ├── pbc
    └── symmetry
```

Recommended parameters

| Dataset                | Shape     | Description             |
| ---------------------- | --------- | ----------------------- |
| `lattice_vectors`      | (3,3)     | Lattice matrix (Å)      |
| `positions`            | (N,3)     | Cartesian positions (Å) |
| `fractional_positions` | (N,3)     | Fractional coordinates  |
| `pbc`                  | (3,) bool | Periodic boundary flags |
| `symmetry/spacegroup`  | string    | Spacegroup label        |

## Observables

These are outputs of the calculation. Try to use the same names present here. Error bars should be represented by a postfix of `_error`, not `_err` etc.

```
/observables
    ├── energy
    ├── forces
    ├── stress
    └── magnetic_moments
```

Example data sets

| Dataset            | Shape  | Description            |
| ------------------ | ------ | ---------------------- |
| `total_energy`     | (nstates,) | Total electronic energy (eV)      |
| `total_energy_error`     | (nstates,) | Error in the total energy (eV)      |
| `kinetic_energy`     | (nstates,) | Total energy (eV)      |
| `kinetic_energy_error`     | (nstates,) | Error in the total energy (eV)      |
| `forces`           | (nstates, N,3)  | Forces (eV/Å)          |
| `forces_error`           | (nstates, N,3)  | Forces (eV/Å)          |
| `stress`           | (3,3)  | Stress tensor (GPa)    |
| `structure_factor`           | (Nk,3)  | $ S(k) $    |
| `structure_factor_error`           | (Nk,3)  | $ S(k) $    |
| `structure_factor_ks`           | (Nk,3)  | K-points    |
| `magnetic_moments` | (nstates, N,)   | Local magnetic moments |

## Provenance

In this section, enter data that may be necessary to 

```
/provenance
    ├── dependent_uuids
    └── other data
```
| Dataset            | Shape  | Description            |
| ------------------ | ------ | ---------------------- |
| `dependent_uuids`     | string [] | uuids of reference DFT calculations, for example, or  |

---

## Upload, View, and Restore HDF5 Files

This repository includes shared bucket utility scripts at the repository root:

- `upload_h5_to_aws.sh`
- `osn_check_and_restore_h5.sh`

These are intended for all contributors uploading schema-compliant HDF5 files to:

- `s3://phy240060/QMCHAMM/`

### Prerequisites

1. Python with `awscli` module available:
   - scripts call `python3 -m awscli ...`
2. Network access to OSN endpoint:
   - `https://uri.osn.mghpcc.org`
3. OSN credentials:
   - RW key/secret for upload
   - RO key/secret for list/check/restore

### Upload all `.h5` files

Run from repository root:

```bash
cd /projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/database_work/qmchamm_database
bash upload_h5_to_aws.sh
```

The script will prompt:

- `Source directory for .h5 upload [...]`

Press Enter to use the default (repo root), or provide any folder path containing `.h5` files.

Default behavior:

- source directory: script directory
- destination: `s3://phy240060/QMCHAMM/`
- dry-run: enabled by default (`DO_DRYRUN=1`)

Optional environment override:

```bash
SOURCE_DIR=/path/to/h5 \
S3_BUCKET=s3://phy240060/QMCHAMM \
OSN_ENDPOINT=https://uri.osn.mghpcc.org \
DO_DRYRUN=1 \
bash upload_h5_to_aws.sh
```

Set `DO_DRYRUN=0` to skip the dry-run.

### List/check/restore files from bucket

```bash
bash osn_check_and_restore_h5.sh list
bash osn_check_and_restore_h5.sh check
bash osn_check_and_restore_h5.sh restore-missing
bash osn_check_and_restore_h5.sh restore-file <uuid>.h5
```

Optional override example:

```bash
SOURCE_DIR=/path/to/local_h5 \
S3_BUCKET=s3://phy240060/QMCHAMM \
OSN_ENDPOINT=https://uri.osn.mghpcc.org \
bash osn_check_and_restore_h5.sh check
```

### Security notes

- Credentials are requested interactively and are not stored in scripts.
- Credentials are exported only in the current process and unset at script end.
The keys for the RW and RO credentials can be accessed here https://coldfront.osn.mghpcc.org/project/329/