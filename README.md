

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
