import h5py
import glob
import numpy as np
import uuid
import datetime
from run_dft import parse_xyz_file
from run_casci import define_system_to_active_space_map
from gather_observables_for_database import (
    gather_casci_data,
    gather_qmc_data,
)


def fill_attrs(
    f,
    system,
    nat,
    calculation_type,
):
    formulas = {
        "FeAlN_32_atoms": "FeAl15N16",
        "CrAlN_32_atoms": "CrAl15N16",
        "NVdiamond_31_atoms": "C30N1",
        "SiVdiamond_31_atoms": "C30Si1",
    }
    species = {
        "FeAlN": ["Fe", "Al", "N"],
        "CrAlN": ["Cr", "Al", "N"],
        "NVdiamond": ["N", "C"],
        "SiVdiamond": ["Si", "C"],
    }
    methods = {
        "CASCI": "CASCI",
        "QMC": "VMC",
    }
    f["system"] = system
    f["formula"] = formulas[f"{system}_{nat}_atoms"]
    f["natoms"] = nat
    f["species"] = species[system]
    f["calculation_type"] = calculation_type
    f["method"] = methods[calculation_type]
    f["author"] = ["Kevin G. Kleiner", "Lucas K. Wagner"]
    f["creation_date"] = datetime.datetime.now().isoformat()
    f["uuid"] = str(uuid.uuid4())


def fill_parameters(f, system, nat, basis, calculation_type):
    charge = {
        "FeAlN": 0,
        "CrAlN": 1,
        "NVdiamond": -1,
        "SiVdiamond": 0,
    }
    sz = {
        "FeAlN": 2.5,
        "CrAlN": 1,
        "NVdiamond": 1,
        "SiVdiamond": 1,
    }
    parameters = f.create_group("parameters")
    dft = parameters.create_group("dft")
    casci = parameters.create_group("casci")
    dft["calculation_type"] = "KROKS"
    dft["sz"] = sz[system]
    dft["charge"] = charge[system]
    dft["xc"] = ["pbe", "pbe0"]
    dft["ecp"] = "ccecp"
    dft["basis"] = f"ccecpccp{basis}"
    dft["exp_to_discard"] = 0.1
    dft["max_cycle"] = 400
    dft["precision"] = 1e-08
    dft["conv_tol"] = 1e-06
    dft["kpoint_grid"] = [1, 1, 1]
    dft["init_guess"] = "atomic"
    cas_mo_indices, nup, ndown = define_system_to_active_space_map()[
        f"{system}_{nat}_atoms"
    ]
    casci["ncas"] = len(cas_mo_indices)
    casci["cas_mo_indices"] = cas_mo_indices - 1
    casci["nelecas"] = [nup, ndown]
    if calculation_type == "QMC":
        ground_state_jastrow_optimization = parameters.create_group(
            "ground_state_jastrow_optimization"
        )
        ensemble_optimization = parameters.create_group("ensemble_optimization")
        vmc = parameters.create_group("vmc")
        ground_state_jastrow_optimization["nconfig"] = 3000
        ground_state_jastrow_optimization["nblocks"] = 10
        ground_state_jastrow_optimization["geminal"] = False
        ground_state_jastrow_optimization["max_iterations"] = 30
        ensemble_optimization["determinants"] = True
        ensemble_optimization["orbitals"] = True
        ensemble_optimization["nconfig"] = 100000
        ensemble_optimization["nblocks"] = 4
        ensemble_optimization["number_to_opt"] = 200
        ensemble_optimization["tau"] = 0.05
        vmc["determinants"] = True
        vmc["orbitals"] = True
        vmc["nconfig"] = 60000
        vmc["nblocks"] = 150
        vmc["include_rdm2"] = False
        vmc["defect_basis_only"] = True


def fill_code(f, calculation_type):
    code = f.create_group("code")
    dft_casci = code.create_group("dft_casci")
    dft_casci["name"] = "pyscf"
    dft_casci["dependencies"] = ["numpy", "scipy", "h5py", "setuptools"]
    dft_casci["system"] = "Illinois Campus Cluster"
    dft_casci["version"] = "2.12.0"
    dft_casci["modules"] = [
        "module load anaconda3",
        "module load openmpi",
        "module load gcc",
    ]
    dft_casci["OMP_NUM_THREADS"] = 8
    dft_casci["MKL_NUM_THREADS"] = 8
    with open("submit_dft.py", "r", encoding="utf-8") as submit_f:
        dft_casci["dft_submit_script"] = submit_f.read()
    with open("run_dft.py", "r", encoding="utf-8") as input_f:
        dft_casci["dft_python_script"] = input_f.read()
    with open(glob.glob("dft_*.out")[0], "r", encoding="utf-8") as output_f:
        dft_casci["dft_stdout"] = output_f.read()
    with open("submit_casci.py", "r", encoding="utf-8") as submit_f:
        dft_casci["casci_submit_script"] = submit_f.read()
    with open("run_casci.py", "r", encoding="utf-8") as input_f:
        dft_casci["casci_python_script"] = input_f.read()
    with open(glob.glob("casci_*.out")[0], "r", encoding="utf-8") as output_f:
        dft_casci["casci_stdout"] = output_f.read()
    if calculation_type == "QMC":
        qmc = code.create_group("qmc")
        ground_state_jastrow_optimization = qmc.create_group(
            "ground_state_jastrow_optimization"
        )
        ensemble_optimization = qmc.create_group("ensemble_optimization")
        vmc = qmc.create_group("vmc")
        qmc["name"] = "pyqmc"
        qmc["dependencies"] = [
            "numpy",
            "scipy",
            "h5py",
            "setuptools",
            "pyscf",
            "pandas",
            "python-dateutil",
            "pytz",
            "tzdata",
            "six",
        ]
        qmc["version"] = "0.6.0"
        qmc["parallelization"] = "mpi4py"
        qmc["OMP_NUM_THREADS"] = 1
        qmc["MKL_NUM_THREADS"] = 1
        qmc["NUMEXPR_NUM_THREADS"] = 1
        ground_state_jastrow_optimization["system"] = "Illinois Campus Cluster"
        ensemble_optimization["system"] = "Aurora"
        vmc["system"] = "Aurora"
        ground_state_jastrow_optimization["modules"] = [
            "module load anaconda3",
            "module load openmpi",
            "module load gcc",
        ]
        ensemble_optimization["modules"] = [
            "module load frameworks",
            "module use /soft/modulefiles",
        ]
        vmc["modules"] = ["module load frameworks", "module use /soft/modulefiles"]
        ground_state_jastrow_optimization["ncores_per_node"] = 128
        ground_state_jastrow_optimization["nnodes"] = 1
        ensemble_optimization["ncores_per_node"] = 102
        ensemble_optimization["nnodes_per_defect"] = 22
        vmc["ncores_per_node"] = 102
        vmc["nnodes_per_state"] = 17
        with open(
            "submit_optimize_jastrow_ground_state.py", "r", encoding="utf-8"
        ) as submit_f:
            ground_state_jastrow_optimization["submit_script"] = submit_f.read()
        with open("optimize_jastrow_ground_state.py", "r", encoding="utf-8") as input_f:
            ground_state_jastrow_optimization["python_script"] = input_f.read()
        with open(glob.glob("jastrow_*.out")[0], "r", encoding="utf-8") as output_f:
            ground_state_jastrow_optimization["stdout"] = output_f.read()
        with open(
            "submit_optimize_ensemble_aurora.py", "r", encoding="utf-8"
        ) as submit_f:
            ensemble_optimization["submit_script"] = submit_f.read()
        with open("optimize_ensemble.py", "r", encoding="utf-8") as input_f:
            ensemble_optimization["python_script"] = input_f.read()
        with open(
            glob.glob("ensemble_opt_*.out")[0], "r", encoding="utf-8"
        ) as output_f:
            ensemble_optimization["stdout"] = output_f.read()
        with open(
            "submit_optimized_wf_vmc_aurora.py", "r", encoding="utf-8"
        ) as submit_f:
            vmc["submit_script"] = submit_f.read()
        with open("run_optimized_wf_vmc.py", "r", encoding="utf-8") as input_f:
            vmc["python_script"] = input_f.read()
        with open(glob.glob("vmc_*.out")[0], "r", encoding="utf-8") as output_f:
            vmc["stdout"] = output_f.read()
        with open("average_qmc_data.py", "r", encoding="utf-8") as average_f:
            code["average_qmc_python_script"] = average_f.read()
        with open("generate_symmops.py", "r", encoding="utf-8") as symmops_f:
            code["generate_symmops_python_script"] = symmops_f.read()


def fill_structure(f, system, nat):
    spacegroups = {
        "FeAlN": "P63mc",
        "CrAlN": "P63mc",
        "NVdiamond": "Fd3m",
        "SiVdiamond": "Fd3m",
    }
    pointgroups = {
        "FeAlN": "c3v",
        "CrAlN": "c3v",
        "NVdiamond": "c3v",
        "SiVdiamond": "d3d",
    }
    lattice_vectors, atom = parse_xyz_file(
        f"geometry_files/{system}/supercell_{nat}_atoms_relaxed.xyz"
    )
    positions = np.array(
        [
            list(map(float, line.split()[1:4]))
            for line in atom.splitlines()
            if line.strip()
        ]
    )
    structure = f.create_group("structure")
    structure["lattice_vectors"] = lattice_vectors
    structure["positions"] = positions
    structure["fractional_positions"] = positions @ np.linalg.inv(lattice_vectors)
    structure["pbc"] = [True, True, True]
    symmetry = structure.create_group("symmetry")
    symmetry["spacegroup"] = spacegroups[system]
    symmetry["pointgroup"] = pointgroups[system]


def fill_observables_casci(f, system, nat, basis):
    observables = f.create_group("observables")
    for xc in ["pbe", "pbe0"]:
        df = gather_casci_data(system, nat, xc, basis)
        wf_group = observables.create_group(xc)
        wf_group["eigenstate"] = df["eigenstate"].values
        wf_group["total_energy"] = df["total_energy"].values
        wf_group["excitation_energy"] = df["excitation_energy"].values
        wf_group["rdm1"] = np.stack(df["rdm1"].values)
        wf_group["spin_squared"] = df["spin_squared"].values
        wf_group["irrep"] = df["irrep"].values


def fill_observables_qmc(f, system, nat, basis):
    observables = f.create_group("observables")
    for xc, wf_type in [
        ["pbe", "casci_j2"],
        ["pbe0", "casci_j2"],
        ["pbe0", "optimized_j2"],
        ["pbe0", "optimized_j2_determinants"],
        ["pbe0", "optimized_j2_determinants_orbitals"],
        ["pbe0", "casci_jg"],
    ]:
        df = gather_qmc_data(system, nat, xc, basis, wf_type=wf_type)
        wf_group = observables.create_group(f"{xc}_{wf_type}")
        wf_group["eigenstate"] = df["eigenstate"].values
        wf_group["total_energy"] = df["total_energy"].values
        wf_group["total_energy_error"] = df["total_energy_error"].values
        wf_group["excitation_energy"] = df["excitation_energy"].values
        wf_group["excitation_energy_error"] = df["excitation_energy_error"].values
        wf_group["rdm1"] = np.stack(df["rdm1"].values)
        wf_group["rdm1_error"] = np.stack(df["rdm1_error"].values)
        wf_group["spin_squared"] = df["spin_squared"].values
        wf_group["irrep"] = df["irrep"].values


def generate_h5(system, nat):
    fill_observables_function = {
        "CASCI": fill_observables_casci,
        "QMC": fill_observables_qmc,
    }
    basis = "vqz"
    for calculation_type in ["CASCI", "QMC"]:
        fname = f"{system}_{nat}_atoms_{calculation_type}.h5"
        print(f"Generating {fname}")
        with h5py.File(fname, "w") as f:
            fill_attrs(f, system, nat, calculation_type)
            fill_code(f, calculation_type)
            fill_structure(f, system, nat)
            fill_parameters(f, system, nat, basis, calculation_type)
            fill_observables_function[calculation_type](f, system, nat, basis)


if __name__ == "__main__":
    for system, nat in [
        ["FeAlN", 32],
        ["CrAlN", 32],
        ["NVdiamond", 31],
        ["SiVdiamond", 31],
    ]:
        generate_h5(system, nat)
