import numpy as np
import pandas as pd
import h5py
import glob
import itertools
import pyscf.fci
import pyqmc.api as pyq
from pyqmc.wf.determinant_tools import reformat_binary_dets


def gather_casci_rdms(mc, casci_index):
    (rdm1_up, rdm1_down), (
        rdm2_upup,
        rdm2_updown,
        rdm2_downdown,
    ) = pyscf.fci.direct_spin1.make_rdm12s(
        fcivec=mc.ci[casci_index],
        norb=mc.ncas,
        nelec=mc.nelecas,
        reorder=False,
    )
    rdm2_downup = rdm2_updown.transpose(2, 3, 0, 1)
    rdm1 = np.sum([rdm1_up, rdm1_down], axis=0)
    rdm2 = np.sum([rdm2_upup, rdm2_updown, rdm2_downdown, rdm2_downup], axis=0)
    return rdm1, rdm2


def select_irrep_c3v(symmetry_expectation_values):
    c3_expectation_value = symmetry_expectation_values["c3"]
    sigma_v_expectation_value = symmetry_expectation_values["sigma_v"]
    spatial_irrep = np.select(
        [
            (c3_expectation_value > 0) & (sigma_v_expectation_value > 0),
            (c3_expectation_value > 0) & (sigma_v_expectation_value < 0),
            c3_expectation_value < 0,
        ],
        [r"A_{1}", r"A_{2}", "E"],
        default=np.array([0], dtype=np.float64),
    )[0]
    return spatial_irrep


def select_irrep_d3d(symmetry_expectation_values):
    c3_expectation_value = symmetry_expectation_values["c3"]
    sigma_d_expectation_value = symmetry_expectation_values["sigma_d"]
    inversion_expectation_value = symmetry_expectation_values["inv"]
    spatial_irrep = np.select(
        [
            (c3_expectation_value > 0) & (sigma_d_expectation_value > 0),
            (c3_expectation_value > 0) & (sigma_d_expectation_value < 0),
            c3_expectation_value < 0,
        ],
        [r"A_{1g}", r"A_{2g}", r"E_{g}"],
        default=np.array([0], dtype=np.float64),
    )[0]
    if inversion_expectation_value < 0:
        spatial_irrep = spatial_irrep.replace("g", "u")
    return spatial_irrep


def generate_symmops_determinant_basis(symmops_orbital_basis, mc):
    symmops_determinant_basis = {}
    for symmop in symmops_orbital_basis:
        symmop_determinant_basis = []
        for spin_channel in [0, 1]:
            symmop_determinant_basis.append(
                compute_symmop_determinant_basis(
                    symmops_orbital_basis[symmop],
                    mc.ncas,
                    dim=int(mc.nelecas[spin_channel]),
                )
            )
        symmops_determinant_basis[symmop] = symmop_determinant_basis
    return symmops_determinant_basis


def reorder_symmop_orbital_basis(symmop_orbital_basis, Di_orbs, Dj_orbs):
    dim = len(Di_orbs)
    Dj_orbs_repeated = np.repeat(Di_orbs, dim)
    Dj_orbs_tiled = np.tile(Dj_orbs, dim)
    symmop_orbital_basis_new = symmop_orbital_basis[
        Dj_orbs_repeated, Dj_orbs_tiled
    ].reshape(dim, dim)
    return symmop_orbital_basis_new


def define_determinant_label(Di_orbs, Dj_orbs):
    key = (
        "".join(str(orb) for orb in Di_orbs)
        + "_"
        + "".join(str(orb) for orb in Dj_orbs)
    )
    return key


def compute_symmop_determinant_basis(symmop_orbital_basis, norb, dim):
    symmop_determinant_basis = {}
    orbital_combinations = list(itertools.combinations(range(norb), dim))
    for Di_orbs in orbital_combinations:
        for Dj_orbs in orbital_combinations:
            determinant_label = define_determinant_label(Di_orbs, Dj_orbs)
            symmop_orbital_basis_new = reorder_symmop_orbital_basis(
                symmop_orbital_basis, Di_orbs, Dj_orbs
            )
            symmop_determinant_basis[determinant_label] = np.linalg.det(
                symmop_orbital_basis_new
            )
    return symmop_determinant_basis


def evaluate_casci_symmetry_expectation_values(
    mc,
    casci_index,
    symmetry_operators_fname,
):
    symmops_orbital_basis = {}
    with h5py.File(symmetry_operators_fname, "r") as f:
        for symmop in f.keys():
            symmops_orbital_basis[symmop] = f[symmop][...]
    symmops_determinant_basis = generate_symmops_determinant_basis(
        symmops_orbital_basis, mc
    )
    D_binary = pyscf.fci.addons.large_ci(
        mc.ci[casci_index], mc.ncas, mc.nelecas, tol=1e-2
    )
    D_list = reformat_binary_dets(D_binary, tol=1e-2)
    symmetry_expectation_values = {}
    for symmop in symmops_determinant_basis:
        symmetry_expectation_value = 0
        for wi, Di_orbs in D_list:
            for wj, Dj_orbs in D_list:
                spin_up_determinant_label = define_determinant_label(
                    Di_orbs[0], Dj_orbs[0]
                )
                spin_down_determinant_label = define_determinant_label(
                    Di_orbs[1], Dj_orbs[1]
                )
                symmetry_expectation_value += (
                    wi
                    * wj
                    * symmops_determinant_basis[symmop][0][spin_up_determinant_label]
                    * symmops_determinant_basis[symmop][1][spin_down_determinant_label]
                )
        symmetry_expectation_values[symmop] = symmetry_expectation_value
    return symmetry_expectation_values


def define_defect_mo_indices_map():
    defect_mo_indices_map = {
        "FeAlN": [68, 69, 70, 71, 72],
        "CrAlN": [68, 69, 70, 71, 72],
        "NVdiamond": [57, 61, 62, 63],
        "SiVdiamond": [44, 48, 58, 59, 61, 62],
    }
    return defect_mo_indices_map


def gather_vmc_rdm1(f_vmc, system):
    defect_mo_indices = define_defect_mo_indices_map()[system]
    norb = f_vmc["rdm1_up_pbe0value"][...].shape[0]
    rdm1 = np.zeros((norb, norb))
    rdm1_error_squared = np.zeros((norb, norb))
    for spin_channel in ["up", "down"]:
        rdm1 += f_vmc[f"rdm1_{spin_channel}_pbe0value"][...]
        rdm1_error_squared += f_vmc[f"rdm1_{spin_channel}_pbe0value_err"][...] ** 2
    rdm1_error = np.sqrt(rdm1_error_squared)
    rdm1 = rdm1[np.ix_(defect_mo_indices, defect_mo_indices)]
    rdm1_error = rdm1_error[np.ix_(defect_mo_indices, defect_mo_indices)]
    return rdm1, rdm1_error


def gather_spin_squared_irrep_from_casci(system, nat, mc, casci_index):
    point_groups = {
        "FeAlN": "c3v",
        "CrAlN": "c3v",
        "NVdiamond": "c3v",
        "SiVdiamond": "d3d",
    }
    irrep_selector = {
        "c3v": select_irrep_c3v,
        "d3d": select_irrep_d3d,
    }
    symmetry_expectation_values = evaluate_casci_symmetry_expectation_values(
        mc,
        casci_index,
        f"symmetry_data/symmops_{system}_{nat}_atoms.hdf5",
    )
    rdm1, rdm2 = gather_casci_rdms(mc, casci_index)
    spin_squared = 0.5 * (mc.ncas + 2) * np.einsum("ii->", rdm1) - 0.25 * (
        2 * np.einsum("ijji->", rdm2) + np.einsum("iijj->", rdm2)
    )
    s = np.sqrt(spin_squared + 0.25) - 0.5
    s = np.round(s * 2) / 2
    spatial_irrep = irrep_selector[point_groups[system]](symmetry_expectation_values)
    irrep = rf"$^{{{int(2 * s + 1)}}}{spatial_irrep}$"
    return spin_squared, irrep


def gather_casci_data(system, nat, xc, basis):
    df = pd.DataFrame({})
    eV_per_har = 27.2114
    mf_chkfile = f"{system}_{nat}_atoms/dft/kroks_{xc}_{basis}_unc_False_cart_False.chk"
    ci_chkfile = (
        f"{system}_{nat}_atoms/casci/casci_{xc}_{basis}_unc_False_cart_False.chk"
    )
    _, _, mc = pyq.recover_pyscf(mf_chkfile, ci_checkfile=ci_chkfile)
    for casci_index in range(mc.energy.shape[0]):
        d = {}
        d["casci_index"] = casci_index
        d["total_energy"] = mc.energy[casci_index] * eV_per_har
        d["excitation_energy"] = (mc.energy[casci_index] - mc.energy[0]) * eV_per_har
        rdm1, _ = gather_casci_rdms(mc, casci_index)
        d["rdm1"] = rdm1
        d["spin_squared"], d["irrep"] = gather_spin_squared_irrep_from_casci(
            system, nat, mc, casci_index
        )
        df = pd.concat([df, pd.DataFrame([d])], sort=False)
    df["eigenstate"] = range(len(df.index))
    return df


def gather_qmc_data(system, nat, xc, basis, wf_type):
    df = pd.DataFrame({})
    data_folders = {
        "casci_j2": f"casci_jastrow_{xc}_{basis}_unc_False_cart_False_averaged",
        "optimized_j2": f"optimized_determinants_False_orbitals_False_{xc}_{basis}_unc_False_cart_False_averaged",
        "optimized_j2_determinants": f"optimized_determinants_True_orbitals_False_{xc}_{basis}_unc_False_cart_False_averaged",
        "optimized_j2_determinants_orbitals": f"optimized_determinants_True_orbitals_True_{xc}_{basis}_unc_False_cart_False_averaged",
        "casci_jg": f"casci_jastrow_{xc}_{basis}_unc_False_cart_False_geminal_ns_4_alpha_s_0.1_averaged",
    }
    data_folder = data_folders[wf_type]
    eV_per_har = 27.2114
    mf_chkfile = f"{system}_{nat}_atoms/dft/kroks_{xc}_{basis}_unc_False_cart_False.chk"
    ci_chkfile = (
        f"{system}_{nat}_atoms/casci/casci_{xc}_{basis}_unc_False_cart_False.chk"
    )
    _, _, mc = pyq.recover_pyscf(mf_chkfile, ci_checkfile=ci_chkfile)
    for vmc_averaged_chkfile in glob.glob(
        f"{system}_{nat}_atoms/qmc/{data_folder}/vmc_*.chk"
    ):
        d = {}
        casci_index = int(
            vmc_averaged_chkfile.split("/")[-1].split("_")[-1].replace(".chk", "")
        )
        d["casci_index"] = casci_index
        with h5py.File(vmc_averaged_chkfile, "r") as f_vmc:
            d["total_energy"] = f_vmc["energytotal"][...] * eV_per_har
            d["total_energy_error"] = f_vmc["energytotal_err"][...] * eV_per_har
            d["excitation_energy"] = f_vmc["excitation_energy"][...] * eV_per_har
            d["excitation_energy_error"] = (
                f_vmc["excitation_energy_err"][...] * eV_per_har
            )
            if "rdm1_up_pbe0value" in f_vmc.keys():
                d["rdm1"], d["rdm1_error"] = gather_vmc_rdm1(f_vmc, system)
        d["spin_squared"], d["irrep"] = gather_spin_squared_irrep_from_casci(
            system, nat, mc, casci_index
        )
        df = pd.concat([df, pd.DataFrame([d])], sort=False)
    df = df.sort_values(by="excitation_energy")
    df["eigenstate"] = range(len(df.index))
    for rdm_entry in ["rdm1", "rdm1_error"]:
        df[rdm_entry] = df[rdm_entry].apply(
            lambda x: (
                x if isinstance(x, np.ndarray) else np.full((mc.ncas, mc.ncas), np.nan)
            )
        )
    return df
