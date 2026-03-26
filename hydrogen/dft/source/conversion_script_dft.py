import os
import re
import json
import glob
import h5py
import numpy as np
from ase.io import read

def parse_filename(filename):
    """
    Extracts pressure, temperature as strings (to preserve padding), 
    and the appropriate output folder from the trajectory filename.
    """
    # Regex to find -pXXX and -tXXX in the filename
    p_match = re.search(r'-p(\d+)', filename)
    t_match = re.search(r'-t(\d+)', filename)
    
    # Keep as strings to preserve formatting like '050'
    p_str = p_match.group(1) if p_match else None
    t_str = t_match.group(1) if t_match else None
    
    # Determine the target folder based on filename clues
    if 'vdw' in filename.lower():
        target_folder = 'vdw'
    elif 'pbe' in filename.lower():
        target_folder = 'pbe'
    else:
        target_folder = 'pbe' # Default fallback
        
    return p_str, t_str, target_folder

def write_dataset_with_attrs(group, name, data, schema_info):
    """
    Helper function to create a dataset and attach description/units as attributes.
    """
    if data is None:
        return
        
    dset = group.create_dataset(name, data=data)
    if 'units' in schema_info:
        dset.attrs['units'] = schema_info['units']
    if 'description' in schema_info:
        dset.attrs['description'] = schema_info['description']

def main():
    # Define paths relative to the script's location
    base_dir = os.path.dirname(os.path.abspath(__file__))
    source_dir = base_dir # Since script is in 'source'
    parent_dir = os.path.dirname(base_dir) # The 'dft' folder
    
    # Find all .traj files in subdirectories
    traj_files = glob.glob(os.path.join(source_dir, '**', '*.traj'), recursive=True)
    
    for traj_path in traj_files:
        filename = os.path.basename(traj_path)
        # Extract the immediate parent folder name (e.g., 'f30c' or 'f36a')
        source_folder = os.path.basename(os.path.dirname(traj_path))
        
        print(f"Processing: {os.path.join(source_folder, filename)}")
        
        # Parse filename for P, T, and output directory
        p_str, t_str, target_folder = parse_filename(filename)
        
        if p_str is None or t_str is None:
            print(f"  Skipping {filename}: Could not parse P and T from filename.")
            continue
            
        # Read the trajectory (all frames)
        traj = read(traj_path, index=':')
        nframes = len(traj)
        
        if nframes == 0:
            print(f"  Skipping {filename}: Trajectory is empty.")
            continue
            
        # Reference first frame for static properties (parameters, code, root attrs)
        frame0 = traj[0]
        calc0 = frame0.calc
        
        # Prepare lists for dynamic data (nframes)
        # Structure
        lattice_vectors = []
        positions = []
        
        # Observables
        potential_energy = []
        forces = []
        stress_voigt = []
        
        # MD State
        md_tot_energy = []
        md_pot_energy = []
        md_kin_energy = []
        md_stress_tensor = []
        md_temperature = []
        
        # Provenance
        iconf_list = []
        
        # Loop over all frames to extract dynamic data
        for atoms in traj:
            info = atoms.info
            calc = atoms.calc
            
            # Structure
            lattice_vectors.append(atoms.get_cell()[:])
            positions.append(atoms.get_positions())
            
            # Observables (DFT)
            if calc and 'energy' in calc.results:
                potential_energy.append(calc.results['energy'])
            if calc and 'forces' in calc.results:
                forces.append(calc.results['forces'])
            if calc and 'stress' in calc.results:
                stress_voigt.append(calc.results['stress'])
                
            # MD State (from info dict)
            md_tot_energy.append(info.get('Etot[eV]', np.nan))
            md_pot_energy.append(info.get('Epot[eV]', np.nan))
            md_kin_energy.append(info.get('Ekin[eV]', np.nan))
            md_temperature.append(info.get('T[K]', np.nan))
            
            # Construct 3x3 MD stress tensor from Voigt components in info
            S = np.zeros((3, 3))
            S[0, 0] = info.get('Pxx[GPa]', np.nan)
            S[1, 1] = info.get('Pyy[GPa]', np.nan)
            S[2, 2] = info.get('Pzz[GPa]', np.nan)
            S[0, 1] = S[1, 0] = info.get('Pxy[GPa]', np.nan)
            S[0, 2] = S[2, 0] = info.get('Pxz[GPa]', np.nan)
            S[1, 2] = S[2, 1] = info.get('Pyz[GPa]', np.nan)
            md_stress_tensor.append(S)
            
            # Provenance
            iconf_list.append(info.get('iconf', -1))

        # Setup Output HDF5 path using the new naming convention
        out_dir = os.path.join(parent_dir, target_folder)
        os.makedirs(out_dir, exist_ok=True)
        
        # e.g., f36a_p050_t600.h5
        out_h5_name = f"{source_folder}_p{p_str}_t{t_str}.h5"
        out_path = os.path.join(out_dir, out_h5_name)
        
        # Write to HDF5
        with h5py.File(out_path, 'w') as f:
            
            # ================= 1. Root Attributes (.attrs) =================
            f.attrs['system'] = frame0.get_chemical_formula() + f" P={float(p_str)} T={float(t_str)}"
            f.attrs['formula'] = frame0.get_chemical_formula()
            f.attrs['natoms'] = len(frame0)
            
            # Save string arrays using h5py variable-length strings
            species = list(set(frame0.get_chemical_symbols()))
            f.attrs.create('species', data=np.array(species, dtype=object), dtype=h5py.string_dtype(encoding='utf-8'))
            
            f.attrs['atomic_numbers'] = frame0.get_atomic_numbers()
            
            # Store numerical values in attributes for DB querying
            f.attrs['temperature'] = float(t_str)
            f.attrs['pressure'] = float(p_str)
            f.attrs['uuid'] = frame0.info.get('uuid', 'unknown')

            # ================= 2. Structure =================
            grp_struct = f.create_group('structure')
            write_dataset_with_attrs(grp_struct, 'lattice_vectors', np.array(lattice_vectors), {'units': 'Å', 'description': 'Lattice matrix'})
            write_dataset_with_attrs(grp_struct, 'positions', np.array(positions), {'units': 'Å', 'description': 'Cartesian positions of atoms'})
            write_dataset_with_attrs(grp_struct, 'pbc', frame0.get_pbc(), {'description': 'Periodic boundary condition flags'})

            # ================= 3. Observables =================
            grp_obs = f.create_group('observables')
            if potential_energy:
                write_dataset_with_attrs(grp_obs, 'potential_energy', np.array(potential_energy), {'units': 'eV', 'description': 'Potential energy of the configurations'})
            if forces:
                write_dataset_with_attrs(grp_obs, 'forces', np.array(forces), {'units': 'eV/Å', 'description': 'Forces on atoms'})
            if stress_voigt:
                write_dataset_with_attrs(grp_obs, 'stress', np.array(stress_voigt), {'units': 'eV/Å^3', 'description': 'Stress (virial) of the configurations'})

            # ================= 4. Parameters =================
            grp_param = f.create_group('parameters')
            
            if calc0 and hasattr(calc0, 'parameters'):
                params = calc0.parameters
                inp_data = params.get('input_data', {})
                
                if 'kpts' in params:
                    write_dataset_with_attrs(grp_param, 'kpoint_grid', np.array(params['kpts']), {'description': 'K-point sampling grid'})
                if 'koffset' in params:
                    write_dataset_with_attrs(grp_param, 'koffset', np.array(params['koffset']), {'description': 'K-point grid offset'})
                
                write_dataset_with_attrs(grp_param, 'ecutwfc', inp_data.get('ecutwfc', np.nan), {'units': 'Ry', 'description': 'Wavefunction cutoff energy'})
                write_dataset_with_attrs(grp_param, 'ecutrho', inp_data.get('ecutrho', np.nan), {'units': 'Ry', 'description': 'Density cutoff energy'})
                
                if 'smearing' in inp_data:
                    grp_param.create_dataset('smearing', data=inp_data['smearing'])
                write_dataset_with_attrs(grp_param, 'sigma', inp_data.get('degauss', np.nan), {'units': 'Ry', 'description': 'Smearing width'})
                
                if 'input_dft' in inp_data:
                    grp_param.create_dataset('xc_functional', data=inp_data['input_dft'])
                
                grp_param.create_dataset('vdw', data=target_folder.upper())
                
                if 'pseudopotentials' in params:
                    pseudo_str = json.dumps(params['pseudopotentials'])
                    grp_param.create_dataset('pseudopotentials', data=pseudo_str)

            # ================= 5. Code =================
            grp_code = f.create_group('code')
            code_name = calc0.name if calc0 else "unknown"
            grp_code.create_dataset('name', data=code_name)

            # ================= 6. MD_state =================
            grp_md = f.create_group('MD_state')
            write_dataset_with_attrs(grp_md, 'total_energy', np.array(md_tot_energy), {'units': 'eV', 'description': 'Total energy during the MD simulation'})
            write_dataset_with_attrs(grp_md, 'potential_energy', np.array(md_pot_energy), {'units': 'eV', 'description': 'Potential energy during the MD simulation'})
            write_dataset_with_attrs(grp_md, 'kinetic_energy', np.array(md_kin_energy), {'units': 'eV', 'description': 'Kinetic energy during the MD simulation'})
            write_dataset_with_attrs(grp_md, 'stress', np.array(md_stress_tensor), {'units': 'GPa', 'description': 'Stress tensor components during MD'})
            write_dataset_with_attrs(grp_md, 'temperature', np.array(md_temperature), {'units': 'K', 'description': 'Temperature during the MD simulation step'})

            # ================= 7. Provenance =================
            grp_prov = f.create_group('provenance')
            write_dataset_with_attrs(grp_prov, 'iconf', np.array(iconf_list), {'description': 'Configuration index from trajectory'})

        print(f"  Successfully wrote -> {out_path}")

if __name__ == "__main__":
    main()
