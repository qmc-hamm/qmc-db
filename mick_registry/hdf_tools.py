'''
Used for generating atomic coordinates
'''
import ase
import h5py
import pandas as pd
import numpy as np

def get_basis(d, disregistry, zshift='CM', c=0, a=2.46):

    '''
    `disregistry` is unitless, defined such that the distance to translate from AB to AB again is 1.0.
    In angstrom, this distance corresponds to 3*bond_length = 3/sqrt(3)*lattice_constant = sqrt(3)*lattice_constant
    See Fig. 1 in [paper]
    '''
    disregistry_ang = 3**0.5*a*disregistry
    orig_basis = np.array([
        [0, 0, 0],
        [0, a/3**0.5, 0],
        [0, a/3**0.5 + disregistry_ang, d],
        [a/2, a/(2*3**0.5) + disregistry_ang, d]
        ])

    # for open boundary condition in the z-direction
    # move the first layer to the middle of the cell
    if zshift == 'first_layer':
        z = c/2
    # or move the center of mass to the middle of the cell
    elif zshift == 'CM':
        z = c/2 - d/2
    shift_vector = np.array([0, 0, z])
    shifted_basis = orig_basis + shift_vector
    return shifted_basis.tolist()

def get_lattice_vectors(a=2.46):
    '''
    Defines lattice vectors for bilayer graphene
    '''
    return [
        [a, 0, 0],
        [1/2*a, 1/2*3**0.5*a, 0],
        [0, 0, 0]
        ]

def make_hdf(csv_name='qmc.csv', hdf_name='qmc.hdf', ntiling=3):
    '''
    Creates HDF file with data from the specified `csv_name`
    Coordinates `coords` are generated from the supercell `(ntiling, ntiling, 1)`
    '''
    df = pd.read_csv(csv_name)
    data = {
        'coords': [],
        'latvecs': [],
        'basis': [],
        'atomic_numbers': [],
        'ntiling': []
    }

    for i, row in df.iterrows():
        latvecs = get_lattice_vectors()
        basis = get_basis(row['d'], row['disregistry'])
        atoms = ase.atoms.Atoms('C4', positions=basis, cell=latvecs, pbc=(1, 1, 0))
        atoms = atoms.repeat((ntiling, ntiling, 1))
        coords = atoms.get_positions()
        atomic_numbers = atoms.get_atomic_numbers()
        print(atomic_numbers)

        data['coords'].append(coords)
        data['latvecs'].append(latvecs)
        data['basis'].append(basis)
        data['atomic_numbers'].append(atomic_numbers)
        data['ntiling'].append(ntiling)

    with h5py.File(hdf_name, 'w') as hdf:
        for k, v in data.items():
            v = np.array(v)
            hdf.create_dataset(k, data=v, dtype=v.dtype)
        hdf.create_dataset('energy', data=df['energy'])
        hdf.create_dataset('energy_err', data=df['energy_err'])
        hdf.create_dataset('disregistry', data=df['disregistry'])

def read_hdf(hdf_name):
    with h5py.File(hdf_name) as hdf:
        latvecs = hdf['latvecs'][...]
        basis = hdf['basis'][...]
        coords = hdf['coords'][...]
        atomic_numbers = hdf['atomic_numbers'][...]
        energy = hdf['energy'][...]
        energy_err = hdf['energy_err'][...]

    # test extracting only the first configuration
    print(energy[0], energy_err[0])
    print(atomic_numbers[0])
    print(coords[0])
    atoms = ase.atoms.Atoms(numbers=atomic_numbers[0], positions=coords[0], pbc=(1, 1, 0))
    atoms.write('test_read.xyz')

if __name__ == '__main__':
    make_hdf(hdf_name='qmc.hdf')
    # read_hdf(hdf_name='qmc.hdf')
