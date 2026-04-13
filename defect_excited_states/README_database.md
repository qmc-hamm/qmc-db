# QMC database for defects

This database contains CASCI and QMC excited state data for defects, including nitrogen-vacancy and silicon-vacancy centers in diamond as well as iron and chromium impurities in aluminum nitride
Each HDF5 file corresponds to a particular defect supercell calculation and contains structure, parameters, code, and observables.

## File Naming

Each file is named:

`<system>_<natoms>_atoms_<calculation_type>.h5`

## HDF5 Layout

`/system/` : defect system name
`/species/` : atomic species in the defect system
`/formula/` : chemical formula with charge for the defect system
`/natoms/` : number of atoms in the supercell for the defect system
`/calculation_type/` : type of calculation done for the defect ('CASCI' or 'QMC')

`/observables/<wf_type>/` : computed ground and excited state observables for each wave function type.
    - `/observables/<wf_type>/eigenstate` : eigenstate indices
    - `/observables/<wf_type>/total_energy` : eigenstate total energy means
    - `/observables/<wf_type>/total_energy_error` : eigenstate total energy standard errors
    - `/observables/<wf_type>/excitation_energy` : eigenstate excitation energy means
    - `/observables/<wf_type>/excitation_energy_error` : eigenstate excitation energy standard errors
    - `/observables/<wf_type>/rdm1` : one-body reduced density matrix means in the defect orbital subspace
    - `/observables/<wf_type>/rdm1_error` : one-body reduced density matrix standard errors in the defect orbital subspace
    - `/observables/<wf_type>/spin_squared` : spin squared expectation values
    - `/observables/<wf_type>/irrep` : spin-spatial irreducible representations for the states

`/structure/` : defect supercell and atomic structure input data.
    - `lattice_vectors` : supercell lattice vectors
    - `positions` : atomic positions within the supercell in Cartesian coordinates (A)
    - `fractional_positions` : atomic positions within the supercell in fractional coordinates
    - `pbc` : periodic boundary condition settings in each dimension
    - `symmetry/pointgroup` : point group for the defect
    - `symmetry/spacegroup` : space group for the host crystal

`/parameters/<calculation stage>/` : input parameters for each stage of the excited state calculation workflow.
    - `/parameters/dft/basis` : basis set used in the DFT calculation of orbitals
    - `/parameters/dft/ecp` : effective core potentials used in DFT
    - `/parameters/dft/calculation_type` : type of DFT calculation (KROKS)
    - `/parameters/dft/charge` : defect charge used in DFT
    - `/parameters/dft/conv_tol` : SCF convergence tolerance used in DFT
    - `/parameters/dft/kpoint_grid` : k-point grid resolution used in DFT
    - `/parameters/dft/max_cycle` : maximum number of SCF cycles allowed in DFT
    - `/parameters/dft/precision` : precision required for the integral evaluations in DFT
    - `/parameters/dft/sz` : spin-z for the target ground state in DFT
    - `/parameters/dft/xc` : exchange-correlation functionals used in DFT
    - `/parameters/dft/exp_to_discard` : exponent cutoff for diffuse functions in the basis set
    - `/parameters/casci/cas_mo_indices` : DFT MO indices to use for active space in CASCI calculation of determinant coefficients
    - `/parameters/casci/ncas` : number of orbitals in the active space in CASCI
    - `/parameters/casci/nelecas` : number of spin-up and spin-down electrons in the active space in CASCI
    - `/parameters/ground_state_jastrow_optimization/nconfig` : number of Monte Carlo walkers used in the ground state Jastrow optimization
    - `/parameters/ground_state_jastrow_optimization/nblocks` : number of blocks/walker used in the Jastrow optimization
    - `/parameters/ground_state_jastrow_optimization/max_iterations` : maximum number of iterations in the Jastrow optimization
    - `/parameters/ground_state_jastrow_optimization/geminal` :  specifies whether a Geminal Jastrow was included in the Jastrow optimization
    - `/parameters/ensemble_optimization/nconfig` : number of Monte Carlo walkers used in the ensemble optimization
    - `/parameters/ensemble_optimization/nblocks` : number of blocks/walker used used in the ensemble optimization
    - `/parameters/ensemble_optimization/determinants` : specifies whether determinant coefficients were included in the ensemble optimization
    - `/parameters/ensemble_optimization/orbitals` : specifies whether orbitals were included in the ensemble optimization
    - `/parameters/ensemble_optimization/number_to_opt` : number of orbital parameters per spin-orbital included in the ensemble optimization
    - `/parameters/ensemble_optimization/tau` : step size for parameter updates in the ensemble optimization
    - `/parameters/vmc/nconfig` : number of Monte Carlo walkers used in the VMC evaluations of energies and density matrices
    - `/parameters/vmc/nblocks` : number of blocks/walker used used in VMC
    - `/parameters/vmc/determinants` : specifies whether to use trial wave functions with optimized determinant coefficients in VMC
    - `/parameters/ensemble_optimization/orbitals` : specifies whether to use trial wave functions with optimized orbitals in VMC
    - `/parameters/ensemble_optimization/defect_basis_only` : specifies whether to restrict the density matrix calculations to the defect orbital subspace
    - `/parameters/ensemble_optimization/include_rdm2` : specifies whether the two-body density matrices were included in the VMC  
        
`/code/` : scripts and outputs in the calculation workflow    
        - `/code/generate_symmops_python_script` : Python script used to generate symmetry representations in the defect orbitals
        - `/code/average_qmc_python_script` : Python script used to average the QMC data
        - `/code/gather_observables_python_script` : Python code containing functionals used to gather the QMC observables
        - `/code/write_h5_python_script` : Python script used to generate the HDF5 file
        - `/code/dft_casci/modules` : list of modules used for DFT and CASCI
        - `/code/dft_casci/system` : computing system used for DFT and CASCI
        - `/code/dft_casci/name` : code name used for DFT and CASCI
        - `/code/dft_casci/version` : code version used for DFT and CASCI
        - `/code/dft_casci/dependencies` : code dependencies required for DFT and CASCI
        - `/code/dft_casci/MKL_NUM_THREADS` : number of threads to use with MKL code in pyscf
        - `/code/dft_casci/OMP_NUM_THREADS` : number of threads to use with OpenMP code in pyscf
        - `/code/dft_casci/dft_submit_script` : submit script used for DFT
        - `/code/dft_casci/dft_python_script` : samplePython script used for DFT
        - `/code/dft_casci/dft_stdout` : sample standard output from DFT
        - `/code/dft_casci/casci_submit_script` : submit script used for CASCI
        - `/code/dft_casci/casci_python_script` : sample Python script used for CASCI
        - `/code/dft_casci/casci_stdout` : sample standard output from CASCI
        - `/code/qmc/name` : code name used for QMC
        - `/code/qmc/parallelization` : parallelization library used for QMC
        - `/code/qmc/version` : pyqmc version used for QMC
        - `/code/qmc/MKL_NUM_THREADS` : number of threads to use with MKL code in pyqmc
        - `/code/qmc/OMP_NUM_THREADS` : number of threads to use with OpenMP code in pyqmc  
        - `/code/qmc/ground_state_jastrow_optimization/modules` : list of modules used for Jastrow optimization
        - `/code/qmc/ground_state_jastrow_optimization/nnodes` : number of compute nodes used for Jastrow optimization
        - `/code/qmc/ground_state_jastrow_optimization/ncores_per_node` : number of cores per nodes used for Jastrow optimization
        - `/code/qmc/ground_state_jastrow_optimization/system` : computing system used for Jastrow optimization
        - `/code/qmc/ground_state_jastrow_optimization/submit_script` : submit script used for Jastrow optimization
        - `/code/qmc/ground_state_jastrow_optimization/python_script` : samplePython script used for Jastrow optimization
        - `/code/qmc/ground_state_jastrow_optimization/stdout` : sample standard output from Jastrow optimization
        - `/code/qmc/ensemble_optimization/modules` : list of modules used for ensemble optimization
        - `/code/qmc/ensemble_optimization/nnodes_per_defect` : number of compute nodes used per defect for ensemble optimization
        - `/code/qmc/ensemble_optimization/ncores_per_node` : number of cores per nodes used for ensemble optimization
        - `/code/qmc/ensemble_optimization/system` : computing system used for ensemble optimization
        - `/code/qmc/ensemble_optimization/submit_script` : submit script used for ensemble optimization
        - `/code/qmc/ensemble_optimization/python_script` : samplePython script used for ensemble optimization
        - `/code/qmc/ensemble_optimization/stdout` : sample standard output from ensemble optimization
        - `/code/qmc/vmc/modules` : list of modules used for VMC
        - `/code/qmc/vmc/nnodes_per_state` : number of compute nodes used per state for VMC
        - `/code/qmc/vmc/ncores_per_node` : number of cores per nodes used for VMC
        - `/code/qmc/vmc/system` : computing system used for VMC
        - `/code/qmc/vmc/submit_script` : submit script used for VMC
        - `/code/qmc/vmc/python_script` : samplePython script used for VMC
        - `/code/qmc/vmc/stdout` : sample standard output from VMC


