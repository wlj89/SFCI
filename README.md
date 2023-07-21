# SFCI

This repository conatins the codes and the electron-electron integrals required to reproduce the results in _Sparse Full Configuration Interaction_.
!A package for sparse full configuration interaction (SFCI). 

## Usage

1. Enter `src/` and compile with `make`. You might need to specify the flag relating to cpu architecture. You also need to install Eigen and passing its path via `-I` flag.
2. Specify the required parameters in **sfci.input**. One should first specify the task type with `task=ground`. Besides, the following parameters are required: 

    a) `Nd`: max number of determinant contained in each Krylov basis vector
    
    b) `Nl`: dimensionality of the effective Hamiltonian 
    
    c) `num_restart`: number of restarting 
    
    d) `num_orb`: number of spin orbitals 
    
    e) `num_e`: number of electrons
    
    f) `num_thread`: number of cores to be used
    
    g) `filename`: file name of input integral. See below for more details.

    h) `init_hf`: bitstring of the inital HF determinant written in decimal form. 
    
    i) `bucket_unit`(optional): bucket memory allocation unit size. Since there is no obvious way to predict the size of buckets in a SpMSpVM, each bucket's memory is allocated generously. The default value is 60000, which is good for the most cases. But it might still lead to overflow. If this happens, increase bucket_unit for remedy. 

3. run with 

    `./sfci`

## Reformatting Integrals

Both RHF and UHF integrals in FCIDUMP format are accpeted, but requires some reformatting. Use `utils/convert_fcidump.py` to reformat FCIDUMP file:

    python convert_fcidump.py -r sys1-rhf_fcidump sys2-rhf_fcidump ...

or

	python convert_fcidump.py -u sys1-uhf_fcidump sys2-uhf_fcidump ...

The reformated integral files will all end with \_IN, and are now recognizable by the SFCI main program. 
