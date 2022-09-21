# SFCI

A package for sparse full configuration interaction (SFCI). 

## Usage

1. Compile with `make`. You might need to specify the flag relating to cpu architecture. You also need to install Eigen and passing its path via `-I` flag.
2. Specify the required parameters in **sfci.input**. They are: 

    a) `Nd`: size of each the Krylov vectors
    
    b) `Nl`: dimensionality of the effective Hamiltonian 
    
    c) `num_block`: maximal number of effective Hamiltonian to be diaogonalized if the calculation does not converge 
    
    d) `num_orb`: number of spin orbitals 
    
    e) `nume_e`: number of electrons
    
    f) `num_thread`: number of cores to be used
    
    g) `filename`: file name of input integral. See below for more details. 
    
    h) `bucket_unit`(optional): bucket memory allocation unit size. Since there is no obvious way to predict the size of buckets in a SpMSpVM, each bucket's memory is allocated generously. The default value is 60000, which is good for most cases. But it might still lead to overflow. If this happens, increase bucket_unit for remedy. 

3. run with 

    `./sfci`

## Integral Generation

Only FCIDUMP formatt . This . We will soon fix it.  
