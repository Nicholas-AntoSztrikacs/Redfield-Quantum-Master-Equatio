# Redfield-Quantum-Master-Equation
This repository contains a set of MATLAB scripts to solve the Redfield Quantum Master Equation for open quantum systems. The implementation allows for the treatment of strong coupling and non-Markovian effects through the use of the reaction coordinate method. The code is flexible and allows for the implementation of generic classes of Hamiltonians, and flexible spectral density functions. Furthermore, the calculation of equilibrium and non-equilibrium steady state observables such as heat currents, reduced density matrix elements, and non-Markovianity wittnesses is possible. 

## Getting Started
If you would like to use this repository for simulating open quantum systems, or would like to get inspired by this repository to make your own, you can either clone or download it

```
git clone https://github.com/Nicholas-AntoSztrikacs/Redfield-Quantum-Master-Equation
```

The main file to run to simulate the dynamics is 'redf_DRC.m'. Cusomizing the Hamiltonians is done via 'Hamiltonian_DRC.m' where any Hamiltonian with a discrete spectrum may be contructed. Furthermore, the initialization of the spectral density functions and the computation of the Redfield tensor is done via 'tensor2_DRC.m'. Lastly, the computation of observables is done in either of the scripts 'SS_current_DRC.m', 'measures.m'. Furthermore, "Ptrace.m" computes the partial trace of the extended system to allow for the computation of the reduced density matrix. Lastly, the scripts "BLP.m" and "RHP.m" are scripts used to compute two measures of non-Markovianity for the exactly solvable pure dephasing model, where "Decoherence.m" produces the exact dynamics for the model. 
